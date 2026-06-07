function [fval, estimates_out, x_upd, x_std, macro_fit, yfit_n, yfit_r, surv_infexp_fit, surv_gdpexp_fit, surv_tbexp_fit, term_prem_n, term_prem_r, irp, fval_t, stderr_all, Hess] = ...
    GGNR_model(macro, yields_n, mats_n, yields_r, mats_r, surv_infexp, hstep_s, ...
    surv_gdpexp, hstep_g, surv_tbexp, hstep_t, tT,  x_0, par_ind_full, eval, max_estimation_iter)
% Code with 3-month T-bill expectations fitted by mid-year values
% No autocorrelated monetary policy shock in z_t (u_t) but
% Time-varying kappa


warning off all

r_lb = []; A_X_for = []; B_X_for = []; s_n = []; Sigma2_X = []; Gamma_m0 = []; Gamma_m = [];
A_X_exp = []; B_X_cum = []; A_X_for_pi = []; B_X_pi = []; A_X_exp_pi = []; B_X_cum_pi = []; Conv_pi = [];
A_X_exp_pi_P = []; B_X_pi_P = []; Conv_pi_P = []; Gamma_s0 = []; Gamma_s = []; ind_r = [];
Gamma_g0 = []; Gamma_g = []; Mu = []; Phi = []; output_message = []; yfit_n = []; yfit_r = []; yfit_n_exp = []; yfit_r_exp = [];
J_smooth = []; x_pred = []; A_X_for_P = []; B_X_for_P = []; A_X_exp_P = []; A_X_for_pi_P = []; B_X_cum_P = []; B_X_cum_pi_P = [];
Phi_j = []; I_Phi_inv = []; delta_1 = []; delta_0 = 0;

if nargin < 15 || isempty(eval)
    eval = 'estimate';
end
if nargin < 16 || isempty(max_estimation_iter)
    max_estimation_iter = length(par_ind_full);
end

if nargout>=15
    stderr_calc = 1;
else
    stderr_calc = 0;
end
verbose_fit = nargout < 14;
%stderr_calc = 0;

if ismember(26,par_ind_full)
    ind_sig_w = 1;
else
    ind_sig_w = 0;
end

% Indicator signalling problems with sigmas (e.g. non-positive or complex):
ind_error_sigma = 0;

freq = 1; % frequency of observations
maxmat=max([mats_n mats_r]*freq);
ind_r = tT<datetime('2004-01-01') | (tT>=datetime('2008-01-01') & tT<datetime('2010-01-01')) ...
     | (tT>=datetime('2020-03-01') & tT<datetime('2021-03-01')); %higher liquidity index


K = 13; %number of primary state variables: [s_t s_{t-1} kappa g z_t z_{t-1} pi_st r_st pi_bar eps_{u,t} eps_{u,t-1} eps_{u,t-2} w]'
[T, N_m] = size(macro);
[T_n, N_n] = size(yields_n);
[T_r, N_r] = size(yields_r);
[T_s, N_s] = size(surv_infexp);
[T_g, N_g] = size(surv_gdpexp);
[T_t, N_t] = size(surv_tbexp);
y_obs_all = [macro surv_infexp surv_gdpexp surv_tbexp yields_n yields_r];


% Bounds on parameter values
x_bounds = x_bounds_vals;

% Transformed parameters by sigmoid function:
x_1 = pars_trans_in(x_0, x_bounds); 



if (T~=T_n) || (T~=T_r) || (T~=T_s) || (T~=T_g) || (~isempty(surv_tbexp) && T~=T_t)
    error('The number of observations must be the same for each series')
end


maxlik_full = @(x_1,startpoint)sum(maxlik_full_t(x_1,startpoint));


allpars(1:length(x_1)) = x_1;


if length(par_ind_full)>length(x_0)
    error('The number of provided starting values is smaller than the number of parameters to estimate')
end

tic

startpoint = x_1(par_ind_full);
fval_t=-maxlik_full_t(startpoint,allpars); %leave it, otherwise est_X is not read (it's needed in estimation
fval = sum(fval_t);
estimates = startpoint;
fval0 = fval;

if verbose_fit
    fprintf('Initial likelihood :       %8.6f;\n', fval0);
end

disp_opt = 'iter-detailed'; % 'off'; %
options_fminsearch = optimset('Display',disp_opt, 'MaxIter', max_estimation_iter);
options_fminunc = optimset('Display',disp_opt);
if strcmp(eval,'eval')==0
    [estimates]=estimate(estimates);
    fval = -maxlik_full(estimates,allpars);
end

estimates_out = x_1;
estimates_out(par_ind_full) = estimates;
% Converting parameters to raw values (deltas for standard errors):
[estimates_out, deltas] = pars_trans_out(estimates_out, x_bounds);



if stderr_calc
    startpoint = x_1(par_ind_full);
    if any(par_ind_full==4) % excluding sigma_i from estimated value
        par_ind_full(par_ind_full==4) = [];
        startpoint(par_ind_full==4) = [];
    end
	fprintf('Computing Hessian-based standard errors at likelihood point %8.6f...\n', -maxlik_full(startpoint, allpars))
    options1=optimset('Display','off','MaxFunEval',0);
    [estimates,fval,~,~,~,Hess]=fminunc(@(startpoint)maxlik_full(startpoint, allpars), startpoint, options1);
    fprintf('Hessian calculation complete.\n')
    VarH = inv(Hess);
    stderr = sqrt(abs(diag(diag(deltas(par_ind_full))'*VarH*diag(deltas(par_ind_full)))))';

    % The fminunc Hessian is used here for the reported standard errors.

    stderr_all = nan(size(x_0));
    stderr_all(par_ind_full) = stderr;
end

%
if ~isempty(output_message)
    fprintf(2,output_message);
end

if verbose_fit
    toc
end

if strcmp(eval,'eval') || isfile('stop.txt') 
    print_fit();
end

term_prem_n = yfit_n - yfit_n_exp;
term_prem_r = yfit_r - yfit_r_exp;




%% Log-likelihood function (negative)
function loglik_t=maxlik_full_t(pars,v)


v(par_ind_full)=pars;

v(6) = pars_trans_in(nanmean(macro(:,1)), x_bounds(6,:)); % mu_g = mean(d_gdp);
x_1(6) = v(6);
v(12) = pars_trans_in(nanmean(macro(:,2)), x_bounds(12,:)); % mu_st = mean(inf);
x_1(12) = v(12);

[A, Mu_til, Phi_til, Sigma_til, Lambda_0, Lambda_1, pars] = VAR_coeffs(v,ind_sig_w, x_bounds);

x_1(16) = pars_trans_in(pars.sigma_kappa, x_bounds(16,:));
x_1(26) = pars_trans_in(pars.sigma_w, x_bounds(26,:));
% To prevent the code from crashing and yet restricting the values:
if ~isreal(pars.sigma_kappa) %|| ~isreal(pars.sigma_w)
    pars.sigma_kappa = real(pars.sigma_kappa);
    pars.sigma_w = real(pars.sigma_w);
    ind_error_sigma = 1;
else
    ind_error_sigma = 0;
end

r_lb = pars.r_lb;


% Measurement errors:
pars.sigma_gdp = pars.sigma_o/100;
x_1(31) = pars_trans_in(pars.sigma_gdp, x_bounds(31,:));
pars.sigma_inf = pars.sigma_o/100;
x_1(32) = pars_trans_in(pars.sigma_inf, x_bounds(32,:));
pars.sigma_ptr = 5*pars.sigma_s;
x_1(33) = pars_trans_in(pars.sigma_ptr, x_bounds(33,:));

sig_ind = 0.25; % relation of sigma_tb to sigma_s

pars.sigma_y = 5*pars.sigma_s; % GDP growth
x_1(35) = pars_trans_in(pars.sigma_y, x_bounds(35,:));
pars.sigma_t = sig_ind*pars.sigma_s; % 3Tbill
x_1(36) = pars_trans_in(pars.sigma_t, x_bounds(36,:));



Mu = A\Mu_til;
Phi = A\Phi_til;

I_Phi_inv = eye(K)/(eye(K) - Phi);

% Calculating Phi^j to avoid repetitions
%hstep_all = unique([hstep_g hstep_t]);
hstep_max = max([hstep_s hstep_g hstep_t]);
Phi_j = zeros(K,K,hstep_max);
Phi_j(:,:,1) = Phi;
Phi_j_g = zeros(K,K,N_g);
Phi_j_t = zeros(K,K,N_t);
k_g = 1;
k_t = 1;
for j = 2:hstep_max
    Phi_j(:,:,j) = Phi_j(:,:,j-1)*Phi;
    if ismember(j, hstep_g)
        Phi_j_g(:,:,k_g) = Phi_j(:,:,j);
        k_g = k_g + 1;
    end
    if ismember(j, hstep_t)
        Phi_j_t(:,:,k_t) = Phi_j(:,:,j);
        k_t = k_t + 1;
    end
end


% forecast variance decomposition for pi_bar - contribution of eps_i
% see Pesaran (24.20)
h_fwd = 12; %forecast horizon for variance error decomposition
theta_pi_i = 0.05; % required forecast variance contribution at the h_fwd horizon
ii = 9; % pi variable's variance decomposition
jj = 1; % i shock

Phi_j_trans = permute(Phi_j(:,:,1:h_fwd-1),[2,1,3]);

Sigma_X_ii = A\Sigma_til(:,2:end);
Sigma2_X_ii = Sigma_X_ii*Sigma_X_ii';
Phi_Sigma_ii_Phi = sum(pagemtimes(pagemtimes(Phi_j(ii,:,1:h_fwd-1),Sigma2_X_ii), Phi_j_trans(:,ii,1:h_fwd-1)));

pars.sigma_i = sqrt(theta_pi_i*(Sigma2_X_ii(ii,ii) + Phi_Sigma_ii_Phi) / ( sum([Phi_j(ii,1,1:h_fwd-1)].^2) - theta_pi_i*sum(Phi_j(9,1,1:h_fwd-1).^2) ));

x_1(4) = pars_trans_in(pars.sigma_i, x_bounds(4,:));
% To prevent the code from crashing when there is no real positive solution:
if pars.sigma_i<=0 || ~isreal(pars.sigma_i)
    pars.sigma_i = real(pars.sigma_i);
    ind_error_sigma = 1;
    output_message = 'sigma_i not real (no solution for 5% forecast error variance) \n';
end
Sigma_til(1,1) = pars.sigma_i;
Sigma_X = A\Sigma_til;
Sigma2_X = Sigma_X*Sigma_X';


% Risk-neutral dynamics
Mu_q = Mu + Sigma_X*Lambda_0;
Phi_q = Phi + Sigma_X*Lambda_1;


delta_1 = [1 zeros(1,K-1)]';
delta_0 = 0;
delta_pi_1 = [zeros(1,8) 1 pars.a_u' 0]';
delta_pi_0 = 0;

 % for inflation expectations:
[A_X_for_P, B_X_for_P, A_X_exp_P, B_X_cum_P, A_X_for_pi_P, B_X_pi_P, A_X_exp_pi_P, B_X_cum_pi_P, Conv_pi_P] = ...
            affine_coefs(Phi, Mu, Sigma2_X, delta_0, delta_1, delta_pi_0, delta_pi_1);

% for fitting yields:
[A_X_for, B_X_for, A_X_exp, B_X_cum, A_X_for_pi, B_X_pi, A_X_exp_pi, B_X_cum_pi, Conv_pi] = ...
            affine_coefs(Phi_q, Mu_q, Sigma2_X, delta_0, delta_1, delta_pi_0, delta_pi_1);


s2_n = cumsum(diag(B_X_for'*(Sigma2_X)*B_X_for)); %
s_n = sqrt(s2_n);



% Kalman terms:
Gamma_m = zeros(N_m,K); % measurement loadings for macro variables
Gamma_m(1,[4:6]) = [1 1 -1]; % Observable: GDP growth
Gamma_m(2,9:9+length(pars.a_u)) = [1 pars.a_u']; % Observable: CPI inflation
Gamma_m(3,5) = 1; % Observable: Output gap
Gamma_m(4,7) = 1; % Perceived Target Rate
Gamma_m0 = zeros(N_m,1);

Gamma_s = zeros(N_s, K);
Gamma_s0 = zeros(N_s,1);
for j = 1:N_s
    Gamma_s(j,:) = sum(B_X_pi_P(:,1:hstep_s(j)),2)'/hstep_s(j); % Observable: expected inflation
    Gamma_s0(j) = sum(A_X_exp_pi_P(1:hstep_s(j))/hstep_s(j) + 0.5*Conv_pi_P(1:hstep_s(j))/hstep_s(j)^2); % Intercept 
end
Gamma_g = nan(N_g, K);
Gamma_g0 = zeros(N_g,1);
for j = 1:N_g
    Gamma_g(j,:) = Gamma_m(1,:)*Phi*I_Phi_inv*(eye(K) - Phi_j_g(:,:,j))/hstep_g(j); % Observable: expected GDP growth
    Gamma_g0(j) = Gamma_m(1,:)*I_Phi_inv*(hstep_g(j)*eye(K) - Phi*I_Phi_inv*(eye(K) - Phi_j_g(:,:,j)))*Mu/hstep_g(j); % Intercept
end



% Time-dependent measurement error for real yields:
sigma_r_t = ones(T,1)*pars.sigma_r;
sigma_r_t(ind_r) = pars.sigma_r + pars.sigma_r_liq;
% Measurement errors:
R_m = [pars.sigma_gdp^2; pars.sigma_inf^2; pars.sigma_o^2; pars.sigma_ptr^2];  % Macro (GDP, inflation, output gap, PTR)
R_s = ones(N_s,1)*pars.sigma_s^2; % Inflation surveys
R_g = ones(N_g,1)*pars.sigma_y^2; % GDP surveys
R_t = ones(N_t,1)*pars.sigma_t^2; % T-bill surveys
R_n = ones(N_n,1)*pars.sigma_n^2; % Nominal yields
R_r = ones(N_r,1)*sigma_r_t(1)^2; % Real yields
R = diag([R_m; R_s; R_g; R_t; R_n; R_r]);

Q = Sigma2_X;
P_pred = reshape((eye(K*K) - kron(Phi,Phi))\Q(:),K,K); % P_{1|0}
P_pred = Phi*P_pred*Phi' + Q;


x_upd = zeros(T,K);
x_pred = zeros(K,T);
%x_pred(:,1) = Mu_unc;% = (eye(K) - Phi)\Mu; % x_{1|0}
x_pred(:,1) = (eye(K) - Phi)\Mu; % x_{1|0}
x_pred(:,1) = Mu + Phi*x_pred(:,1);
x_std = zeros(T,K);
loglik_t = zeros(T,1);
Z_all = eye(N_m + N_s + N_g + N_t + N_n + N_r); % a basis for selection non-missing observations


% Real and Nominal yields
[yfit_all_n, yfit_all_r, JJ_n, JJ_r] = ...
    y_fitting(x_pred(:,1), A_X_for, B_X_for, A_X_exp, r_lb, s_n, A_X_for_pi, B_X_pi, Sigma2_X, B_X_cum, B_X_cum_pi);
y_pred_n = yfit_all_n(mats_n*freq); % predicted nominal yields
y_pred_r = yfit_all_r(mats_r*freq); % predicted real yields
Gamma_n = JJ_n(mats_n*freq,:);
Gamma_r = JJ_r(mats_r*freq,:);


% Expected 3-month T-bill rates
y_pred_t = nan(1,N_t);
Gamma_t = zeros(N_t, K);
obs_t = 1;
if ~isnan(surv_tbexp(obs_t,:))
    yfit_tb_all = zeros(max(hstep_t(~isnan(surv_tbexp(obs_t,:)))),1);
    %yfit_tb_all(1) = yfit_all_n(3);
    Gamma_tb_all = zeros(max(hstep_t(~isnan(surv_tbexp(obs_t,:)))), K);
    %Gamma_tb_all(1,:) = JJ_n(3,:);
    x_pred_h = x_pred(:,obs_t);
    for j = 1:max(hstep_t(~isnan(surv_tbexp(obs_t,:)))-1)
        x_pred_h = Mu + Phi*x_pred_h;
        if mod(j-5,12)==0
            [yfit_all_t, ~, JJ_t] = ...
                y_fitting(x_pred_h, A_X_for(1:3), B_X_for(:,1:3), A_X_exp(1:3), r_lb, s_n(1:3), A_X_for_pi(1:3), B_X_pi(:,1:3), Sigma2_X, B_X_cum(:,1:3), B_X_cum_pi(:,1:3));
            yfit_tb_all(j+1,1) = yfit_all_t(3); % predicted 3-month T-bill
            Gamma_tb_all(j+1,:) = JJ_t(3,:)*Phi_j(:,:,j);
        end
    end
    for j = hstep_t(~isnan(surv_tbexp(obs_t,:)))
        y_pred_t(hstep_t == j) = mean(yfit_tb_all(1:j),'omitnan');
        Gamma_t(hstep_t == j,:) = mean(Gamma_tb_all(1:j,:),'omitnan');
    end
end



Gamma_full = [Gamma_m; Gamma_s; Gamma_g; Gamma_t; Gamma_n; Gamma_r];

y_obs = y_obs_all(1,:)';
y_pred = [(Gamma_m0 + Gamma_m*x_pred(:,1)); (Gamma_s0 + Gamma_s*x_pred(:,1)); ...
    (Gamma_g0 + Gamma_g*x_pred(:,1)); y_pred_t'; y_pred_n'; y_pred_r'];


% index of non-missing observable variables
ind_obs = ~isnan(y_obs) & any(Gamma_full' ~= 0)';
Z_obs = Z_all(ind_obs,:);

Gamma_obs = Z_obs*Gamma_full; % measurement loadings for all variables
Delta = Gamma_obs*P_pred*Gamma_obs' + Z_obs*R*Z_obs';

Delta_inv = eye(sum(ind_obs))/Delta;

v = y_obs(ind_obs) - y_pred(ind_obs); % innovations
loglik_t(1) = - sum(ind_obs)*0.5*log(2*pi) - 0.5*log(det(Delta)) - 0.5*v'*Delta_inv*v;

x_upd(1,:) = (x_pred(:,1) + P_pred*Gamma_obs'*Delta_inv*v)';
P_upd = P_pred - P_pred*Gamma_obs'/(Gamma_obs*P_pred*Gamma_obs' + Z_obs*R*Z_obs')*Gamma_obs*P_pred;
x_std(1,:) = sqrt(diag(P_upd))';


J_smooth = zeros(K,K,T-1);

cond_num = 0;
%surv_tbexp_fit = nan(T,N_t);

for t = 2:T
    if any(isinf(Delta(:))) || min(real(eig(Delta))) <= 1e-16
        cond_num = 1;
        break
    end
    x_pred(:,t) = Mu + Phi*x_upd(t-1,:)';
    P_pred = Phi*P_upd*Phi' + Q;
    J_smooth(:,:,t-1) = P_upd*Phi'*pinv(P_pred); % J

        
    [yfit_all_n, yfit_all_r, JJ_n, JJ_r] = ...
        y_fitting(x_pred(:,t), A_X_for, B_X_for, A_X_exp, r_lb, s_n, A_X_for_pi, B_X_pi, Sigma2_X, B_X_cum, B_X_cum_pi);
    y_pred_n = yfit_all_n(mats_n*freq); % predicted nominal yields
    y_pred_r = yfit_all_r(mats_r*freq); % predicted real yields
    Gamma_n = JJ_n(mats_n*freq,:);
    Gamma_r = JJ_r(mats_r*freq,:);
    
    y_pred_t = nan(1,N_t);
    Gamma_t = zeros(N_t, K);
    %
    if any(~isnan(surv_tbexp(t,:)))
        yfit_tb_all = nan(max(hstep_t(~isnan(surv_tbexp(t,:)))),1);
        %yfit_tb_all = nan(max(hstep_t),1);
        %yfit_tb_all(1) = yfit_all_n(3);
        Gamma_tb_all = nan(max(hstep_t(~isnan(surv_tbexp(t,:)))), K);
        %Gamma_tb_all = nan(max(hstep_t), K);
        x_pred_h = x_pred(:,t);
        for j = 1:max(hstep_t(~isnan(surv_tbexp(t,:)))-1)
            x_pred_h = Mu + Phi*x_pred_h;
            if mod(j-5,12)==0
                [yfit_all_t, ~, JJ_t] = ...
                    y_fitting(x_pred_h, A_X_for(1:3), B_X_for(:,1:3), A_X_exp(1:3), r_lb, s_n(1:3), A_X_for_pi(1:3), B_X_pi(:,1:3), Sigma2_X, B_X_cum(:,1:3), B_X_cum_pi(:,1:3));
                yfit_tb_all(j+1,1) = yfit_all_t(3); % predicted 3-month T-bill
                Gamma_tb_all(j+1,:) = JJ_t(3,:)*Phi_j(:,:,j);
            end
        end
        y_pred_t = nan(size(hstep_t));
        Gamma_t = zeros(length(hstep_t), K);
        for j = 1:length(hstep_t(~isnan(surv_tbexp(t,:))))
            y_pred_t(j) = mean(yfit_tb_all(1:hstep_t(j)),'omitnan');
            Gamma_t(j,:) = mean(Gamma_tb_all(1:hstep_t(j),:),'omitnan');
        end
    end
    Gamma_full = [Gamma_m; Gamma_s; Gamma_g; Gamma_t; Gamma_n; Gamma_r];

    y_obs = y_obs_all(t,:)';
    y_pred = [(Gamma_m0 + Gamma_m*x_pred(:,t)); (Gamma_s0 + Gamma_s*x_pred(:,t)); ...
        (Gamma_g0 + Gamma_g*x_pred(:,t)); y_pred_t'; y_pred_n'; y_pred_r'];
    %ind_obs = ~isnan(y_obs);
    ind_obs = ~isnan(y_obs) & any(Gamma_full' ~= 0)';
    Z_obs = Z_all(ind_obs,:);

    Gamma_obs = Z_obs*Gamma_full; % measurement loadings for all variables
    v = y_obs(ind_obs) - y_pred(ind_obs); % innovations

    R_r = ones(N_r,1)*sigma_r_t(t)^2;
    R = diag([R_m; R_s; R_g; R_t; R_n; R_r]);
    Delta = Gamma_obs*P_pred*Gamma_obs' + Z_obs*R*Z_obs';
    Delta_inv = eye(sum(ind_obs))/Delta;

    loglik_t(t) = - sum(ind_obs)*0.5*log(2*pi) - 0.5*log(det(Delta)) - 0.5*v'*Delta_inv*v;

    
    x_upd(t,:) = (x_pred(:,t) + P_pred*Gamma_obs'*Delta_inv*v)';
    P_upd = P_pred - P_pred*Gamma_obs'/(Gamma_obs*P_pred*Gamma_obs' + Z_obs*R*Z_obs')*Gamma_obs*P_pred;
    x_std(t,:) = sqrt(diag(P_upd))';

end


loglik_t = -loglik_t;
%[(1:length(loglik_t))' loglik_t]
output_message = [];
if cond_num
    %fprintf(2,'Numerical condition error \n');
    output_message = 'Numerical condition error \n';
    %loglik_t = 1e+36;
elseif ind_error_sigma    
    %loglik_t = 1e+36;
end


end


%%
    function [estimates, fval] = estimate(x_00)
        fval_diff=1; estimates=x_00; k=0; fval=-fval0;
        %delete output.txt
        estimates_out = x_0;
        estimates_out(par_ind_full) = estimates;
        dlmwrite('output.txt',pars_trans_out(estimates_out, x_bounds), '-append', 'precision', 18); %saving starting values to a file

        while ~isfile('stop.txt') &&  fval_diff>0.00 %&& k<1
            k=k+1;
            fval_old = fval;
            x_00 = estimates;
            
            [estimates, fval] = fminsearch(@(x_1)maxlik_full(x_1,allpars), x_00, options_fminsearch);
            x_00 = estimates;
            %[estimates, fval] = fminunc(@(x_1)maxlik_full(x_1,allpars), x_00, options_fminunc);
            fval_diff = fval_old-fval;
            
            fprintf('%4.0f; ', k)
            fprintf('Updated likelihood :       %8.6f;\t', -fval);
            fprintf(';\t Log-lik. improvement:   %4.8f;', fval_diff);
            fprintf('\n');
            
            % Saving the output
            estimates_out = x_1;
            estimates_out(par_ind_full) = estimates;
            dlmwrite('output.txt',pars_trans_out(estimates_out, x_bounds), '-append', 'precision', 18); %saving results to a file
            % use: "output = dlmread('output.txt');" to read the results
            print_fit();
        end
    end

    function [yfit_all_n, yfit_all_r, JJ_n, JJ_r] = ...
            y_fitting_xxx(X0, A_X_for, B_X_for, A_X_exp, r_lb, s_n, A_X_for_pi, B_X_pi, Sigma2_X, B_X_cum, B_X_cum_pi)
        T0 = size(X0,2);
        Mm = length(A_X_for);

        y_fit_short = A_X_exp(1) + X0'*B_X_for(:,1); % shadow short rate
        Probs_short = ones(T0,1);
        Probs_short(y_fit_short < r_lb) = 0;
        y_fit_short(y_fit_short < r_lb) = r_lb; % observed short rate, T x 1
        
        mu = A_X_exp(2:end)' + X0'*B_X_for(:,2:end) - r_lb;  % T x maxmat-1
        z_n = ( mu ./ s_n(1:end-1)'); % T x maxmat-1
        Probs_long = mex_cdf(z_n); % T x maxmat-1
        pdf_z = normpdf(z_n); % T x maxmat-1
        
        % Nominal forward rates
        f_fit_long = r_lb + mu .* Probs_long + pdf_z .* s_n(1:end-1)' + Probs_long.*(A_X_for(2:end)'-A_X_exp(2:end)');
        %f_fit_long = r_lb + mu .* Probs_long + pdf_z .* s_n(1:end-1)' - 0.5*Probs_long.*sum(B_X_cum(:,2:end).*(Sigma2_X*B_X_cum(:,2:end)));
        f_fit = [y_fit_short f_fit_long]; % forward rates
        yfit_all_n = cumsum(f_fit, 2) ./ (1:Mm);
        
        % Expected inflation
        pi_fit = A_X_for_pi' + X0'*B_X_pi; % negative expected inflation

        % Real forward rates
        BXcumSig2BXcum_pi = sum(B_X_cum.*(Sigma2_X*B_X_cum_pi)); % 1 x maxmat
        r_fit = f_fit - pi_fit + [Probs_short Probs_long].*BXcumSig2BXcum_pi;
        yfit_all_r = cumsum(r_fit, 2) ./ (1:Mm);
        
        if nargout > 2
            % Jacobian
            JJ_f = zeros(K*T0,Mm);
            if A_X_exp(1) + X0'*B_X_for(:,1) > r_lb
                JJ_f(:,1) = B_X_for(:,1);
            end
            if Mm>1
                %JJ_f(2:end,:) = repmat(Probs_long,1,K) .* B_X_bar(:,2:end)'; % Jacobian for forwards
                JJ_f(:,2:end) = repelem(Probs_long,K,1) .* repmat(B_X_for(:,2:end),T0,1); % T*K x maxmat
            end
            JJ_n = cumsum(JJ_f,2)' ./ (1:Mm)'; %Jacobian for nominal yields
            JJ_r = cumsum((JJ_f - repmat(B_X_pi, T0, 1) + ...
                [zeros(T0*K,1) repelem(pdf_z,K,1).*repmat(B_X_for(:,2:end),T0,1).*BXcumSig2BXcum_pi(2:end) ] ), 2)' ./ (1:Mm)';
            %JJ_r = cumsum((JJ_f - repmat(B_X_pi, T0, 1)), 2)' ./ (1:Mm)'; % previous, incorrect version
        end
    end

    function [A_X_for, B_X_for, A_X_exp, B_X_cum, A_X_for_pi, B_X_pi, A_X_exp_pi, B_X_cum_pi, Conv_pi] = ...
            affine_coefs(Phi_q, Mu_q, Sigma2_X, delta_0, delta_1, delta_pi_0, delta_pi_1)
        
        B_X_for=zeros(K,maxmat); % forward rate loadings (starting from f0=y1) x freq
        B_X_for(:,1) = delta_1;
        A_X_exp = zeros(maxmat,1); % intercept for forward rates (starting from f0=y1) x freq
        A_X_exp(1) = delta_0;
        A_X_for = zeros(maxmat,1); % intercept for forward rates (starting from f0=y1) x freq
        A_X_for(1) = A_X_exp(1);
        B_X_cum = zeros(K,maxmat); % cumulative loadings
        
        for j=2:maxmat
            B_X_for(:,j) = Phi_q'*B_X_for(:,j-1);
            A_X_exp(j) = A_X_exp(j-1) + B_X_for(:,j-1)'*Mu_q;
            B_X_cum(:,j) = B_X_cum(:,j-1) + B_X_for(:,j-1);
            A_X_for(j) = A_X_exp(j) - 0.5*B_X_cum(:,j)'*Sigma2_X*B_X_cum(:,j)/freq;
        end
        
        B_X_pi = zeros(K,maxmat); % forward rate loadings (starting from f0=y1) x freq
        B_X_pi(:,1) = Phi_q'*delta_pi_1;
        A_X_exp_pi = zeros(maxmat,1); % intercept for forward rates (starting from f0=y1) x freq
        A_X_exp_pi(1) = delta_pi_0 + delta_pi_1'*Mu_q;
        B_X_cum_pi = zeros(K,maxmat); % cumulative of all B_X_for's starting from 0 
        B_X_cum_pi(:,1) = delta_pi_1; %
        Conv_pi = zeros(maxmat,1); % Convexity term for forward rates
        Conv_pi(1) = B_X_cum_pi(:,1)'*Sigma2_X*B_X_cum_pi(:,1)/freq;
        A_X_for_pi = zeros(maxmat,1); % intercept for forward rates (starting from f0=y1) x freq
        A_X_for_pi(1) = A_X_exp_pi(1) + 0.5*Conv_pi(1);
        
        for j=2:maxmat
            B_X_pi(:,j) = Phi_q'*B_X_pi(:,j-1);
            A_X_exp_pi(j) = A_X_exp_pi(j-1) + B_X_pi(:,j-1)'*Mu_q;
            B_X_cum_pi(:,j) = B_X_cum_pi(:,j-1) + B_X_pi(:,j-1);
            Conv_pi(j) = B_X_cum_pi(:,j)'*Sigma2_X*B_X_cum_pi(:,j)/freq;
            A_X_for_pi(j) = A_X_exp_pi(j) + 0.5*Conv_pi(j);
        end
    end
% Phi_q = Phi; Mu_q = Mu;
% B_X_pi_P = B_X_pi; A_X_exp_pi_P = A_X_exp_pi; Conv_pi_P = Conv_pi;

    function print_fit()
        x_smooth = x_upd';

        %
        surv_tbexp_fit = nan(T,N_t);
        y_pred_t = nan(1,N_t);
        obs_t = T;
        if ~isnan(surv_tbexp(obs_t,:))
            yfit_tb_all = zeros(max(hstep_t(~isnan(surv_tbexp(obs_t,:)))),1);
            x_pred_h = x_smooth(:,obs_t);
            [yfit_all_t] = ...
                y_fitting(x_pred_h, A_X_for(1:3), B_X_for(:,1:3), A_X_exp(1:3), r_lb, s_n(1:3), A_X_for_pi(1:3), B_X_pi(:,1:3), Sigma2_X, B_X_cum(:,1:3), B_X_cum_pi(:,1:3));
            %yfit_tb_all(1) = yfit_all_t(3);
            for j = 1:max(hstep_t(~isnan(surv_tbexp(obs_t,:)))-1)
                x_pred_h = Mu + Phi*x_pred_h;
                if mod(j-5,12)==0
                    [yfit_all_t] = ...
                        y_fitting(x_pred_h, A_X_for(1:3), B_X_for(:,1:3), A_X_exp(1:3), r_lb, s_n(1:3), A_X_for_pi(1:3), B_X_pi(:,1:3), Sigma2_X, B_X_cum(:,1:3), B_X_cum_pi(:,1:3));
                    yfit_tb_all(j+1,1) = yfit_all_t(3); % predicted 3-month T-bill
                end
            end
            for j = hstep_t(~isnan(surv_tbexp(obs_t,:)))
                y_pred_t(hstep_t == j) = mean(yfit_tb_all(1:j),'omitnan');
            end
        end
        surv_tbexp_fit(obs_t,:) = y_pred_t;
        %}
        for t=T-1:-1:1
            %x_smooth(:,t) = x_upd(t,:)' + J_smooth(:,:,t)*(x_smooth(:,t+1) - x_pred(:,t+1));
            x_smooth(:,t) = x_upd(t,:)';

            y_pred_t = nan(1,N_t);
            obs_t = t;
            if any(~isnan(surv_tbexp(obs_t,:)))
                yfit_tb_all = nan(max(hstep_t(~isnan(surv_tbexp(obs_t,:)))),1);
                x_pred_h = x_smooth(:,obs_t);
                [yfit_all_t] = ...
                    y_fitting(x_pred_h, A_X_for(1:3), B_X_for(:,1:3), A_X_exp(1:3), r_lb, s_n(1:3), A_X_for_pi(1:3), B_X_pi(:,1:3), Sigma2_X, B_X_cum(:,1:3), B_X_cum_pi(:,1:3));
                %yfit_tb_all(1) = yfit_all_t(3);
                for j = 1:max(hstep_t(~isnan(surv_tbexp(obs_t,:)))-1)
                    x_pred_h = Mu + Phi*x_pred_h;
                    if mod(j-5,12)==0
                        [yfit_all_t] = ...
                            y_fitting(x_pred_h, A_X_for(1:3), B_X_for(:,1:3), A_X_exp(1:3), r_lb, s_n(1:3), A_X_for_pi(1:3), B_X_pi(:,1:3), Sigma2_X, B_X_cum(:,1:3), B_X_cum_pi(:,1:3));
                        yfit_tb_all(j+1,1) = yfit_all_t(3); % predicted 3-month T-bill
                    end
                end
                for j = hstep_t(~isnan(surv_tbexp(obs_t,:)))
                    y_pred_t(hstep_t == j) = mean(yfit_tb_all(1:j),'omitnan');
                end
            end
            surv_tbexp_fit(obs_t,:) = y_pred_t;
        end


        x_upd = x_smooth';

        macro_fit = Gamma_m0' + x_upd*Gamma_m';
        surv_infexp_fit = Gamma_s0' + x_upd*Gamma_s';
        surv_gdpexp_fit = Gamma_g0' + x_upd*Gamma_g';

        % Fitted yields:
        [yfit_all_n, yfit_all_r] = ...
            y_fitting(x_upd', A_X_for, B_X_for, A_X_exp, r_lb, s_n, A_X_for_pi, B_X_pi, Sigma2_X, B_X_cum, B_X_cum_pi);
        yfit_n = yfit_all_n(:,mats_n*freq);
        yfit_r = yfit_all_r(:,mats_r*freq);

        % Expected yields:
        [yfit_all_n_exp, yfit_all_r_exp] = ...
            y_fitting(x_upd', A_X_for_P, B_X_for_P, A_X_exp_P, r_lb, s_n, A_X_for_pi_P, B_X_pi_P, Sigma2_X, B_X_cum_P, B_X_cum_pi_P);
        yfit_n_exp = yfit_all_n_exp(:,mats_n*freq);
        yfit_r_exp = yfit_all_r_exp(:,mats_r*freq);

        % ILS:
        Gamma_s_all_Q = cumsum(B_X_pi,2)./(1:maxmat); % Observable: expected inflation
        Gamma_s0_all_Q = cumsum(A_X_exp_pi)./(1:maxmat)' + 0.5*cumsum(Conv_pi)./(1:maxmat)'.^2; % Intercept 
        ils = x_upd*Gamma_s_all_Q + Gamma_s0_all_Q';

        % Expected inflation:
        Gamma_s_all_P = cumsum(B_X_pi_P,2)./(1:maxmat); % Observable: expected inflation
        Gamma_s0_all_P = cumsum(A_X_exp_pi_P)./(1:maxmat)' + 0.5*cumsum(Conv_pi_P)./(1:maxmat)'.^2; % Intercept 
        expinf = x_upd*Gamma_s_all_P + Gamma_s0_all_P';

        % Inflation risk premium:
        irp = ils - expinf;

        if verbose_fit
            u = yields_n - yfit_n;
            u2=u(2:end,:).^2;
            rmse=sqrt(nanmean(u2))*10000*12;
            fprintf('Maturities:              \t\t'); fprintf('%8.0f', mats_n); fprintf('\n');
            fprintf('Nominal yields RMSE (bp):\t\t ');
            fprintf('%8.2f', rmse);
            fprintf(';\t Average RMSE: %3.2f', mean(rmse))
            fprintf('\n');
    
            ind_r;
            u = yields_r - yfit_r;
            u2=u(2:end,:).^2;
            rmse=sqrt(nanmean(u2(~ind_r(2:end),:)))*10000*12;
            fprintf('Real yields - Normal liquidity RMSE (bp):\t\t ');
            fprintf('%8.2f', rmse);
            fprintf(';\t Average RMSE: %3.2f', mean(rmse))
            fprintf('\n');
    
    
            rmse=sqrt(nanmean(u2(ind_r(2:end),:)))*10000*12;
            fprintf('Real yields - High liquidity RMSE (bp):\t\t\t ');
            fprintf('%8.2f', rmse);
            fprintf(';\t Average RMSE: %3.2f', mean(rmse))
            fprintf('\n');
    
            rmse=sqrt(nanmean(u2))*10000*12;
            fprintf('Average real RMSE (bp):\t\t\t\t\t\t\t ');
            fprintf('%8.2f', rmse);
            fprintf(';\t Average RMSE: %3.2f', mean(rmse))
            fprintf('\n');
    
    
            u = surv_infexp - surv_infexp_fit;
            u2=u(2:end,:).^2;
            rmse=sqrt(nanmean(u2))*10000*12;
            fprintf('Fitting exp.inf. RMSE (bp):\t\t\t    ');
            fprintf('%8.2f\t\t\t\t\t\t', rmse(1));fprintf(' %8.2f', rmse(2));% fprintf('%8.2f', rmse(3));
            fprintf(';\t Average RMSE: %3.2f', mean(rmse))
            fprintf('\n');
    
    
            u = surv_gdpexp - surv_gdpexp_fit;
            u2=u(2:end,:).^2;
            rmse=sqrt(nanmean(u2))*10000*12;
            fprintf('Fitting exp.GDP. RMSE (bp):\t\t\t    ');
            fprintf('%8.2f\t\t\t\t\t\t', rmse(1)); fprintf(' %8.2f', rmse(2));
            fprintf(';\t Average RMSE: %3.2f', mean(rmse))
            fprintf('\n');
    
            u = surv_tbexp - surv_tbexp_fit;
            u2=u(2:end,:).^2;
            rmse=sqrt(nanmean(u2))*10000*12;
            fprintf('Fitting exp.3-m rate RMSE (bp):\t\t\t');
            fprintf('%8.2f\t\t\t\t\t\t', rmse(1)); fprintf(' %8.2f', rmse(2));
            fprintf(';\t Average RMSE: %3.2f', mean(rmse))
            fprintf('\n');
        end

    end

end
