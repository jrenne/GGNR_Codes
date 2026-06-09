function cond_corr_simulation(x0, x_upd, N, H, hstep_s, tT)

% computes conditional correlation between real rates and expected
% inflation (Mundell-Tobin effect) for maturities 3 months and 12 months

% requires having different coefficiens

%N = 10000; % number of simulations
%H = 12; % forecast horizon

%hstep_s = [3 12 120]; % Maturities for conditional correlations

K = 13;
r_lb = 0;
x_bounds = x_bounds_vals;
x1 = pars_trans_in(x0, x_bounds);

[A, Mu_til, Phi_til, Sigma_til, Lambda_0, Lambda_1, pars] = VAR_coeffs(x1, 0, x_bounds);

Mu = A\Mu_til;
Phi = A\Phi_til;
Sigma_X = A\Sigma_til;
Sigma2_X = Sigma_X*Sigma_X';


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

Gamma_s = zeros(length(hstep_s), K);
Gamma_s0 = zeros(length(hstep_s),1);
for j = 1:length(hstep_s)
    Gamma_s(j,:) = sum(B_X_pi_P(:,1:hstep_s(j)),2)'/hstep_s(j); % Observable: expected inflation
    Gamma_s0(j) = sum(A_X_exp_pi_P(1:hstep_s(j))/hstep_s(j) + 0.5*Conv_pi_P(1:hstep_s(j))/hstep_s(j)^2); % Intercept 
end


cond_corr = zeros(length(x_upd), length(hstep_s)); % conditional correlations

%cond_corr_3 = zeros(length(x_upd),1); %cond. corr. for 3-month maturities
%cond_corr_12 = zeros(length(x_upd),1); %cond. corr. for 12-month maturities

eps = randn(8,N,H);
for t = 1:length(x_upd)
    X0 = x_upd(t,:)';
    
    X_h = ones(K, N).*X0;
    for h = 1:H
        %eps = randn(8,N);
        X_h = Mu + Phi*X_h + Sigma_X*eps(:,:,h);
    end
    
    [yfit_all_n, yfit_all_r] = ...
                y_fitting(X_h, A_X_for, B_X_for, A_X_exp, r_lb, s_n, A_X_for_pi, B_X_pi, Sigma2_X, B_X_cum, B_X_cum_pi);
    
    ExpInf = Gamma_s0 + Gamma_s*X_h;
    
    for j = 1:length(hstep_s)
        rho = corr([ExpInf(j,:)' yfit_all_r(:,hstep_s(j))]);
        cond_corr(t,j) = rho(1,2);
    end
    %C_3 = corr([ExpInf(1,:)' yfit_all_r(:,3)]);
    %cond_corr_3(t) = C_3(1,2);
    %C_12 = corr([ExpInf(2,:)' yfit_all_r(:,12)]);
    %cond_corr_12(t) = C_12(1,2);
end

legend_txt = {};
for j = 1:length(hstep_s)
    legend_txt = [legend_txt, [num2str(hstep_s(j)) '-month maturity']];
end

figure;
plot(tT, cond_corr, 'LineWidth', 1.2)
grid on
legend(legend_txt,'Location','southwest')
%legend('3-month maturity', '12-month maturity','Location','southwest')
%title([{'Conditional correlation between real rates and expected inflation'}, {'at 12-month forecast horizon'}],'Interpreter', 'Latex')

fig = gcf;
%print(fig, "-dpdf", "-painters", "../Figures/fig_cond_corr.eps")
save_pdf_figure(fig, "../Figures/fig_cond_corr.pdf")

