function [irf_n, irf_r] = imp_res_cond1(x_0, h_max, Gamma, b, delta)


[N_r, K] = size(Gamma); % Number of restrictions

r_lb = 0;
x_bounds = x_bounds_vals;
v = pars_trans_in(x_0, x_bounds);

[A, Mu_til, Phi_til, Sigma_til, Lambda_0, Lambda_1, pars] = VAR_coeffs(v, 0, x_bounds);

N_e = size(Sigma_til,2); % Number of shocks

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


% Affine coefficients for fitting yields:
[A_X_for, B_X_for, A_X_exp, B_X_cum, A_X_for_pi, B_X_pi, A_X_exp_pi, B_X_cum_pi, Conv_pi] = ...
            affine_coefs(Phi_q, Mu_q, Sigma2_X, delta_0, delta_1, delta_pi_0, delta_pi_1);

s2_n = cumsum(diag(B_X_for'*(Sigma2_X)*B_X_for)); %
s_n = sqrt(s2_n); %sigma_n 


% Unconditional distribution
E_X = (eye(K) - Phi)\Mu; % Should be the same as Mu_unc
V_X = (eye(K*K) - kron(Phi, Phi))\Sigma2_X(:);

% Conditional forecasts
E_X_cond1 = zeros(K,h_max+1);
E_X_cond0 = zeros(K,h_max+1);
V_X_cond = zeros(K^2, h_max+1);
E_X_aux1 = [Sigma_X Phi*reshape(V_X,K,K)*Gamma']*pinv([eye(N_e) zeros(N_e,N_r); zeros(N_r,N_e) Gamma*reshape(V_X,K,K)*Gamma'])*...
    ([delta; b] - [zeros(N_e,1); Gamma*E_X]);
E_X_aux0 = [Sigma_X Phi*reshape(V_X,K,K)*Gamma']*pinv([eye(N_e) zeros(N_e,N_r); zeros(N_r,N_e) Gamma*reshape(V_X,K,K)*Gamma'])*...
    ([zeros(N_e,1); b] - [zeros(N_e,1); Gamma*E_X]);
V_X_aux = [Sigma_X Phi*reshape(V_X,K,K)*Gamma']*pinv([eye(N_e) zeros(N_e,N_r); zeros(N_r,N_e) Gamma*reshape(V_X,K,K)*Gamma'])*...
    [Sigma_X Phi*reshape(V_X,K,K)*Gamma']';

irf_n = zeros(length(A_X_for), h_max+1);
irf_r = zeros(length(A_X_for), h_max+1);

for j = 1:h_max+1
    % Conditional forecasts
    E_X_cond1(:,j) = E_X + E_X_aux1;
    E_X_cond0(:,j) = E_X + E_X_aux0;
    V_X_cond(:,j) = V_X - V_X_aux(:);

    % Real and Nominal yields with shock
    [yfit_all_n, yfit_all_r, ~, ~, JJ_nn, JJ_rr] = ...
        y_fitting(E_X_cond1(:,j), A_X_for, B_X_for, A_X_exp, r_lb, s_n, A_X_for_pi, B_X_pi, Sigma2_X, B_X_cum, B_X_cum_pi);
    q_n1 = yfit_all_n + 0.5*V_X_cond(:,j)'*JJ_nn';
    q_r1 = yfit_all_r + 0.5*V_X_cond(:,j)'*JJ_rr';

    % Real and Nominal yields without shock
    [yfit_all_n, yfit_all_r, ~, ~, JJ_nn, JJ_rr] = ...
        y_fitting(E_X_cond0(:,j), A_X_for, B_X_for, A_X_exp, r_lb, s_n, A_X_for_pi, B_X_pi, Sigma2_X, B_X_cum, B_X_cum_pi);
    q_n0 = yfit_all_n + 0.5*V_X_cond(:,j)'*JJ_nn';
    q_r0 = yfit_all_r + 0.5*V_X_cond(:,j)'*JJ_rr';

    % Impulse reponse functions
    irf_n(:,j) = q_n1' - q_n0';
    irf_r(:,j) = q_r1' - q_r0';

    E_X_aux1 = Phi*E_X_aux1;
    E_X_aux0 = Phi*E_X_aux0;
    V_X_aux = Phi*V_X_aux*Phi';
end


end