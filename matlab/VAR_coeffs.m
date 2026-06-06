function [A, Mu_til, Phi_til, Sigma_til, Lambda_0, Lambda_1, pars] = VAR_coeffs(v, ind_sig_w, x_bounds)

if nargin<2
    ind_sig_w = 0;
end

K = 13;
v = pars_trans_out(v, x_bounds);

% Parameters 
% Equation i (12):
pars.rho_i = v(1);
pars.alpha_pi = v(2);
pars.alpha_z = v(3);
pars.sigma_i = v(4);

% Equation g (15):
pars.rho_g = v(5);
pars.mu_g = v(6);
pars.sigma_g = v(7);

% Equation z (19):
pars.alpha = v(8);
pars.rho_z = v(9);
pars.sigma_z = v(10);

% Equation pi_st (14):
pars.rho_st = v(11);
pars.mu_st = v(12);
pars.sigma_st = v(13);

% Equation r_st (13):
pars.mu_kappa = v(14);
pars.rho_kappa = v(15);

%sigma_kappa = v(16);
pars.sigma_kappa = sqrt(0.00056^2 - pars.sigma_g^2);
v(16) = pars.sigma_kappa;

pars.theta = v(17);

% Equation pi_bar (17):
pars.rho_bar = v(18);
pars.beta = v(19);
pars.sigma_pi = v(20);

% Equation pi (18)
a0 = v(21);
a1 = v(22);
pars.a_u = [a0 a1 0]';

pars.r_lb = v(23);

% Equation w (20)
%pars.rho_w = v(24);
pars.rho_w = v(24);
pars.mu_w = v(25); %0;
if ind_sig_w==0
    pars.sigma_w = sqrt(1-pars.rho_w^2);
end
v(26) = pars.sigma_w;

pars.sigma_wg = v(27);
pars.sigma_wz = v(28);
pars.sigma_wpi = v(29);

% Measurement errors
pars.sigma_o = v(30); % Output gap
pars.sigma_gdp = v(31); % GDP growth
pars.sigma_inf = v(32); % Inflation
pars.sigma_ptr = v(33); % Perceived target rate
pars.sigma_s = v(34); % Inflation survey
pars.sigma_y = v(35); % GDP growth survey
pars.sigma_t = v(36); % Expected inflation survey
pars.sigma_n = v(37); % Nominal yields
pars.sigma_r = v(38); % Real yields
pars.sigma_r_liq = v(39); % Real yields - liquidity premium

% Equation lambda (21)
pars.Lambda_0 = v(40:47)';
pars.Lambda_z = v(48:55)';
pars.Lambda_pi = v(56:63)';
pars.Lambda_w = v(64:71)';

Lambda_1 = [zeros(8,4) pars.Lambda_z zeros(8,3) pars.Lambda_pi zeros(8,3) pars.Lambda_w];
Lambda_0 = pars.Lambda_0;

A = eye(K);
A(1,5) = -pars.rho_i*pars.alpha_z;
A(1,7) = -pars.rho_i*(1-pars.alpha_pi);
A(1,8) = -pars.rho_i;
A(1,9) = -pars.rho_i*pars.alpha_pi;
A(1,10:10+length(pars.a_u)-1) = -pars.rho_i*pars.alpha_pi*pars.a_u';
A(8,3) = -1;
A(9,7) = -(1 - pars.rho_bar);
%A(end,7) = 1;

Mu_til = zeros(K,1);
Mu_til(3) = (1-pars.rho_kappa)*pars.mu_kappa;
Mu_til(4) = (1-pars.rho_g)*pars.mu_g;
Mu_til(5) = pars.alpha*(1-pars.rho_bar)*(1-pars.rho_st)*pars.mu_st;
Mu_til(7) = (1-pars.rho_st)*pars.mu_st;
Mu_til(8) = pars.mu_kappa + pars.theta*(1-pars.rho_g)*pars.mu_g;
Mu_til(end) = (1-pars.rho_w)*pars.mu_w;

Phi_til = zeros(K);
Phi_til(1,1) = 1-pars.rho_i;
Phi_til(2,1) = 1;
Phi_til(3,3) = pars.rho_kappa;
Phi_til(4,4) = pars.rho_g;
Phi_til(5,1) = -pars.alpha;
Phi_til(5,5) = pars.alpha*pars.beta + pars.rho_z;
Phi_til(5,7) = pars.alpha*(1-pars.rho_bar)*pars.rho_st;
Phi_til(5,8) = pars.alpha;
Phi_til(5,9) = pars.alpha*pars.rho_bar;
Phi_til(5,10:10+length(pars.a_u)-1) = pars.alpha*[pars.a_u(2:end)' 0];
Phi_til(6,5) = 1;
Phi_til(7,7) = pars.rho_st;
Phi_til(8,4) = pars.theta*pars.rho_g;
Phi_til(9,5) = pars.beta;
Phi_til(9,9) = pars.rho_bar;
Phi_til(10:10+length(pars.a_u)-1,10:10+length(pars.a_u)-1) = diag(ones(length(pars.a_u)-1,1),-1);
Phi_til(end, end) = pars.rho_w;

% Vector of unconditional means:
Mu_unc = zeros(K,1);
Mu_unc(1) = pars.mu_kappa + pars.theta*pars.mu_g + pars.mu_st;
Mu_unc(2) = Mu_unc(1);
Mu_unc(3) = pars.mu_kappa;
Mu_unc(4) = pars.mu_g;
Mu_unc(7) = pars.mu_st;
Mu_unc(8) = pars.mu_kappa + pars.theta*pars.mu_g;
Mu_unc(9) = pars.mu_st;
Mu_unc(K) = pars.mu_w;


Sigma_til = zeros(K,8);
Sigma_til(1,1) = pars.sigma_i;
Sigma_til(3,8) = pars.sigma_kappa;
Sigma_til(4,2) = pars.sigma_g;
Sigma_til(5,3) = pars.sigma_z;
Sigma_til(7,4) = pars.sigma_st;
Sigma_til(8,2) = pars.theta*pars.sigma_g;
Sigma_til(9,5) = pars.sigma_pi;
Sigma_til(10,6) = 1;
Sigma_til(end,2) = pars.sigma_wg;
Sigma_til(end,3) = pars.sigma_wz;
Sigma_til(end,5) = pars.sigma_wpi;
Sigma_til(end,7) = pars.sigma_w;

end
