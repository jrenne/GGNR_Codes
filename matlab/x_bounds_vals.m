function x_bounds = x_bounds_vals()
% bounds for transformed parameters

x_bounds = [...
    0 0.995; ... % rho_i
    1 3; ... % alpha_pi
    0 0.1; ... % alpha_z
    0 1; ... % sigma_i
    0 0.995; ... % rho_g
    0 1; ... % mu_g
    0 0.00056; ... % sigma_g
    0 1; ... % alpha
    0 0.995; ... % rho_z
    0 1; ... % sigma_z
    0 0.995; ... % rho_st
    -100 100; ... % mu_st
    0 1; ... % sigma_st
    0 1; ... % mu_kappa
    0 0.995; ... % rho_kappa
    -100 100; ... % sigma_kappa
    -100 100; ... % theta
    -0.5 0.995; ... % rho_bar
    0 1; ... % beta
    0 0.005; ... % sigma_pi
    0 1; ... % a0
    0 0.005; ... % a1
    -100 100; ... % r_lb
    0 1; ... % rho_w
    -100 100; ... % mu_w
    -100 100; ... % sigma_w
    -0.1 0.1; ... % sigma_wg
    -0.1 0.1; ... % sigma_wz
    -0.1 0.1; ... % sigma_wpi
    0 0.05; ... % sigma_o
    0 0.05; ... % sigma_gdp
    0 0.05; ... % sigma_inf
    0 0.05; ... % sigma_ptr
    0 0.05; ... % sigma_s
    0 0.05; ... % sigma_y
    0 0.05; ... % sigma_t
    0 0.05; ... % sigma_n
    0 0.05; ... % sigma_r
    0 0.05; ... % sigma_r_liq
    ones(32,1)*[-10 10]]; % lambdas
