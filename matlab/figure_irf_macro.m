function [imp_res] = figure_irf_macro(x_0, lags)

% Impulse response function for the GGNR model
% v - vector of parameters
% X0 - state vector only for interest rates
% lags - number of lags
% k - which shock?
x_bounds = x_bounds_vals;
v = pars_trans_in(x_0, x_bounds);

%%
K = 13; %number of primary state variables: [s_t s_{t-1} kappa g z_t z_{t-1} pi_st r_st pi_bar eps_{u,t} eps_{u,t-1} eps_{u,t-2} w]'

[A, ~, Phi_til, Sigma_til, ~, ~, pars] = VAR_coeffs(v,0, x_bounds);

Phi = A\Phi_til;
Sigma_X = A\Sigma_til;

Id = eye(size(Sigma_X,2)); % 8 possible shocks: e_i, e_g, e_z, e_st, e_pi, e_u, e_w, e_kappa

imp_res_template = nan(K,8,lags+1);
imp_res_template(:,:,1) = Sigma_X;
for j = 1:lags
    imp_res_template(:,:,j+1) = Phi*imp_res_template(:,:,j);
end

set(0, 'DefaultFigureRenderer', 'painters'); figure('Position', [150 100 950 850]);

subplot(5,3,1);
plot(0:lags, squeeze(imp_res_template(9,1,:))*1200, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$\bar{\pi}_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon_{i,t}$', 'Interpreter','latex')

subplot(5,3,2);
plot(0:lags, squeeze(imp_res_template(5,1,:))*1200, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$z_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon_{i,t}$', 'Interpreter','latex')

subplot(5,3,3);
plot(0:lags, squeeze(imp_res_template(1,1,:))*1200, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$s_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon_{i,t}$', 'Interpreter','latex')


subplot(5,3,4);
plot(0:lags, squeeze(imp_res_template(9,2,:))*1200, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$\bar{\pi}_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon_{g,t}$', 'Interpreter','latex')

subplot(5,3,5);
plot(0:lags, squeeze(imp_res_template(5,2,:))*1200, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$z_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon_{g,t}$', 'Interpreter','latex')

subplot(5,3,6);
plot(0:lags, squeeze(imp_res_template(1,2,:))*1200, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$s_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon_{g,t}$', 'Interpreter','latex')



subplot(5,3,7);
plot(0:lags, squeeze(imp_res_template(9,3,:))*100, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$\bar{\pi}_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon_{z,t}$', 'Interpreter','latex')

subplot(5,3,8);
plot(0:lags, squeeze(imp_res_template(5,3,:))*100, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$z_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon_{z,t}$', 'Interpreter','latex')

subplot(5,3,9);
plot(0:lags, squeeze(imp_res_template(1,3,:))*100, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$s_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon_{z,t}$', 'Interpreter','latex')


subplot(5,3,10);
plot(0:lags, squeeze(imp_res_template(9,4,:))*1200, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$\bar{\pi}_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon^*_{t}$', 'Interpreter','latex')

subplot(5,3,11);
plot(0:lags, squeeze(imp_res_template(5,4,:))*1200, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$z_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon^*_{t}$', 'Interpreter','latex')

subplot(5,3,12);
plot(0:lags, squeeze(imp_res_template(1,4,:))*1200, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$s_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon^*_{t}$', 'Interpreter','latex')


subplot(5,3,13);
plot(0:lags, squeeze(imp_res_template(9,5,:))*1200, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$\bar{\pi}_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon_{\pi,t}$', 'Interpreter','latex')

subplot(5,3,14);
plot(0:lags, squeeze(imp_res_template(5,5,:))*1200, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$z_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon_{\pi,t}$', 'Interpreter','latex')

subplot(5,3,15);
plot(0:lags, squeeze(imp_res_template(1,5,:))*1200, 'LineWidth', 1.5)
grid on; xlim([0 lags])
ylabel('$s_t$', 'Rotation',0, 'Interpreter','latex')
title('$\varepsilon_{\pi,t}$', 'Interpreter','latex')


fig = gcf;
%print(fig, "-dpdf", "-painters", "../Figures/fig_irf_macro.eps")
save_pdf_figure(fig, "../Figures/fig_irf_macro.pdf")




