function figure_prob_contourplot(x_0, b, var_txt)

if length(var_txt)~=2 || length(b) ~=2
    error('Make sure there are 2 conditioning variables!')
end

Gamma = [];
for j = 1:length(var_txt)
    if strcmp(var_txt{j}, 'r^*_t')
        Gamma = [Gamma; [zeros(1,7) 1 zeros(1,5)]];
        multiplier = 1200;
    elseif strcmp(var_txt{j}, '\pi^*_t')
        Gamma = [Gamma; [zeros(1,6) 1 zeros(1,6)]];
        multiplier = 1200;
    else
        error(['What is the conditioning variable ', var_txt{j}, '?'])
    end
end

K = 13;
x_bounds = x_bounds_vals;
x_1 = pars_trans_in(x_0, x_bounds);

r_lb = 0;
[A, Mu_til, Phi_til, Sigma_til] = VAR_coeffs(x_1, 0, x_bounds);

Mu = A\Mu_til;
Phi = A\Phi_til;
Sigma_X = A\Sigma_til;
Sigma2_X = Sigma_X*Sigma_X';

delta_1 = [1 zeros(1,K-1)]';
delta_0 = 0;

E_X = (eye(K) - Phi)\Mu; % Should be the same as Mu_unc
V_X = reshape((eye(K*K) - kron(Phi, Phi))\Sigma2_X(:), K,K);

var_s = delta_1'*V_X*delta_1 - delta_1'*(V_X*Gamma')/(Gamma*V_X*Gamma')*(Gamma*V_X)*delta_1;

[XX, YY] = meshgrid(b{1}/multiplier, b{2}/multiplier);

e_s = delta_0 + delta_1'*E_X ...
    + delta_1'*(V_X*Gamma')/(Gamma*V_X*Gamma')*[1; 0]*XX ...
    + delta_1'*(V_X*Gamma')/(Gamma*V_X*Gamma')*[0; 1]*YY...
    - delta_1'*(V_X*Gamma')/(Gamma*V_X*Gamma')*Gamma*E_X;

Probs = normcdf((r_lb - e_s) / sqrt(var_s));
contour_lines = [0.30 0.25 0.20 0.15 0.10 0.05 0.005 0.001];

set(0, 'DefaultFigureRenderer', 'painters'); figure('Position', [560 430 630 423]);
colors = get(gca, 'colororder');
contour(XX*multiplier, YY*multiplier, Probs, contour_lines, 'color',colors(1,:), 'ShowText','on', 'LineWidth',1.5)
grid on
ylabel('Inflation target, $\pi^*$', 'Interpreter','latex')
xlabel('Natural rate of interest, $r^*$', 'Interpreter','latex')
%title('Probabilities of hitting lower bound conditional on $r^*$ and $\pi^*$', 'Interpreter','latex')

fig = gcf;
%print(fig, "-dpdf", "-painters", "../Figures/fig_probs.eps")
save_pdf_figure(fig, "../Figures/fig_probs.pdf")


