
x_bounds = x_bounds_vals();
v = pars_trans_in(x0, x_bounds);
[~, ~, ~, ~, ~, ~, pars] = VAR_coeffs(v,0, x_bounds);

set(0, 'DefaultFigureRenderer', 'painters');
figure('Position', [100 50 700 700]);
%figure('Position', [150 100 850 850]);

colors_def = get(gca,'colororder');


subplot(2,2,1); 
fill([tT' flip(tT)'],...
    [[(surv_infexp_fit(:,1) - 2*pars.sigma_s)*1200]' [flip(surv_infexp_fit(:,1) + 2*pars.sigma_s)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
p = plot(tT(~isnan(surv_infexp_int(:,1))), surv_infexp_int(~isnan(surv_infexp_int(:,1)),1)*1200, ...
    tT, surv_infexp_fit(:,1)*1200, 'LineWidth',1); hold off
%plot(tT, (surv_infexp_fit(:,1) + 2*x0(26))*1200, 'k--', tT, (surv_infexp_fit(:,1) - 2*x0(26))*1200, 'k--', 'LineWidth',0.1); hold off
grid on; ylim_1y = ylim;
title('SPF expected inflation: 1-year', 'Interpreter','Latex')
ylabel('%', 'Rotation', 0)
legend(p, 'Observed', 'Fitted', 'Interpreter','Latex')

subplot(2,2,2);
fill([tT' flip(tT)'],...
    [[(surv_infexp_fit(:,2) - 2*pars.sigma_s)*1200]' [flip(surv_infexp_fit(:,2) + 2*pars.sigma_s)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
plot(tT(~isnan(surv_infexp_int(:,2))), surv_infexp_int(~isnan(surv_infexp_int(:,2)),2)*1200, ...
    tT, surv_infexp_fit(:,2)*1200, 'LineWidth',1); hold off
grid on; ylim(ylim_1y)
ylabel('%', 'Rotation', 0)
title('SPF expected inflation: 10-year', 'Interpreter','Latex')
%{
subplot(2,2,3);
fill([tT' flip(tT)'],...
    [[(surv_gdpexp_fit(:,1) - 2*pars.sigma_y)*1200]' [flip(surv_gdpexp_fit(:,1) + 2*pars.sigma_y)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
plot(tT(~isnan(surv_gdpexp(:,1))), surv_gdpexp(~isnan(surv_gdpexp(:,1)),1)*1200, ...
    tT, surv_gdpexp_fit(:,1)*1200, 'LineWidth',1); hold off
grid on; 
ylim_1y = ylim;
title('SPF GDP growth: 1-year', 'Interpreter','Latex')
ylabel('%', 'Rotation', 0)
%legend('Observed', 'Fitted', 'Interpreter','Latex')

subplot(3,2,4);
fill([tT' flip(tT)'],...
    [[(surv_gdpexp_fit(:,2) - 2*pars.sigma_y)*1200]' [flip(surv_gdpexp_fit(:,2) + 2*pars.sigma_y)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
plot(tT(~isnan(surv_gdpexp(:,2))), surv_gdpexp(~isnan(surv_gdpexp(:,2)),2)*1200, ...
    tT, surv_gdpexp_fit(:,2)*1200, 'LineWidth',1); hold off
grid on; 
ylim(ylim_1y)
ylabel('%', 'Rotation', 0)
title('SPF GDP growth: 10-year', 'Interpreter','Latex')
%legend('Observed', 'Fitted', 'Interpreter','Latex')
%}
subplot(2,2,3);
fill([tT(~isnan(surv_tbexp(:,1)))' flip(tT(~isnan(surv_tbexp(:,1))))'],...
    [[(surv_tbexp_fit(~isnan(surv_tbexp(:,1)),1) - 2*pars.sigma_t)*1200]' [flip(surv_tbexp_fit(~isnan(surv_tbexp(:,1)),1) + 2*pars.sigma_t)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
plot(tT(~isnan(surv_tbexp(:,1))), surv_tbexp(~isnan(surv_tbexp(:,1)),1)*1200, ...
    tT, surv_tbexp_fit(:,1)*1200, 'LineWidth',1); hold off
grid on; ylim_1y = ylim;
ylabel('%', 'Rotation', 0)
title('SPF 3-month T-bill rate: 1-year', 'Interpreter','Latex')
%legend('Observed', 'Fitted', 'Interpreter','Latex')

subplot(2,2,4);
fill([tT(~isnan(surv_tbexp_fit(:,2)))' flip(tT(~isnan(surv_tbexp_fit(:,2))))'],...
    [[(surv_tbexp_fit(~isnan(surv_tbexp_fit(:,2)),2) - 2*pars.sigma_t)*1200]' [flip(surv_tbexp_fit(~isnan(surv_tbexp_fit(:,2)),2) + 2*pars.sigma_t)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
plot(tT(~isnan(surv_tbexp_int(:,2))), surv_tbexp_int(~isnan(surv_tbexp_int(:,2)),2)*1200, ...
    tT, surv_tbexp_fit(:,2)*1200, 'LineWidth',1); hold off
grid on; ylim(ylim_1y)
ylabel('%', 'Rotation', 0)
title('SPF 3-month T-bill rate: 10-year', 'Interpreter','Latex')


fig = gcf;
%print(fig, "-dpdf", "-painters", "../Figures/fig_surv_revision.eps")
save_pdf_figure(fig, "../Figures/fig_surv_revision.pdf")
