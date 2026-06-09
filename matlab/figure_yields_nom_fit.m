
v = pars_trans_in(x0, x_bounds);
[~, ~, ~, ~, ~, ~, pars] = VAR_coeffs(v,0, x_bounds);

set(0, 'DefaultFigureRenderer', 'painters'); figure('Position', [100 50 700 700]);

subplot(3,2,1); 
fill([tT' flip(tT)'],...
    [[(yfit_n(:,1) - 2*pars.sigma_n)*1200]' [flip(yfit_n(:,1) + 2*pars.sigma_n)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
p = plot(tT(~isnan(yields_n(:,1))), yields_n(~isnan(yields_n(:,1)),1)*1200, tT, yfit_n(:,1)*1200,'LineWidth',1.); hold off
grid on
title('Nominal 3-month yield', 'Interpreter','Latex')
ylabel('%', 'Rotation', 0)
legend(p, 'Observed', 'Fitted', 'Interpreter','Latex')

subplot(3,2,2); 
fill([tT' flip(tT)'],...
    [[(yfit_n(:,2) - 2*pars.sigma_n)*1200]' [flip(yfit_n(:,2) + 2*pars.sigma_n)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
plot(tT(~isnan(yields_n(:,2))), yields_n(~isnan(yields_n(:,2)),2)*1200, tT, yfit_n(:,2)*1200,'LineWidth',1.); hold off
grid on
ylabel('%', 'Rotation', 0)
title('Nominal 1-year yield', 'Interpreter','Latex')
%legend('Observed', 'Fitted', 'Interpreter','Latex')

subplot(3,2,3); 
fill([tT' flip(tT)'],...
    [[(yfit_n(:,3) - 2*pars.sigma_n)*1200]' [flip(yfit_n(:,3) + 2*pars.sigma_n)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
plot(tT(~isnan(yields_n(:,3))), yields_n(~isnan(yields_n(:,3)),3)*1200, tT, yfit_n(:,3)*1200,'LineWidth',1.); hold off
grid on
ylabel('%', 'Rotation', 0)
title('Nominal 2-year yield', 'Interpreter','Latex')
%legend('Observed', 'Fitted', 'Interpreter','Latex')

subplot(3,2,4); 
fill([tT' flip(tT)'],...
    [[(yfit_n(:,4) - 2*pars.sigma_n)*1200]' [flip(yfit_n(:,4) + 2*pars.sigma_n)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
plot(tT(~isnan(yields_n(:,4))), yields_n(~isnan(yields_n(:,4)),4)*1200, tT, yfit_n(:,4)*1200,'LineWidth',1.); hold off
grid on
ylabel('%', 'Rotation', 0)
title('Nominal 3-year yield', 'Interpreter','Latex')
%legend('Observed', 'Fitted', 'Interpreter','Latex')

subplot(3,2,5); 
fill([tT' flip(tT)'],...
    [[(yfit_n(:,5) - 2*pars.sigma_n)*1200]' [flip(yfit_n(:,5) + 2*pars.sigma_n)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
plot(tT(~isnan(yields_n(:,5))), yields_n(~isnan(yields_n(:,5)),5)*1200, tT, yfit_n(:,5)*1200,'LineWidth',1.); hold off
grid on
ylabel('%', 'Rotation', 0)
title('Nominal 5-year yield', 'Interpreter','Latex')

subplot(3,2,6); 
fill([tT' flip(tT)'],...
    [[(yfit_n(:,7) - 2*pars.sigma_n)*1200]' [flip(yfit_n(:,7) + 2*pars.sigma_n)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
plot(tT(~isnan(yields_n(:,7))), yields_n(~isnan(yields_n(:,7)),7)*1200, tT, yfit_n(:,7)*1200,'LineWidth',1.); hold off
grid on
ylabel('%', 'Rotation', 0)
title('Nominal 10-year yield', 'Interpreter','Latex')



fig = gcf;
%print(fig, "-depsc", "-painters", "../Figures/fig_term_prem.eps")
%print(fig, "-dpdf", "-painters", "../Figures/fig_yields_nom_revision.eps")
save_pdf_figure(fig, "../Figures/fig_yields_nom_revision.pdf")

