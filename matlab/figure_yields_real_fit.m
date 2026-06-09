
v = pars_trans_in(x0, x_bounds);
[~, ~, ~, ~, ~, ~, pars] = VAR_coeffs(v,0, x_bounds);

%ind_r = tT<datetime('2004-01-01') | (tT>=datetime('2008-01-01') & tT<datetime('2010-01-01')); %higher liquidity index
ind_r = tT<datetime('2004-01-01') | (tT>=datetime('2008-01-01') & tT<datetime('2010-01-01')) ...
     | (tT>=datetime('2020-03-01') & tT<datetime('2021-03-01')); %higher liquidity index
sigma_r_t = ones(length(tT),1)*abs(x0(38));
sigma_r_t(ind_r) = abs(pars.sigma_r) + abs(pars.sigma_r_liq);

set(0, 'DefaultFigureRenderer', 'painters'); figure('Position', [100 50 700 700]);

subplot(2,2,1);
fill([tT' flip(tT)'],...
    [[(yfit_r(:,1) - 2*sigma_r_t)*1200]' [flip(yfit_r(:,1) + 2*sigma_r_t)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
p = plot(tT(~isnan(yields_r(:,1))), yields_r(~isnan(yields_r(:,1)),1)*1200, tT, yfit_r(:,1)*1200,'LineWidth',1.); 
hold off
%plot(tT, (yfit_r(:,1) + 2*sigma_r_t)*1200, 'k--', tT, (yfit_r(:,1) - 2*sigma_r_t)*1200, 'k--', 'LineWidth',0.1)
grid on
title('Real 2-year yield', 'Interpreter','Latex')
ylabel('%', 'Rotation', 0)
%legend(p, 'Observed', 'Fitted', 'Interpreter','Latex')

subplot(2,2,2);
fill([tT' flip(tT)'],...
    [[(yfit_r(:,2) - 2*sigma_r_t)*1200]' [flip(yfit_r(:,2) + 2*sigma_r_t)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
plot(tT(~isnan(yields_r(:,2))), yields_r(~isnan(yields_r(:,2)),2)*1200, tT, yfit_r(:,2)*1200,'LineWidth',1.); hold off
%plot(tT, (yfit_r(:,2) + 2*sigma_r_t)*1200, 'k--', tT, (yfit_r(:,2) - 2*sigma_r_t)*1200, 'k--', 'LineWidth',0.1)
grid on
ylabel('%', 'Rotation', 0)
title('Real 3-year yield', 'Interpreter','Latex')
%legend('Observed', 'Fitted', 'Interpreter','Latex')

subplot(2,2,3);
fill([tT' flip(tT)'],...
    [[(yfit_r(:,3) - 2*sigma_r_t)*1200]' [flip(yfit_r(:,3) + 2*sigma_r_t)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
plot(tT(~isnan(yields_r(:,3))), yields_r(~isnan(yields_r(:,3)),3)*1200, tT, yfit_r(:,3)*1200,'LineWidth',1.); hold off
%plot(tT, (yfit_r(:,3) + 2*sigma_r_t)*1200, 'k--', tT, (yfit_r(:,3) - 2*sigma_r_t)*1200, 'k--', 'LineWidth',0.1)
grid on
ylabel('%', 'Rotation', 0)
title('Real 5-year yield', 'Interpreter','Latex')
%legend('Observed', 'Fitted', 'Interpreter','Latex')

subplot(2,2,4);
fill([tT' flip(tT)'],...
    [[(yfit_r(:,4) - 2*sigma_r_t)*1200]' [flip(yfit_r(:,4) + 2*sigma_r_t)*1200]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
p = plot(tT(~isnan(yields_r(:,4))), yields_r(~isnan(yields_r(:,4)),4)*1200, tT, yfit_r(:,4)*1200,'LineWidth',1.); hold off
grid on
ylabel('%', 'Rotation', 0)
title('Real 10-year yield', 'Interpreter','Latex')
legend(p, 'Observed', 'Fitted', 'Interpreter','Latex', 'Location', 'SW')



fig = gcf;
%print(fig, "-depsc", "-painters", "../Figures/fig_term_prem.eps")
%print(fig, "-dpdf", "-painters", "../Figures/fig_yields_real_revision.eps")
save_pdf_figure(fig, "../Figures/fig_yields_real_revision.pdf")


