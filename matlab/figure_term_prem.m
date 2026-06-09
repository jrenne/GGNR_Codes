%set(0, 'DefaultFigureRenderer', 'painters'); figure('Position', [100 50 700 700]);
set(0, 'DefaultFigureRenderer', 'painters'); 
%figure('Position', [150 100 950 850]);
figure('Position', [150 200 850 650]);

colors = get(gca, 'colororder');

subplot(3,1,1); 
plot(tT, term_prem_n(:,end)*1200, tT, term_prem_10y_KW, tT,term_prem_10y_ACM, ...
    tT(~isnan(surv_tbexp(:,2))), (yields_n(~isnan(surv_tbexp(:,2)),end) - surv_tbexp(~isnan(surv_tbexp(:,2)),2))*1200,'x', ...
    'LineWidth',1.)
line([tT(1) tT(end)], [0 0],'color','k','Linewidth',0.01)
grid on;
recessionplot;
title('Nominal 10y yield term premium', 'Interpreter','Latex')
legend('GGNR', 'KW', 'ACM','Model free')
ylabel('%', 'Rotation', 0)
%{
subplot(3,2,2); plot(tT(tT>'Dec-2006'), term_prem_n(tT>'Dec-2006',end)*1200, ...
    tT(tT>'Dec-2006'), term_prem_10y_KW(tT>'Dec-2006'), ...
    tT(tT>'Dec-2006'), term_prem_10y_ACM(tT>'Dec-2006'), ...
    tT(tT>'Dec-2006' & ~isnan(surv_tbexp(:,2))), (yields_n(tT>'Dec-2006' & ~isnan(surv_tbexp(:,2)),end) - surv_tbexp(tT>'Dec-2006' & ~isnan(surv_tbexp(:,2)),2))*1200,'x', ...
    'LineWidth',1.)
line([tT(tT=='Dec-2006') tT(end)], [0 0],'color','k','Linewidth',0.01)
grid on;
recessionplot;
title('Nominal 10y yield term premium, 2007 - 2023', 'Interpreter','Latex')
ylabel('%', 'Rotation', 0)
%legend('1-y term prem.', '5-y term prem.', '10-y term prem.', 'Interpreter','Latex')
%legend('10-y yield', '1-y term prem.', '3-y term prem.', '5-y term prem.', '10-y term prem.', 'Interpreter','Latex')
%annotation('textbox', [0.45 0.78 0.1 0.2], 'String','Method 1', 'EdgeColor','none')
%}

subplot(3,1,2); 
p1 = plot(tT(~isnan(yields_r(:,end))), term_prem_r(~isnan(yields_r(:,end)),end)*1200, 'LineWidth', 1.); hold on
p2 = plot(tT, term_prem_r_10y_DKW, 'LineWidth', 1.);
p3 = plot(tT(~isnan(surv_tbexp(:,2))), (yields_r(~isnan(surv_tbexp(:,2)),end) - surv_tbexp(~isnan(surv_tbexp(:,2)),2) + surv_infexp(~isnan(surv_tbexp(:,2)),2))*1200, ...
    'x', 'color', colors(4,:), 'LineWidth', 1.); hold off
line([tT(1) tT(end)], [0 0],'color','k','Linewidth',0.01)
grid on;
ylabel('%', 'Rotation', 0)
ylim([-2 6])
recessionplot;
legend([p1 p2 p3], 'GGNR', 'DKW', 'Model free')
title('Real 10y yield term premium', 'Interpreter','Latex')
%{
subplot(3,2,4); 
plot(tT(tT>'Dec-2006'), term_prem_r(tT>'Dec-2006',end)*1200, ...
    tT(tT>'Dec-2006'), term_prem_r_10y_DKW(tT>'Dec-2006',end), 'LineWidth', 1.); hold on
plot(tT(tT>'Dec-2006' & ~isnan(surv_tbexp(:,2))), (yields_r(tT>'Dec-2006' & ~isnan(surv_tbexp(:,2)),end) - surv_tbexp(tT>'Dec-2006' & ~isnan(surv_tbexp(:,2)),2) + surv_infexp(tT>'Dec-2006' & ~isnan(surv_tbexp(:,2)),2))*1200, ...
    'x', 'color', colors(4,:), 'LineWidth', 1.); hold off
line([tT(tT=='Dec-2006') tT(end)], [0 0],'color','k','Linewidth',0.01);
grid on;
recessionplot;
title('Real 10y yield term premium, 2007 - 2023', 'Interpreter','Latex')
ylabel('%', 'Rotation', 0)
%legend('2-y term prem.', '5-y term prem.', '7-y term prem.', '10-y term prem.', 'Interpreter','Latex')
%legend('10-y yield', '2-y term prem.', '5-y term prem.', '7-y term prem.', '10-y term prem.', 'Interpreter','Latex')
%annotation('textbox', [0.45 0.31 0.1 0.2], 'String','Method 2', 'EdgeColor','none')
%}

subplot(3,1,3);
plot(tT, irp(:,end)*1200, tT, inf_prem_r_10y_DKW, 'LineWidth', 1.); hold on
plot(tT(~isnan(surv_infexp(:,2))), (yields_n(~isnan(surv_infexp(:,2)),end) - yields_r(~isnan(surv_infexp(:,2)),end) - surv_infexp(~isnan(surv_infexp(:,2)),2))*1200, ...
    'x', 'color', colors(4,:), 'LineWidth', 1.); hold off
line([tT(1) tT(end)], [0 0],'color','k','Linewidth',0.01)
ylim([-2 4])
grid on;
recessionplot;
legend('GGNR', 'DKW', 'Model free')
title('ILS 10y inflation premium', 'Interpreter','Latex')
ylabel('%', 'Rotation', 0)
%{
subplot(3,2,6); plot( ...
    tT(tT>'Dec-2006'), irp(tT>'Dec-2006',end)*1200, tT(tT>'Dec-2006'), ...
    inf_prem_r_10y_DKW(tT>'Dec-2006',end), 'LineWidth', 1.); hold on
plot(tT(tT>'Dec-2006' & ~isnan(surv_infexp(:,2))), (yields_n(tT>'Dec-2006' & ~isnan(surv_infexp(:,2)),end) - yields_r(tT>'Dec-2006' & ~isnan(surv_infexp(:,2)),end) - surv_infexp(tT>'Dec-2006' & ~isnan(surv_infexp(:,2)),2))*1200, ...
    'x', 'color', colors(4,:), 'LineWidth', 1.); hold off
line([tT(tT=='Dec-2006') tT(end)], [0 0],'color','k','Linewidth',0.01)
grid on;
recessionplot;
title('ILS 10y inflation premium, 2007 - 2023', 'Interpreter','Latex')
ylabel('%', 'Rotation', 0)
%}

fig = gcf;
%print(fig, "-dpdf", "-painters", "../Figures/fig_term_prem_revision.eps")
save_pdf_figure(fig, "../Figures/fig_term_prem_revision.pdf")
