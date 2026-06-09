%%
set(0, 'DefaultFigureRenderer', 'painters'); 
%figure('Position', [100 50 1200 700]);
figure('Position', [100 50 700 700]);

colors_def = get(gca,'colororder');

subplot(2,2,1); 
plot(tT, macro(:,2)*1200); hold on; plot(tT, x_upd(:,9)*1200, tT, x_upd(:,7)*1200, 'LineWidth',1.5) % inflation persistent component, pi_bar
plot(tT, macro(:,2)*1200,'color', colors_def(1,:)); grid on
line([tT(1) tT(end)], [0 0],'color','k','Linewidth',0.01)
ylabel('%', 'Rotation', 0)
title('Inflation', 'Interpreter','Latex')
legend('Inflation', '$\bar{\pi}$', '$\pi^*$', ...
    'Location','S','Interpreter','Latex')

subplot(2,2,2); 
plot(tT, macro(:,4)*1200,'x', 'LineWidth',1.5); hold on;
plot(tT(1:3:end), x_upd(1:3:end,7)*1200, 'color', colors_def(2,:), 'LineWidth',1.5) % Perceived Inflation Target Rate, pi^8
%line([tT(1) tT(end)], [0 0],'color','k','Linewidth',0.01)
grid on
ylabel('%', 'Rotation', 0)
title('Perceived Inflation Target Rate', 'Interpreter','Latex')
legend('Observed PTR', '$\pi^*$', ...
    'Location','NE','Interpreter','Latex')

subplot(2,2,3); 
plot(tT, macro(:,1)*1200, '-', 'Color',colors_def(1,:), 'LineWidth', 1); hold on
plot(tT, x_upd(:,4)*1200, 'Color', colors_def(2,:), 'LineWidth',1.5); hold off
line([tT(1) tT(end)], [0 0],'color','k','Linewidth',0.01)
grid on; ylim([-12 17])
ylabel('%', 'Rotation', 0)
title('GDP growth', 'Interpreter','Latex')
legend('GDP growth ', '$g_t$', 'Interpreter','Latex', ...
    'Location','S','Interpreter','Latex')

subplot(2,2,4); 
plot(tT, macro(:,3)*100,'x', tT(2:3:end), x_upd(2:3:end,5)*100, 'LineWidth',1.5) % Output gap
line([tT(1) tT(end)], [0 0],'color','k','Linewidth',0.01)
grid on
ylabel('%', 'Rotation', 0)
title('Output gap', 'Interpreter','Latex')
legend('Observed', '$z_t$', 'Location','S', 'Interpreter','Latex')


fig = gcf;
%print(fig, "-dpdf", "-painters", "../Figures/fig_macro_vars_revision.eps")
save_pdf_figure(fig, "../Figures/fig_macro_vars_revision.pdf")

