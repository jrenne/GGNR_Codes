
set(0, 'DefaultFigureRenderer', 'painters'); figure('Position', [460 330 730 423]);
colors = get(gca, 'colororder');

clear p
fill([tT' flip(tT)'],...
    [[(x_upd(:,1) - 1.96*x_std(:,1) )*1200]' [flip((x_upd(:,1) + 1.95*x_std(:,1) )*1200)]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
p(1) = plot(tT, ffr, 'k', 'Linewidth',1.2);
p(2) = plot(tT, x_upd(:,1)*1200, 'Color', colors(2,:), 'Linewidth',1.2);
line([tT(1) tT(end)], [0 0],'color','k','Linewidth',0.01); hold off
grid on; %ylim([-12 17])
ylabel('%', 'Rotation', 0)
legend(p, 'Fed funds rate', '$s_t$', 'Interpreter','Latex')
%title('$r^*_t$', 'Interpreter','Latex')


fig = gcf;
%print(fig, "-dpdf", "-painters", "../Figures/fig_shadow_rate_revision.eps")
save_pdf_figure(fig, "../Figures/fig_shadow_rate_revision.pdf")

