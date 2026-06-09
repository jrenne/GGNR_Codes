

set(0, 'DefaultFigureRenderer', 'painters');
figure('Position', [150 200 850 650]);

subplot(3,1,1)
fill([tT' flip(tT)'],...
    [[(x_upd(:,8) - 1.96*x_std(:,8) )*1200]' [flip((x_upd(:,8) + 1.96*x_std(:,8) )*1200)]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
p = plot(tT, x_upd(:,8)*1200, tT(~isnan(r_LW)), r_LW(~isnan(r_LW)),  tT(~isnan(r_HLW)), r_HLW(~isnan(r_HLW)), 'Linewidth',1.2);
%plot(tT, (x_upd(:,8) + 2*x_std(:,8) )*1200,'k--', tT, (x_upd(:,8) - 2*x_std(:,8) )*1200,'k--', 'Linewidth',0.1); hold off
line([tT(1) tT(end)], [0 0],'color','k','Linewidth',0.01); hold off
grid on; %ylim([-12 17])
ylim([-8 10])
ylabel('%', 'Rotation', 0)
legend(p, 'GGNR', 'LW', 'HLW')
title('$r^*_t$', 'Interpreter','Latex')

subplot(3,1,2)
fill([tT' flip(tT)'],...
    [[(x_upd(:,7) - 1.96*x_std(:,7) )*1200]' [flip((x_upd(:,7) + 1.96*x_std(:,7) )*1200)]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
plot(tT, x_upd(:,7)*1200,'Linewidth',1.2);
%plot(tT, (x_upd(:,7) + 2*x_std(:,7) )*1200,'k--', tT, (x_upd(:,7) - 2*x_std(:,7) )*1200,'k--', 'Linewidth',0.1); hold off
line([tT(1) tT(end)], [0 0],'color','k','Linewidth',0.01)
grid on; %ylim([-12 17])
ylabel('%', 'Rotation', 0)
title('$\pi^*_t$', 'Interpreter','Latex')

subplot(3,1,3)
fill([tT' flip(tT)'],...
    [[(x_upd(:,end) - 1.96*x_std(:,end) )*1]' [flip((x_upd(:,end) + 1.96*x_std(:,end) )*1)]'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
plot(tT, x_upd(:,end)*1,'Linewidth',1.2);
line([tT(1) tT(end)], [0 0],'color','k','Linewidth',0.01)
grid on; %ylim([-12 17])
%ylabel('%', 'Rotation', 0)
title('$w_t$', 'Interpreter','Latex')

fig = gcf;
%print(fig, "-dpdf", "-painters", "../Figures/fig_stars_revision.eps")
save_pdf_figure(fig, "../Figures/fig_stars_revision.pdf")
