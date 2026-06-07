function [irf_n, irf_r]=figure_irf_yields(estimates, mats2plot, lags, delta, b, var_txt)

%mats2plot = [3 24 120];


Gamma = [];
for j = 1:length(var_txt)
    if strcmp(var_txt{j}, 's_t')
        Gamma = [Gamma; [1 zeros(1,12)]*12];
    elseif strcmp(var_txt{j}, 'z_t')
        Gamma = [Gamma; zeros(1,4) 1 zeros(1,8)];
    elseif strcmp(var_txt{j}, 'r^*_t')
        Gamma = [Gamma; [zeros(1,7) 1 zeros(1,5)]*12];
    elseif strcmp(var_txt{j}, 'r^*_t+\pi^*_t')
        Gamma = [Gamma; [zeros(1,6) 1 1 zeros(1,5)]*12];
    else
        error(['What is the conditioning variable ', var_txt{j}, '?'])
    end
end

lines = {"-", "--", "-.", ":"};
set(0, 'DefaultFigureRenderer', 'painters'); figure('Position', [250 300 880 420]);
colors = get(gca, 'colororder');
legend_txt = {};

for j = 1:size(b,2)
    [irf_n_all, irf_r_all] = imp_res_cond1(estimates, lags, Gamma, b(:,j), [delta/12; zeros(7,1)]);
    irf_n = irf_n_all(mats2plot,:);
    irf_r = irf_r_all(mats2plot,:);
    cond_txt = ['$',var_txt{1}, '=', num2str(b(1,j)*100), '\%$'];
    for k = 2:size(b,1)
        cond_txt = [cond_txt, ',$', var_txt{k}, '=', num2str(b(k,j)*100), '\%$'];
    end
    
    subplot(1,2,1);
    for k = 1:length(mats2plot)
        plot(0:lags, irf_n(k,:)*1200, 'LineStyle', lines{j}, 'color', colors(k,:), 'LineWidth',1.5); hold on
        legend_txt = [legend_txt, [num2str(mats2plot(k)) ' months, ' cond_txt]];
    end
    grid on; axis tight;
    title('Nominal yields', 'Interpreter','Latex')

    subplot(1,2,2);
    for k = 1:length(mats2plot)
        plot(0:lags, irf_r(k,:)*1200, 'LineStyle', lines{j}, 'color', colors(k,:), 'LineWidth',1.5); hold on
    end
    grid on; axis tight; 
    title('Real yields', 'Interpreter','Latex')

    %subplot(size(b,2),2,j*2-1); ylim([min([irf_n(:); irf_r(:)]), max([irf_n(:); irf_r(:)])]*1200)
    %subplot(size(b,2),2,j*2); ylim([min([irf_n(:); irf_r(:)]), max([irf_n(:); irf_r(:)])]*1200)
end

% making space for the legend
subplot(1,2,1);
plot([0 lags(end)], [0 0], 'color', [0.7 0.7 0.7], 'LineWidth', 0.5) % zero line
[ylimits1] = ylim;
subplot(1,2,2);
plot([0 lags(end)], [0 0], 'color', [0.7 0.7 0.7], 'LineWidth', 0.5) % zero line
[ylimits2] = ylim;
y_min = min([ylimits1 ylimits2]);
y_max = max([ylimits1 ylimits2]);

if delta>0
    subplot(1,2,1);
    ylim([y_min 0.35])
    subplot(1,2,2);
    ylim([y_min 0.35])
else
    subplot(1,2,1);
    ylim([-0.35 y_max])
    subplot(1,2,2);
    ylim([-0.35 y_max])
end

% Legend in the top right figure
if delta > 0
    subplot(1,2,2); legend(legend_txt, 'Interpreter', 'latex', 'Location', 'NE')
else
    subplot(1,2,2); legend(legend_txt, 'Interpreter', 'latex', 'Location', 'SE')
end

%yt = yticks;
% Standardizing the ticks across all figures
%for j = 1:size(b,2)
%    subplot(size(b,2),2,j*2-1); yticks([yt])
%    subplot(size(b,2),2,j*2);  yticks([yt])
%end

% Saving the figure:
%fig = gcf; print(fig, "-dpdf", "-painters", "../Figures/fig_irf_s=0bp,r_pi.pdf")
%fig = gcf; print(fig, "-dpdf", "-painters", "../Figures/fig_irf_s=25bp,r_pi.pdf")
%fig = gcf; print(fig, "-dpdf", "-painters", "../Figures/fig_irf_s=0,r.pdf")
%fig = gcf; print(fig, "-dpdf", "-painters", "../Figures/fig_irf_s.pdf")
end
