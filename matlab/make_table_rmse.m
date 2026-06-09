% =========================================================================
% Create Latex table of the fit of yields
% =========================================================================

% Save table? 1=yes, 0=no
save_table = 1;

u = yields_n - yfit_n;
u2 = u(2:end,:).^2;
rmse_n = sqrt(nanmean(u2))*10000*12;

ind_r = tT<datetime('2004-01-01') | (tT>=datetime('2008-01-01') & tT<datetime('2010-01-01')) ...
     | (tT>=datetime('2020-03-01') & tT<datetime('2021-03-01')); %higher liquidity index
u = yields_r - yfit_r;
u2=u(2:end,:).^2;
rmse_r_nor = sqrt(nanmean(u2(~ind_r(2:end),:)))*10000*12;
rmse_r_liq = sqrt(nanmean(u2(ind_r(2:end),:)))*10000*12;


u = surv_infexp_int - surv_infexp_fit;
u2 = u(2:end,:).^2;
rmse_s = sqrt(nanmean(u2))*10000*12;

u = surv_gdpexp_int - surv_gdpexp_fit;
u2 = u(2:end,:).^2;
rmse_g = sqrt(nanmean(u2))*10000*12;

u = surv_tbexp_int - surv_tbexp_fit;
u2 = u(2:end,:).^2;
rmse_t = sqrt(nanmean(u2))*10000*12;


% The table:
txt = {'\begin{tabularx}{\textwidth}{>{\hsize=2.0\hsize\raggedright\arraybackslash}X*{8}{>{\hsize=.875\hsize\centering\arraybackslash}X}} \hline '};
txt{2} = 'Maturities & 3m & 1y & 2y & 3y & 5y & 7y & 10y & Average \\ \hline';
txt{end+1} = '\multicolumn{9}{l}{Yields:} \\';

txt{end+1}  = 'Nominal yields ';
txt{end} = [txt{end}, sprintf('& %2.0f ', [rmse_n mean(rmse_n)]), '\\'];

txt{end+1} = 'Real yields - high liq.';
for j = 1:length(rmse_n)
    if ismember(mats_n(j), mats_r)
        txt{end} = [txt{end}, sprintf('& %2.0f ', rmse_r_nor(mats_r==mats_n(j)))];
    else
        txt{end} = [txt{end}, '& '];
    end
end
txt{end} = [txt{end}, sprintf('& %2.0f ', mean(rmse_r_nor)), '\\'];

txt{end+1} = 'Real yields - low liq.';
for j = 1:length(rmse_n)
    if ismember(mats_n(j), mats_r)
        txt{end} = [txt{end}, sprintf('& %2.0f ', rmse_r_liq(mats_r==mats_n(j)))];
    else
        txt{end} = [txt{end}, '& '];
    end
end
txt{end} = [txt{end}, sprintf('& %2.0f ', mean(rmse_r_liq)), '\\'];


txt{end} = [txt{end}, ' \hline'];
txt{end+1} = '\multicolumn{9}{l}{Surveys:} \\';


txt{end+1} = 'Inflation rate';
for j = 1:length(rmse_n)
    if ismember(mats_n(j), hstep_s)
        txt{end} = [txt{end}, sprintf('& %2.0f ', rmse_s(hstep_s==mats_n(j)))];
    else
        txt{end} = [txt{end}, '& '];
    end
end
txt{end} = [txt{end}, sprintf('& %2.0f ', mean(rmse_s)), '\\'];
%{
txt{end+1} = 'GDP growth rate';
for j = 1:length(rmse_n)
    if ismember(mats_n(j), hstep_g)
        txt{end} = [txt{end}, sprintf('& %2.2f ', rmse_g(hstep_g==mats_n(j)))];
    else
        txt{end} = [txt{end}, '& '];
    end
end
txt{end} = [txt{end}, sprintf('& %2.2f ', mean(rmse_g)), '\\'];
%}
txt{end+1} = '3m T-bill rate';
for j = 1:length(rmse_n)
    if ismember(mats_n(j), hstep_t)
        txt{end} = [txt{end}, sprintf('& %2.0f ', rmse_t(hstep_t==mats_n(j)))];
    else
        txt{end} = [txt{end}, '& '];
    end
end
txt{end} = [txt{end}, sprintf('& %2.0f ', mean(rmse_t)), '\\'];



txt{end+1} = '\hline';
txt{end+1} = '\end{tabularx}';

latexTable = sprintf('%s\n', txt{:});

if save_table
    latexFileName = '../Tables/table_rmse_revision.tex';
    fid = fopen(latexFileName, 'w');
    fprintf(fid, '%s', latexTable);
    fclose(fid);
    disp(['LaTeX table of descriptive statistics saved as ' latexFileName]);
end

% edit ../Tables/table_rmse_revision.tex
