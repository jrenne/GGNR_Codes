% =========================================================================
% Create Latex table of model parameterization
% =========================================================================

% Save table? 1=yes, 0=no
save_table = 1;


pars_for_table = x0;
pars.names = cell(length(x0),1);
pars.multiplier = zeros(length(x0),1);% Annualization multiplier
pars.multiplier(:) = 1;

% Model parameters:
pars.names{1} = '\rho_i';
pars.names{2} = '\alpha_{\pi}';
pars.names{3} = '\alpha_z';
pars.multiplier(strcmp(pars.names, '\alpha_z')) = 12;
pars.names{4} = '\sigma_i';
pars.multiplier(strcmp(pars.names, '\sigma_i')) = 1200;
pars.names{5} = '\rho_g';
pars.names{6} = '\mu_g';
pars.multiplier(strcmp(pars.names, '\mu_g')) = 1200;
pars.names{7} = '\sigma_g';
pars.multiplier(strcmp(pars.names, '\sigma_g')) = 1200;
pars.names{8} = '\alpha';
pars.multiplier(strcmp(pars.names, '\alpha')) = 100;
pars.names{9} = '\rho_z';
pars.names{10} = '\sigma_z';
pars.multiplier(strcmp(pars.names, '\sigma_z')) = 100;
pars.names{11} = '\rho^*';
pars.names{12} = '\mu^*';
pars.multiplier(strcmp(pars.names, '\mu^*')) = 1200;
pars.names{13} = '\sigma^*';
pars.multiplier(strcmp(pars.names, '\sigma^*')) = 1200;
pars.names{14} = '\mu_{\kappa}';
pars.multiplier(strcmp(pars.names, '\mu_{\kappa}')) = 1200;
pars.names{15} = '\rho_{\kappa}';
pars.names{16} = '\sigma_{\kappa}';
pars.multiplier(strcmp(pars.names, '\sigma_{\kappa}')) = 1200;
pars.names{17} = '\theta';
pars.names{18} = '\bar{\rho}';
pars.names{19} = '\beta';
pars.multiplier(strcmp(pars.names, '\beta')) = 12;
pars.names{20} = '\sigma_{\pi}';
pars.multiplier(strcmp(pars.names, '\sigma_{\pi}')) = 1200;
pars.names{21} = 'a_{u,1}';
pars.multiplier(strcmp(pars.names, '\a_{u,1}')) = 1200;
pars.names{22} = 'a_{u,2}';
pars.multiplier(strcmp(pars.names, '\a_{u,2}')) = 1200;
pars.names{23} = '\underline{i}';
pars.names{24} = '\rho_w';
pars.names{25} = '\mu_w';
pars.multiplier(strcmp(pars.names, '\mu_w')) = 1200;
pars.names{26} = '\sigma_w';
pars.names{27} = '\sigma_{w,g}';
pars.names{28} = '\sigma_{w,z}';
pars.names{29} = '\sigma_{w,\pi}';

% Measurement errors:
pars.names{30} = '\sigma_o';
pars.multiplier(strcmp(pars.names, '\sigma_o')) = 100;
pars.names{31} = '\sigma_{gdp}';
pars.multiplier(strcmp(pars.names, '\sigma_{gdp}')) = 1200;
pars.names{32} = '\sigma_{inf}';
pars.multiplier(strcmp(pars.names, '\sigma_{inf}')) = 1200;
pars.names{33} = '\sigma_{ptr}';
pars.multiplier(strcmp(pars.names, '\sigma_{ptr}')) = 1200;
pars.names{34} = '\sigma_s';
pars.multiplier(strcmp(pars.names, '\sigma_s')) = 1200;
pars.names{35} = '\sigma_y';
pars.multiplier(strcmp(pars.names, '\sigma_y')) = 1200;
pars.names{36} = '\sigma_t';
pars.multiplier(strcmp(pars.names, '\sigma_t')) = 1200;
pars.names{37} = '\sigma_n';
pars.multiplier(strcmp(pars.names, '\sigma_n')) = 1200;
pars.names{38} = '\sigma_r';
pars.multiplier(strcmp(pars.names, '\sigma_r')) = 1200;
pars.names{39} = '\sigma_{r,liq}';
pars.multiplier(strcmp(pars.names, '\sigma_{r,liq}')) = 1200;

% Price of risk parameters:
pars.names{40} = '\lambda_{0,\varepsilon_i}';
pars.names{41} = '\lambda_{0,\varepsilon_g}';
pars.names{42} = '\lambda_{0,\varepsilon_z}';
pars.names{43} = '\lambda_{0,\varepsilon^*}';
pars.names{44} = '\lambda_{0,\varepsilon_\pi}';
pars.names{45} = '\lambda_{0,6}';
pars.names{46} = '\lambda_{0,7}';
pars.names{47} = '\lambda_{0,8}';
pars.names{48} = '\lambda_{z,1}';
pars.names{49} = '\lambda_{z,2}';
pars.names{50} = '\lambda_{z,3}';
pars.names{51} = '\lambda_{z,4}';
pars.names{52} = '\lambda_{z,5}';
pars.names{53} = '\lambda_{z,6}';
pars.names{54} = '\lambda_{z,7}';
pars.names{55} = '\lambda_{z,8}';
pars.names{56} = '\lambda_{\pi,1}';
pars.names{57} = '\lambda_{\pi,2}';
pars.names{58} = '\lambda_{\pi,3}';
pars.names{59} = '\lambda_{\pi,4}';
pars.names{60} = '\lambda_{\pi,5}';
pars.names{61} = '\lambda_{\pi,6}';
pars.names{62} = '\lambda_{\pi,7}';
pars.names{63} = '\lambda_{\pi,8}';
pars.names{64} = '\lambda_{w,\varepsilon_i}';
pars.names{65} = '\lambda_{w,\varepsilon_g}';
pars.names{66} = '\lambda_{w,\varepsilon_z}';
pars.names{67} = '\lambda_{w,\varepsilon^*}';
pars.names{68} = '\lambda_{w,\varepsilon_\pi}';
pars.names{69} = '\lambda_{w,6}';
pars.names{70} = '\lambda_{w,7}';
pars.names{71} = '\lambda_{w,8}';


%par_ind_est1 = par_ind_est(1:22); % estimated model parameters
par_ind_est1 = par_ind_est(1:19); % estimated model parameters
par_ind_est2 = [par_ind_est(20:23) 34]; % measurement errors
par_ind_est_aux = par_ind_est(24:end); % price of risk parameters


% The table:
txt = {'\begin{tabular}{crc|crrc} \hline '};
txt{2} = 'Parameter & Estimate (S.E.) & Multiplier';
txt{2} = [txt{2}, ' & ', txt{2}, '\\ \hline'];
txt{end+1} = '\multicolumn{4}{l}{Estimated model parameters:} & & \\';
% The Table comes here
txt_line = [];
for j = 1:length(par_ind_est1)
    par_txt = [pars.names{par_ind_est1(j)}];
    est_val = x0(par_ind_est1(j))*pars.multiplier(par_ind_est1(j));
    if ~exist('stderr','var')  || isnan(stderr(par_ind_est1(j)))
        stderr_txt = '-';
    else
        stderr_txt = num2str(stderr(par_ind_est1(j))*pars.multiplier(par_ind_est1(j)), '%.4f');
    end
    txt_line = [txt_line, '& $ ', par_txt, ' $ & $ \underset{', stderr_txt, '}{', num2str(est_val, '%.4f'), ...
        '} $ & $ \times ', num2str(pars.multiplier(par_ind_est1(j))) ' $ '];

    if mod(j,2)==0
        txt{end+1} = [txt_line(2:end), ' \\'];
        txt_line = [];
    end
end
if ~isempty(txt_line)
    txt{end+1} = [txt_line(2:end), ' \\'];
    txt_line = [];
end

txt{end+1} = '\multicolumn{4}{l}{Calibrated model parameters:} & & \\';
par_ind_est1 = setdiff(1:par_ind_est1(end), par_ind_est1);
for j = 1:length(par_ind_est1)
    par_txt = [pars.names{par_ind_est1(j)}];
    est_val = x0(par_ind_est1(j))*pars.multiplier(par_ind_est1(j));
    txt_line = [txt_line, '& $ ', par_txt, ' $ & $ ' , num2str(est_val, '%.4f'), ...
        ' $ & $ \times ', num2str(pars.multiplier(par_ind_est1(j))) ' $ '];

    if mod(j,2)==0
        txt{end+1} = [txt_line(2:end), ' \\'];
        txt_line = [];
    end
end
if ~isempty(txt_line)
    txt{end+1} = [txt_line(2:end), ' \\'];
    txt_line = [];
end

txt{end+1} = '\multicolumn{4}{l}{Measurement errors:} & & \\';
for j = 1:length(par_ind_est2)
    par_txt = [pars.names{par_ind_est2(j)}];
    est_val = x0(par_ind_est2(j))*pars.multiplier(par_ind_est2(j));
    if ~exist('stderr','var')  || isnan(stderr(par_ind_est2(j)))
        stderr_txt = '-';
    else
        stderr_txt = num2str(stderr(par_ind_est2(j))*pars.multiplier(par_ind_est2(j)), '%.4f');
    end
    txt_line = [txt_line, '& $ ', par_txt, ' $ & $ \underset{', stderr_txt, '}{', num2str(est_val, '%.4f'), ...
        '} $ & $ \times ', num2str(pars.multiplier(par_ind_est2(j))) ' $ '];

    if mod(j,2)==0
        txt{end+1} = [txt_line(2:end), ' \\'];
        txt_line = [];
    end
end
txt{end+1} = [txt_line(2:end), '& \multicolumn{2}{l}{$ \qquad \sigma_{ptr} = 5 \times \sigma_s $} & \\ '];

txt{end+1} = [' \multicolumn{2}{l}{$ \qquad \sigma_{gdp} = \sigma_o/100 $} & & ' ...
    ' \multicolumn{2}{l}{$\qquad \sigma_{inf} = \sigma_o/100 $} & \\'];
txt{end+1} = [' \multicolumn{2}{l}{$ \qquad \sigma_{y} = \sigma_s $} & & ' ...
    ' \multicolumn{2}{l}{$\qquad \sigma_{3tb} = \sigma_s/4 $} & \\'];


txt{end+1} = '\multicolumn{4}{l}{Price of risk parameters:} & & \\';
par_ind_est1 = zeros(size(par_ind_est_aux));
par_ind_est1(1:2:end) = par_ind_est_aux(1:ceil(length(par_ind_est1)/2));
par_ind_est1(2:2:end) = par_ind_est_aux(length(par_ind_est1)/2+1:end);
txt_line = [];
for j = 1:length(par_ind_est1)
    par_txt = [pars.names{par_ind_est1(j)}];
    est_val = x0(par_ind_est1(j))*pars.multiplier(par_ind_est1(j));
    if ~exist('stderr','var')  || isnan(stderr(par_ind_est1(j)))
        stderr_txt = '-';
    else
        stderr_txt = num2str(stderr(par_ind_est1(j))*pars.multiplier(par_ind_est1(j)), '%.4f');
    end
    txt_line = [txt_line, '& $ ', par_txt, ' $ & $ \underset{', stderr_txt, '}{', num2str(est_val, '%.4f'), ...
        '} $ &  '];

    if mod(j,2)==0
        txt{end+1} = [txt_line(2:end), ' \\'];
        txt_line = [];
    end
end


txt{end+1} = '\hline';
txt{end+1} = '\end{tabular}';

latexTable = sprintf('%s\n', txt{:});

if save_table
    latexFileName = '../Tables/table_param_revision.tex';
    fid = fopen(latexFileName, 'w');
    fprintf(fid, '%s', latexTable);
    fclose(fid);
    disp(['LaTeX table of descriptive statistics saved as ' latexFileName]);
end
% edit ../Tables/table_param_revision.tex

