%% Replicate the results in Golinski, Guilloux-Nefussi, and Renne
%
% This script evaluates the model at the parameter values used in the paper.

clearvars;
close all;

run_estimation = false; % set to true to continue estimation from x0
max_estimation_iter = 100; % used only when run_estimation is true
 
script_dir = fileparts(mfilename('fullpath'));
cd(script_dir);
addpath(script_dir);

if ~exist('../Figures', 'dir')
    mkdir('../Figures');
end
if ~exist('../Tables', 'dir')
    mkdir('../Tables');
end

load_replication_data('../data/data_JPR.mat')


% Interpolate missing observations.
macro_int = macro;
for j = 3:length(tT)
    if isnan(macro_int(j,3))
        macro_int(j,3:4) = macro_int(j-1,3:4);
    end
end


surv_infexp_int = surv_infexp;
for j = 155:length(tT)
    if isnan(surv_infexp_int(j,1))
        surv_infexp_int(j,1) = surv_infexp_int(j-1,1);
    end
end
for j = 278:length(tT)
    if isnan(surv_infexp_int(j,2))
        surv_infexp_int(j,2) = surv_infexp_int(j-1,2);
    end
end

surv_gdpexp_int = surv_gdpexp;
for j = 2:length(tT)
    if isnan(surv_gdpexp_int(j,1))
        surv_gdpexp_int(j,1) = surv_gdpexp_int(j-1,1);
    end
end
for j = 281:length(tT)
    if isnan(surv_gdpexp_int(j,2))
        surv_gdpexp_int(j,2) = surv_gdpexp_int(j-1,2);
    end
end


surv_tbexp_int = surv_tbexp;
for j = 155:length(tT)
    if isnan(surv_tbexp_int(j,1))
        surv_tbexp_int(j,1) = surv_tbexp_int(j-1,1);
    end
end
for j = 281:length(tT)
    if isnan(surv_tbexp_int(j,2))
        surv_tbexp_int(j,2) = surv_tbexp_int(j-1,2);
    end
end


%% Evaluate the paper specification

surv_gdpexp_int = surv_gdpexp_int*nan;

% Parameter vector at the optimum.
x0 = [0.0104382545811124005,1.5,0.041666666666666699,0.000424503825792636974,0.79041234946036798,0.00221100170331838984,0.000507469265124473024,0.747219154808188013,0.946493694514137007,0.00417761772895170001,0.994999999999999996,0.00335487686007241012,3.07978178192533006e-05,1.14215241622012997e-19,0.994999999999999996,0.000236801488497918004,1,-0.381562308802388006,0.00780528812560161035,1.10528561709930992e-05,0.00245726191949591785,0.00082486187531531543,0,0.928721447221876995,0,0.370778199831249988,0.0933049551681942035,0.000918623576709006969,0.00328229760763929001,0.00499451540382658994,4.99451540382658975e-05,4.99451540382658975e-05,0.00100000000000000002,0.00020000000000000001,0.00100000000000000002,5.00000000000000024e-05,0.000252040803131816022,0.000225037669504208997,0.000781977814498796098,0.0776018851782361019,0.21309572617054201,-1.08650492172258994,0.202307500330795004,2.1269596752413702,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.150750599551518011,0.227559426649380003,-0.43899629971773102,0.544441562310157945,3.28184089816334534,0,0,0];
par_ind_est = [1 [] 5 [] 7:11 [] 13:15 [] 18:22 [] 24 [] 27:29 30 37:39 40:44 64:68];
if run_estimation
    model_mode = 'estimate';
else
    model_mode = 'eval';
end
[fval, estimates, x_upd, x_std, macro_fit, yfit_n, yfit_r, surv_infexp_fit, surv_gdpexp_fit, surv_tbexp_fit, term_prem_n, term_prem_r, irp] = ...
    GGNR_model([macro_int], yields_n, mats_n, yields_r, mats_r, surv_infexp_int, hstep_s,surv_gdpexp_int,hstep_g,surv_tbexp_int,hstep_t, tT, x0, par_ind_est, model_mode, max_estimation_iter);

par_ind_std = [1 [] 5 [] 7:10 [] 13 [] 18:22 [] 24 [] 27:29 [] 37:39 40:44 64:68];

[stderr_g, ~, ~] = stderr_emp(macro_int, yields_n, mats_n, yields_r, mats_r, surv_infexp_int, ...
    hstep_s,surv_gdpexp_int,hstep_g,surv_tbexp_int,hstep_t, tT, x0, par_ind_est, par_ind_std);

stderr = nan(length(x0),1);
stderr(par_ind_std) = stderr_g;

x_bounds = x_bounds_vals();

make_table_pars;
make_table_rmse;

figure_term_prem;
figure_yields_nom_fit;
figure_yields_real_fit;
figure_surv_fit;
figure_stars_vars
figure_macro_vars;
figure_shadow_rate;

figure_irf_macro(estimates, 120);
figure_prob_contourplot(x0, [{-2:0.1:6}, {1:0.1:6}], {'r^*_t', '\pi^*_t'});
cond_corr_simulation(x0, x_upd, 10000, 12, [3 12 120], tT);

figure_irf_yields(estimates, [24 120], 120, 0.0025/x0(4), [zeros(1,2); 0.04 0 ], {'s_t', 'r^*_t+\pi^*_t'});
save_pdf_figure(gcf, "../Figures/fig_irf_s=0bp,r_pi.pdf")
figure_irf_yields(estimates, [24 120], 120, -0.0025/x0(4), [0.0025*ones(1,2); 0.04 0 ], {'s_t', 'r^*_t+\pi^*_t'});
save_pdf_figure(gcf, "../Figures/fig_irf_s=25bp,r_pi.pdf")
figure_irf_yields(estimates, [24 120],  120, 0.0025/x0(4), [ones(1,2)*0.0; 0.02 -0.02], {'s_t', 'r^*_t'});
save_pdf_figure(gcf, "../Figures/fig_irf_s=0,r.pdf")
figure_irf_yields(estimates, [24 120],  120, 0.0025/x0(4), [0.04 -0.02], {'s_t'});
save_pdf_figure(gcf, "../Figures/fig_irf_s.pdf")

disp("Replication complete. Outputs were written to ../Figures and ../Tables.")
