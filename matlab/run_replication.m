%% Replicate the U.S. results in Golinski, Guilloux-Nefussi, and Renne
%
% This script evaluates the model at the parameter values used in the paper
% and reproduces the main tables and figures. It does not run the manual
% estimation/search loop that was used during development.

clearvars;
close all;

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


% Fill missing low-frequency macro and survey observations by carrying the
% latest available value forward. This mirrors the construction used in the
% paper input and in the R data-export script.
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

% GDP-growth survey expectations are not used in the reported specification.
surv_gdpexp_int = surv_gdpexp_int*nan;

% Parameter vector at the optimum used for the paper tables and figures.
% The signs below are adjusted to match the internal parametrization used by
% GGNR_model_v13.
x0 = [0.0104692704306508302,1.5,0.041666666666666699,0.000348014805880652851,0.790118757887699585,0.00221100170331838637,0.000507324858229496662,0.746945578196417714,0.946370617132632153,0.0041925040025360032,0.994999999999999996,0.00335487686007240882,3.08501560552246011e-05,9.6346703182495661e-19,0.994999999999999996,0.000163435759318986129,1,-0.389849461534773623,0.00782197986082663888,9.09179714100563069e-06,0.00246735490060132978,0.000769389970434012326,0,0.928560426437546993,0,0.297461066916042682,-0.093275823925008633,-2.24999997705116472e-05,-0.00327791599019908175,0.005,4.98672879444940427e-05,4.98672879444940427e-05,0.00100000000000000024,0.000200000000000000064,0.00100000000000000024,5.00000000000000227e-05,0.000269099350877976217,0.000202413462872257804,0.000784943899561216609,0.0774982223987592533,0.212811141517988744,-1.08506490685090284,0.192418924815925152,2.04484670902519383,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.153621200043389194,-0.22725553893644701,0.191376414590706645,-0.517871938730337433,-0.54003939127517242,0,0,0];
x0(27:29) = -x0(27:29);
x0(64:end) = -x0(64:end);
par_ind_est = [1 [] 5 [] 7:11 [] 13:15 [] 18:22 [] 24 [] 27:29 30 37:39 40:44 64:68];
[fval, estimates, x_upd, x_std, macro_fit, yfit_n, yfit_r, surv_infexp_fit, surv_gdpexp_fit, surv_tbexp_fit, term_prem_n, term_prem_r, irp] = ...
    GGNR_model_v13([macro_int], yields_n, mats_n, yields_r, mats_r, surv_infexp_int, hstep_s,surv_gdpexp_int,hstep_g,surv_tbexp_int,hstep_t, tT, x0, par_ind_est, 'eval');

par_ind_std = [1 [] 5 [] 7:10 [] 13 [] 18:22 [] 24 [] 27:29 [] 37:39 40:44 64:68];

[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, irp, fval_t, stderr_all, Hess] = ...
    GGNR_model_v13([macro_int], yields_n, mats_n, yields_r, mats_r, surv_infexp_int, hstep_s,surv_gdpexp_int,hstep_g,surv_tbexp_int,hstep_t, tT, x0, par_ind_std, 'eval');


[stderr_g, stderr_h, stderr_w, CovMat_g, CovMat_h, CovMat_w, J_ij] = ...
    stderr_emp(macro_int, yields_n, mats_n, yields_r, mats_r, surv_infexp_int, ...
    hstep_s,surv_gdpexp_int,hstep_g,surv_tbexp_int,hstep_t, tT, x0, par_ind_est, par_ind_std);

% Standard errors used in the paper tables. Keeping these fixed avoids small
% numerical differences across Matlab versions in the finite-difference step.
stderr_g = [0.00422375673990622	0.389298510930572	0.000102394041406234	0.977257594985957	0.0144546549475322	0.000360256414274548	3.53201334278942e-06	1.85839312001177	0.00684118331246009	0.000108630756841868	0.000391713007506593	0.000520508426802116	0.0815567000621122	0.107467733286203	0.434413354203170	0.430253379170576	1.95665014939236e-05	5.36625933592381e-05	0.000312235677405370	0.527666178284888	0.232993301648023	15.5931279754173	4.03417166233649	265.702235421144	0.264222144302864	2.61461067535409	3.73902640543924	1.51204245649152	30.9733095800693];
stderr_w = [0.000281208948195397	0.00798532582274681	2.36732957306201e-06	0.0275153015985763	0.00194486379352145	8.17029203149404e-05	8.93103417927143e-07	0.0141385547736520	0.000223694152248699	1.27664276545052e-05	1.23627938185219e-05	3.29050066446198e-05	0.00224800712410666	0.000865936526552891	0.000643116145854824	0.000289515818334756	1.07966833211745e-06	3.35862879405972e-06	1.26604624526543e-05	0.00117478041208538	0.00742054490393229	0.0226495887708648	0.0191885993396253	0.196803953748753	0.00561013890255147	0.0156627223044626	0.0273215968141118	0.0164129186243947	0.532895750905238];

stderr = nan(length(x0),1);
stderr(par_ind_std) = stderr_g;

make_table_pars;
make_table_rmse;

% Main paper figures.
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

% Additional impulse-response figures for selected initial conditions.
figure_irf_yields(estimates, [24 120], 120, 0.0025/x0(4), [zeros(1,2); 0.04 0 ], {'s_t', 'r^*_t+\pi^*_t'});
print(gcf, "-dpdf", "-painters", "../Figures/fig_irf_s.pdf")
figure_irf_yields(estimates, [24 120], 120, -0.0025/x0(4), [0.0025*ones(1,2); 0.04 0 ], {'s_t', 'r^*_t+\pi^*_t'});
print(gcf, "-dpdf", "-painters", "../Figures/fig_irf_s=25bp,r_pi.pdf")
figure_irf_yields(estimates, [24 120],  120, 0.0025/x0(4), [ones(1,2)*0.0; 0.02 -0.02], {'s_t', 'r^*_t'});
print(gcf, "-dpdf", "-painters", "../Figures/fig_irf_s=0,r.pdf")
figure_irf_yields(estimates, [24 120],  120, 0.0025/x0(4), [0.04 -0.02], {'s_t'});
print(gcf, "-dpdf", "-painters", "../Figures/fig_irf_s=0bp,r_pi.pdf")

disp("Replication complete. Outputs were written to ../Figures and ../Tables.")





