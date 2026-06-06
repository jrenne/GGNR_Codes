function [stderr_g, stderr_h, stderr_w, CovMat_g_x0, CovMat_h_x0, CovMat_w_x0, J_ij] = ...
    stderr_emp(macro_int, yields_n, mats_n, yields_r, mats_r, surv_infexp_int, ...
    hstep_s,surv_gdpexp_int,hstep_g,surv_tbexp_int,hstep_t, tT, x0, par_ind_est, par_ind_std)


[logl0, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, logl_t0] = ...
    GGNR_model_v13([macro_int], yields_n, mats_n, yields_r, mats_r, surv_infexp_int, ...
    hstep_s,surv_gdpexp_int,hstep_g,surv_tbexp_int,hstep_t, tT, x0, par_ind_est, 'eval');

x_bounds = x_bounds_vals;
x1 = pars_trans_in(x0, x_bounds);
[~, deltas] = pars_trans_out(x1, x_bounds);

grad_t = nan(size(logl_t0,1), length(par_ind_std)); % t-time gradient
logl_j = nan(length(par_ind_std),1);
J_ij = nan(length(par_ind_std)); % Hessian

%g_j = nan(length(par_ind_std),1);
%G_ij = nan(length(par_ind_std)); % sum log-likelihood gradient


for j = par_ind_std
    disp(find(par_ind_std==j))

    if ismember(j, [66])
        par_range = 0:0.2*x1(j):2*x1(j);
    elseif ismember(j, [8, 22, 41])
        par_range = 0.5*x1(j):0.1*x1(j):1.5*x1(j);
    elseif ismember(j, [7, 18, 20, 27, 29, 43, 44, 65, 67])
        par_range = 0.8*x1(j):0.04*x1(j):1.2*x1(j);
    elseif ismember(j, [28, 68])
        par_range = x1(j)-5*0.01:0.01:x1(j)+5*0.01;
    elseif abs(x1(j)) >= 0.001
        par_range = 0.9*x1(j):0.02*x1(j):1.1*x1(j);
    else
        par_range = x1(j)-5*0.0002:0.0002:x1(j)+5*0.0002;
    end

    node_distance = 3;
    epsilon_j = abs(par_range(6+node_distance) - par_range(6));
    
    x1_aux_j = x1;
    x1_aux_j(j) = par_range(6+node_distance);
    x0_aux_j = x0;
    x0_aux_j(j) = pars_trans_out(par_range(6-node_distance), x_bounds(j,:));
    [logl_j1, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, logl_tj1] = ...
        GGNR_model_v13([macro_int], yields_n, mats_n, yields_r, mats_r, surv_infexp_int, ...
        hstep_s,surv_gdpexp_int,hstep_g,surv_tbexp_int,hstep_t, tT, x0_aux_j, par_ind_est, 'eval');
    
    x0_aux_j(j) = pars_trans_out(par_range(6+node_distance), x_bounds(j,:));
    [logl_j2, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, logl_tj2] = ...
        GGNR_model_v13([macro_int], yields_n, mats_n, yields_r, mats_r, surv_infexp_int, ...
        hstep_s,surv_gdpexp_int,hstep_g,surv_tbexp_int,hstep_t, tT, x0_aux_j, par_ind_est, 'eval');
    
    logl_tj = sum(logl_tj1+logl_tj2, 2)/2;
    logl_j(j) = logl_j2;

    grad_t(:,find(par_ind_std==j)) = (logl_t0 - logl_tj)/epsilon_j;
    %g_j(find(par_ind_std==j)) = sum(logl_t0 - logl_tj2)/epsilon_j;
%{
    for i = par_ind_std(1:find(par_ind_std==j))
            if ismember(i, [66])
                par_range = 0:0.2*x1_aux_j(i):2*x1_aux_j(i);
            elseif ismember(i, [8, 22, 41])
                par_range = 0.5*x1_aux_j(i):0.1*x1_aux_j(i):1.5*x1_aux_j(i);
            elseif ismember(i, [5, 7, 18, 20, 27, 29, 43, 44, 65, 67])
                par_range = 0.8*x1_aux_j(i):0.04*x1_aux_j(i):1.2*x1_aux_j(i);
            elseif ismember(i, [28, 68])
                par_range = x1_aux_j(i)-5*0.01:0.01:x1_aux_j(i)+5*0.01;
            elseif abs(x1(i)) >= 0.001
                par_range = 0.9*x1_aux_j(i):0.02*x1_aux_j(i):1.1*x1_aux_j(i);
            else
                par_range = x1_aux_j(i)-5*0.0002:0.0002:x1_aux_j(i)+5*0.0002;
            end
        
            epsilon_i = abs(par_range(6+node_distance) - par_range(6));

            x0_aux_i = x0_aux_j;            
            x0_aux_i(i) = pars_trans_out(par_range(6+node_distance), x_bounds(i,:));
            [logl_ij, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, logl_tj3] = GGNR_model_v13([macro_int], yields_n, mats_n, yields_r, mats_r, surv_infexp_int, ...
                hstep_s,surv_gdpexp_int,hstep_g,surv_tbexp_int,hstep_t, tT, x0_aux_i, par_ind_est, 'eval');
        
            J_ij(find(par_ind_std==j),find(par_ind_std==i)) = (logl_j(i) + logl_j(j) - logl0 - logl_ij)/(epsilon_i*epsilon_j);
            J_ij(find(par_ind_std==i),find(par_ind_std==j)) = J_ij(find(par_ind_std==j),find(par_ind_std==i));

            %g_i =  sum(logl_tj2 - logl_tj3)/epsilon_i;
            %G_ij(find(par_ind_std==j),find(par_ind_std==i)) = (g_j(find(par_ind_std==j)) - g_i)/epsilon_i;
            %G_ij(find(par_ind_std==i),find(par_ind_std==j)) = G_ij(find(par_ind_std==j),find(par_ind_std==i));
    end
%}
end

OPG_x1 = grad_t'*grad_t;

CovMat_g_x1 = OPG_x1^(-1);
CovMat_g_x0 = diag(deltas(par_ind_std))*CovMat_g_x1*diag(deltas(par_ind_std));

CovMat_h_x1 = J_ij^(-1);
CovMat_h_x0 = diag(deltas(par_ind_std))*CovMat_h_x1*diag(deltas(par_ind_std));

CovMat_w_x0 = diag(deltas(par_ind_std))*CovMat_h_x1*OPG_x1*CovMat_h_x1*diag(deltas(par_ind_std));

stderr_g = sqrt(diag(CovMat_g_x0))';
stderr_h = sqrt(diag(CovMat_h_x0))';
stderr_w = sqrt(diag(CovMat_w_x0))';
%{
[L, D] = eig(J_ij);
D(D<0) = 0;

H_ij = L'*D*L;
CovMat_h_x1 = pinv(H_ij);
CovMat_w_x0 = diag(deltas(par_ind_std))*CovMat_h_x1*OPG_x1*CovMat_h_x1*diag(deltas(par_ind_std));
stderr_w = sqrt(diag(CovMat_w_x0))';
%}





