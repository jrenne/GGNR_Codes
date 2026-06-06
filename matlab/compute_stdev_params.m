% =========================================================================
% Compute standard deviation of parameter estimates
% =========================================================================

% Choose method to calculate Covariance matrix of sub_parameters
% 1 = Matlab's fminunc Hessian
% 2 = self-calculated Hessian
% 3 = Outer product of gradients
% 4 = Sandwich estimator (of methods 1 and 3)
cov_method = 3; 

if indic_estim_MLE == 1
    epsilon = 0.01; %10^(-10); %epsilon = 10^-5;

    % avoid recalculation of endogenous parameters
    endo_flag = 0;

    % Create vectors of all and of optimised parameters:
    FILTER = FILTER_MLE;
    parameters     = model_sol.param;
    sub_parameters = model_sol.param(FILTER==1);

    if cov_method == 1
    
        Hessian = model_sol.Hessian;
        CovMat = Hessian^(-1);

    elseif cov_method == 2

        f = @(x)compute_logl_aux(x,Data_StateSpace,model_sol,endo_flag);
        
        % Check log-likelihood:
        f0 = f(sub_parameters);
        disp(f0);
        
        fi = nan(sum(FILTER),1);
        for i=1:sum(FILTER)
            sub_parameters_i    = sub_parameters;
            sub_parameters_i(i) = sub_parameters_i(i) + epsilon;
            fi(i) = f(sub_parameters_i);
        end
        
        J = zeros(sum(FILTER),sum(FILTER));
        for i=1:sum(FILTER)
            sub_parameters_i    = sub_parameters;
            sub_parameters_i(i) = sub_parameters_i(i) + epsilon;
            %fi = f(sub_parameters_i);
            for j=1:i
                sub_parameters_ij    = sub_parameters_i;
                sub_parameters_ij(j) = sub_parameters_ij(j) + epsilon;
                fij = f(sub_parameters_ij);
                J(i,j) = (fij - fi(i) - fi(j) + f0)/epsilon^2;
                %J(i,j) = (fij - 2*fi + f0)/epsilon^2;
                J(j,i) = J(i,j);
            end
        end
        
        CovMat = J^(-1);

    elseif cov_method == 3
      
        % Check log-likelihood:
        [~,logl_t] = compute_logl_aux(sub_parameters,Data_StateSpace,model_sol,endo_flag);
        disp(sum(logl_t));
        
        T = length(Data_StateSpace.dataset);
        grad_lokl_t = nan(T,sum(FILTER));
        
        for i = 1:sum(FILTER)
            sub_params_i     = sub_parameters;
            sub_params_i(i)  = sub_params_i(i) + epsilon;
            [~,logl_t_i]     = compute_logl_aux(sub_params_i,Data_StateSpace,model_sol,endo_flag);
            grad_lokl_t(:,i) = (logl_t_i - logl_t)/epsilon;
        end
        
        OPG = grad_lokl_t'*grad_lokl_t;
        CovMat = OPG^(-1);

    elseif cov_method == 4

        Hessian = model_sol.Hessian;

        % Check log-likelihood:
        [~,logl_t] = compute_logl_aux(sub_parameters,Data_StateSpace,model_sol,endo_flag);
        disp(sum(logl_t));
        
        T = length(Data_StateSpace.dataset);
        grad_lokl_t = nan(T,sum(FILTER));
        
        for i = 1:sum(FILTER)
            sub_params_i     = sub_parameters;
            sub_params_i(i)  = sub_params_i(i) + epsilon;
            [~,logl_t_i]     = compute_logl_aux(sub_params_i,Data_StateSpace,model_sol,endo_flag);
            grad_lokl_t(:,i) = (logl_t_i - logl_t)/epsilon;
        end
        
        OPG = grad_lokl_t' * grad_lokl_t;
        CovMat = Hessian^(-1) * OPG * Hessian^(-1);

    else

        disp('Choose method for covariance calculation.')
        return;

    end

    % Matrix G will contain the derivatives of transformed param w.r.t. truely estimated parameters:
    param_transf_baseline = model_sol.param_transf';   
    G = zeros(max(size(param_transf_baseline)),size(CovMat,1));

    % % Numerical derivatives for delta method
    % for i=1:sum(FILTER)
    %     sub_parameters_i        = sub_parameters;
    %     sub_parameters_i(i)     = sub_parameters_i(i) + epsilon;
    %     parameters_i            = parameters;
    %     parameters_i(FILTER==1) = sub_parameters_i;
    % 
    %     [rho_g, rho_z, rho_w, rho_m, rho_k,...
    %         sigma_g, sigma_z, sigma_w, sigma_m, sigma_k,...
    %         rho_pi_star, sigma_pi_star, rho_pi_tilde, sigma_pi_z, mu_pi, mu_c,...
    %         mu_gamma, mu_kappa, delta, rho_gz] = make_param_transf(parameters_i);
    % 
    %     param_transf_temp = [rho_g, rho_z, rho_w, rho_m, rho_k,...
    %         sigma_g, sigma_z, sigma_w, sigma_m, sigma_k,...
    %         rho_pi_star, sigma_pi_star, rho_pi_tilde, sigma_pi_z, mu_pi, mu_c,...
    %         mu_gamma, mu_kappa, delta, rho_gz]';
    % 
    %     % reload all other parameters to avoid numerical imprecisions
    %     indic_param = find(FILTER);
    %     param_transf_i    = param_transf_baseline;
    %     param_transf_i(indic_param(i)) = param_transf_temp(indic_param(i));
    % 
    %     G(:,i) = (param_transf_i - param_transf_baseline)/epsilon;
    % end

    % Analytical derivatives for delta method
    for i=[1:5,11,13,19]
        G(i,sum(FILTER(1:i))) = exp(parameters(i))/(1+exp(parameters(i)))^2;
    end
    for i=[6:10,12,14:15,20] %16 and 17 not included!
        G(i,sum(FILTER(1:i))) = exp(parameters(i));
    end
    for i=18
        G(i,sum(FILTER(1:i))) = 2*exp(parameters(i))/(1+exp(parameters(i)))^2;
    end

    % Get covariance matrix of transformed economic parameters
    CovMat_transf  = G * CovMat * G';
    
    % Obtain StDev for mu_c & mu_gamma via Monte Carlo draws of sub_parameters

    rng(1337);
    n_draws          = 10000;
    sub_param_draws  = mvnrnd(sub_parameters,CovMat,n_draws);
    param_transf_MC  = nan(n_draws,length(parameters));

    % ensure recalculation of endogenous parameters
    endo_flag = 1;
    
    % Loop through each draw
    for d = 1:n_draws
        disp(d);
        sub_parameters_d = sub_param_draws(d,:)';
        model_d = model_sol;
        model_d.param(FILTER==1) = sub_parameters_d;
        model_d = make_model_sol(model_d,endo_flag);
        param_transf_MC(d,:) = model_d.param_transf;
    end

    % Get covariance matrix based on MC draws
    CovMat_transf_MC = cov(param_transf_MC);
    
    % Compare original correlation matrix and MC-based matrix
    % For sub_parameters, the max difference should be small (otherwise increase n_draws)
    CorrMat_transf    = corrcov(CovMat_transf);
    CorrMat_transf_MC = corrcov(CovMat_transf_MC);
    max_error = max(max(abs(CorrMat_transf_MC(FILTER==1,FILTER==1) - CorrMat_transf(FILTER==1,FILTER==1))));
    disp(max_error);

    % Plot histogram for mu_c and mu_gamma
    figure('WindowState','maximized');
    subplot(1,2,1);
    indic_param = strcmp(model_sol.names_param,'mu_c');
    histogram(param_transf_MC(:,indic_param),'Normalization','pdf');
    hold on;
    line([mean(param_transf_MC(:,indic_param)) mean(param_transf_MC(:,indic_param))],ylim, 'Color', 'r', 'LineWidth', 2);
    line([median(param_transf_MC(:,indic_param)) median(param_transf_MC(:,indic_param))],ylim, 'Color', 'b', 'LineWidth', 2);
    line([model_sol.mu_c0 model_sol.mu_c0],ylim, 'Color', 'k', 'LineWidth', 2);
    hold off;
    legend('MC histrogram','MC mean','MC median','point estimate')
    title('pdf of mu\_c');
    subplot(1,2,2);
    indic_param = strcmp(model_sol.names_param,'mu_gamma');
    histogram(param_transf_MC(:,indic_param),'Normalization','pdf');
    hold on;
    line([mean(param_transf_MC(:,indic_param)) mean(param_transf_MC(:,indic_param))],ylim, 'Color', 'r', 'LineWidth', 2);
    line([median(param_transf_MC(:,indic_param)) median(param_transf_MC(:,indic_param))],ylim, 'Color', 'b', 'LineWidth', 2);
    line([model_sol.mu_gamma0 model_sol.mu_gamma0],ylim, 'Color', 'k', 'LineWidth', 2);
    hold off;
    legend('MC histrogram','MC mean','MC median','point estimate')
    title('pdf of mu\_gamma');

    % Save results
    model_sol.CovMat           = CovMat;
    model_sol.CovMat_transf    = CovMat_transf;
    model_sol.CovMat_transf_MC = CovMat_transf_MC;

    if indic_save_model == 1
        save("results/save_MLE_approach.mat","model_sol");
    end
else
    param_transf_baseline = model_sol.param_transf';
    CovMat                = model_sol.CovMat;
    CovMat_transf         = model_sol.CovMat_transf;
    CovMat_transf_MC      = model_sol.CovMat_transf_MC;
end

param_transf_stdev = sqrt(diag(CovMat_transf));

indic_param = strcmp(model_sol.names_param,'mu_c');
param_transf_stdev(indic_param) = sqrt(CovMat_transf_MC(indic_param,indic_param));

indic_param = strcmp(model_sol.names_param,'mu_gamma');
param_transf_stdev(indic_param) = sqrt(CovMat_transf_MC(indic_param,indic_param));

param_transf_lower = param_transf_baseline - 2*param_transf_stdev;
param_transf_upper = param_transf_baseline + 2*param_transf_stdev;

table2print = array2table(round([param_transf_baseline param_transf_stdev param_transf_lower param_transf_upper],5),'RowNames',model_sol.names_param,...
    'VariableNames',{'Estimate','StDev','Est-2xStDev','Est+2xStDev'});
disp(table2print);
