    function [A_X_for, B_X_for, A_X_exp, B_X_cum, A_X_for_pi, B_X_pi, A_X_exp_pi, B_X_cum_pi, Conv_pi] = ...
            affine_coefs(Phi_q, Mu_q, Sigma2_X, delta_0, delta_1, delta_pi_0, delta_pi_1)
        
        maxmat = 120;
        freq = 1;
        K = length(Mu_q);
        B_X_for=zeros(K,maxmat); % forward rate loadings (starting from f0=y1) x freq
        B_X_for(:,1) = delta_1;
        A_X_exp = zeros(maxmat,1); % intercept for forward rates (starting from f0=y1) x freq
        A_X_exp(1) = delta_0;
        A_X_for = zeros(maxmat,1); % intercept for forward rates (starting from f0=y1) x freq
        A_X_for(1) = A_X_exp(1);
        B_X_cum = zeros(K,maxmat); % cumulative loadings
        
        for j=2:maxmat
            B_X_for(:,j) = Phi_q'*B_X_for(:,j-1);
            A_X_exp(j) = A_X_exp(j-1) + B_X_for(:,j-1)'*Mu_q;
            B_X_cum(:,j) = B_X_cum(:,j-1) + B_X_for(:,j-1);
            A_X_for(j) = A_X_exp(j) - 0.5*B_X_cum(:,j)'*Sigma2_X*B_X_cum(:,j)/freq;
        end
        
        B_X_pi = zeros(K,maxmat); % forward rate loadings (starting from f0=y1) x freq
        B_X_pi(:,1) = Phi_q'*delta_pi_1;
        A_X_exp_pi = zeros(maxmat,1); % intercept for forward rates (starting from f0=y1) x freq
        A_X_exp_pi(1) = delta_pi_0 + delta_pi_1'*Mu_q;
        B_X_cum_pi = zeros(K,maxmat); % cumulative of all B_X_for's starting from 0 
        B_X_cum_pi(:,1) = delta_pi_1; %
        Conv_pi = zeros(maxmat,1); % Convexity term for forward rates
        Conv_pi(1) = B_X_cum_pi(:,1)'*Sigma2_X*B_X_cum_pi(:,1)/freq;
        A_X_for_pi = zeros(maxmat,1); % intercept for forward rates (starting from f0=y1) x freq
        A_X_for_pi(1) = A_X_exp_pi(1) + 0.5*Conv_pi(1);
        
        for j=2:maxmat
            B_X_pi(:,j) = Phi_q'*B_X_pi(:,j-1);
            A_X_exp_pi(j) = A_X_exp_pi(j-1) + B_X_pi(:,j-1)'*Mu_q;
            B_X_cum_pi(:,j) = B_X_cum_pi(:,j-1) + B_X_pi(:,j-1);
            Conv_pi(j) = B_X_cum_pi(:,j)'*Sigma2_X*B_X_cum_pi(:,j)/freq;
            A_X_for_pi(j) = A_X_exp_pi(j) + 0.5*Conv_pi(j);
        end
    end
