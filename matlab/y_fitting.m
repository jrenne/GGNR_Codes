function [yfit_all_n, yfit_all_r, JJ_n, JJ_r, JJ_nn, JJ_rr, Probs] = ...
            y_fitting(X0, A_X_for, B_X_for, A_X_exp, r_lb, s_n, A_X_for_pi, B_X_pi, Sigma2_X, B_X_cum, B_X_cum_pi)
        [K, T0] = size(X0);
        Mm = length(A_X_for);

        y_fit_short = A_X_exp(1) + X0'*B_X_for(:,1); % shadow short rate
        Probs_short = ones(T0,1);
        Probs_short(y_fit_short < r_lb) = 0;
        y_fit_short(y_fit_short < r_lb) = r_lb; % observed short rate, T x 1
        
        mu = A_X_exp(2:end)' + X0'*B_X_for(:,2:end) - r_lb;  % T x maxmat-1
        z_n = ( mu ./ s_n(1:end-1)'); % T x maxmat-1
        %Probs_long = mex_cdf(z_n); % T x maxmat-1
        Probs_long = normcdf(z_n); % T x maxmat-1
        Probs = [Probs_short Probs_long];
        pdf_z = normpdf(z_n); % T x maxmat-1
        
        
        % Nominal forward rates
        f_fit_long = r_lb + mu .* Probs_long + pdf_z .* s_n(1:end-1)' + Probs_long.*(A_X_for(2:end)'-A_X_exp(2:end)');
        %f_fit_long = r_lb + mu .* Probs_long + pdf_z .* s_n(1:end-1)' - 0.5*Probs_long.*sum(B_X_cum(:,2:end).*(Sigma2_X*B_X_cum(:,2:end)));
        f_fit = [y_fit_short f_fit_long]; % forward rates
        yfit_all_n = cumsum(f_fit, 2) ./ (1:Mm);
        
        % Expected inflation
        pi_fit = A_X_for_pi' + X0'*B_X_pi; % negative expected inflation

        % Real forward rates
        BXcumSig2BXcum_pi = [0 sum(B_X_cum(:,2:end).*(Sigma2_X*B_X_cum_pi(:,2:end))) ]; % 1 x maxmat
        r_fit = f_fit - pi_fit + Probs.*BXcumSig2BXcum_pi;
        yfit_all_r = cumsum(r_fit, 2) ./ (1:Mm);

        % First derivative of y wrt X
        if nargout > 2 
            % Jacobian
            JJ_f = zeros(K*T0,Mm);
            JJ_r = zeros(K*T0,Mm);
            if A_X_exp(1) + X0'*B_X_for(:,1) > r_lb
                JJ_f(:,1) = repmat(B_X_for(:,1), T0,1);
            end
            JJ_r(:,1) = JJ_f(:,1) - repmat(B_X_pi(:,1), T0, 1);
            if Mm>1
                %JJ_f(2:end,:) = repmat(Probs_long,1,K) .* B_X_bar(:,2:end)'; % Jacobian for forwards
                JJ_f(:,2:end) = repelem(Probs_long,K,1) .* repmat(B_X_for(:,2:end),T0,1); % T*K x maxmat
                JJ_r(:,2:end) = JJ_f(:,2:end) + (-repmat(B_X_pi(:,2:end), T0, 1) + ...
                    repelem(pdf_z./s_n(1:end-1)',K,1).*repmat(B_X_for(:,2:end),T0,1).*BXcumSig2BXcum_pi(2:end) );
            end
            JJ_n = cumsum(JJ_f,2)' ./ (1:Mm)'; %Jacobian for nominal yields
            JJ_r = cumsum(JJ_r,2)' ./ (1:Mm)'; %Jacobian for nominal yields
            %JJ_r = JJ_n + cumsum((-repmat(B_X_pi, T0, 1) + ...
            %    [zeros(T0*K,1) repelem(pdf_z./s_n(1:end-1)',K,1).*repmat(B_X_for(:,2:end),T0,1).*BXcumSig2BXcum_pi(2:end) ] ), 2)' ./ (1:Mm)';

            %JJ_r = cumsum((JJ_f - repmat(B_X_pi, T0, 1)), 2)' ./ (1:Mm)'; % previous, incorrect version
        end

        % Second derivativ eof y wrt X
        if nargout > 4
            BXfor_BXfor = zeros(K^2, Mm); % K^2 x maxmat
            for j = 1:Mm
                BXfor_BXfor(:,j) = reshape(B_X_for(:,j)*B_X_for(:,j)',[],1);
            end
            JJ_nn = cumsum([zeros(K^2*T0,1) repelem(pdf_z./s_n(1:end-1)',K^2,1)] .* repmat(BXfor_BXfor,T0,1), 2)' ./ (1:Mm)'; % T*K x maxmat
            JJ_rr = JJ_nn - ...
                cumsum(([zeros(T0*K^2,1) repelem(pdf_z.*z_n./s_n(1:end-1)'.*BXcumSig2BXcum_pi(2:end),K^2,1).*repmat(BXfor_BXfor(:,2:end),T0,1) ] ), 2)' ./ (1:Mm)';
        end
    end

