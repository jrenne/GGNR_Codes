function [x11] = pars_trans_in(x00, x_bounds)

% sigmoid parameter transformation from original to bounded
% x0: original vector of N parameters
% x_bounds: Nx2 array of lower and upper bounds

[N1, N2] = size(x00);

if N2>N1
    x00 = x00';
end

x11 = - log((x_bounds(:,2) - x_bounds(:,1))./ (x00 - x_bounds(:,1)) - 1);

if N2>N1
    x11 = x11';
end

