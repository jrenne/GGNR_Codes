function [x00, delta] = pars_trans_out(x11, x_bounds)

% sigmoid parameter transformation from bounded to original
% x1: vector of N parameters bounded with sigmoid
% x_bounds: Nx2 array of lower and upper bounds
% delta: gradient of the vector of parameters for calculating standard errors

[N1, N2] = size(x11);

if N2>N1
    x11 = x11';
end

x00 = (x_bounds(:,2) - x_bounds(:,1))./(1+exp(-x11)) + x_bounds(:,1); %transforming the sigmoid to x
if nargout>1
    delta = (x_bounds(:,2) - x_bounds(:,1)) ./ ((1+exp(-x11)).^2) .* exp(-x11);
end

if N2>N1
    x00 = x00';
end
