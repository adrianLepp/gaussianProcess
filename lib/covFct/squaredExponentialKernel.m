function k = squaredExponentialKernel(x1,x2,hyperparameter)
%squared exponential / Gaussian kernel
%   sigmaF l    Hyperparameter
%   x1 x2       Inputs (vectors allowed)

    k = hyperparameter.sigmaF^2 * exp(-  (x1-x2)*(1/(2*hyperparameter.l^2))*(x1-x2).');
end

