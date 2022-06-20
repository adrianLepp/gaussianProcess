function k = squaredExponentialKernel(x1,x2,sigmaF,l)
%squared exponential / Gaussian kernel
%   sigmaF l    Hyperparameter
%   x1 x2       Inputs (vectors allowed)

    k = sigmaF^2 * exp(-  (x1-x2)*(1/(2*l^2))*(x1-x2).');
end