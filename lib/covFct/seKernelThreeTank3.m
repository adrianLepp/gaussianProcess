function k = seKernelThreeTank3(x1,x2,hyperparameter)
%squared exponential / Gaussian kernel
%   sigmaF l    Hyperparameter
%   x1 x2       Inputs (vectors allowed)
    dx13 = x1(1) - x1(3);
    dx32 = x1(3) - x1(2);
    x1 = [dx13,dx32];
    
    dx13 = x2(1) - x2(3);
    dx32 = x2(3) - x2(2);
    x2 = [dx13,dx32];

    k = hyperparameter.sigmaF^2 * exp(-  (x1-x2)*(1/(2*hyperparameter.l^2))*(x1-x2).');
end

