function k = kernel(x1,x2,sigmaF,l)
%squared exponential / Gaussian kernel
% Hyperparameter sigma & l
    k = sigmaF^2 * exp(-  (x1-x2)*(1/(2*l^2))*(x1-x2).');
end

