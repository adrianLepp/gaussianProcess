function [y_mu,y_s2] = GPpredict_V1(K_y,x,y,xs,sigmaF,l)
%Prediction step for GP
%   K_y = (K + I * sigmaN )^-1
%   x           trainingInputs
%   y           trainingOutputs
%   xs          new Test Input
%   sigmaF l    Hyperparameter
    n = length(x);
    ks = zeros(n,1);
    for i = 1 : n
         ks(i) = squaredExponentialKernel(x(i,:),xs,sigmaF,l); 
    end

%    K_y = (K+sigmaN*eye(n))^-1;
    y_mu = (ks.' * K_y * y);
    y_s2 = (squaredExponentialKernel(xs,xs,sigmaF,l) - ks.'* K_y * ks);
end

