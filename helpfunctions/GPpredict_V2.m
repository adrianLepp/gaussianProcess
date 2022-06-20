function [y_mu,y_s2] = GPpredict_V2(K,x,y,xs,sigmaF,sigmaN,l)
%Prediction step for GP V2: should be faster than V1 but maybe not always
%computable
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
    L = chol(K+sigmaN^2*eye(n));
    alpha = L.'\(L\y);
    
    v = L\ks;

    y_mu = ks.' * alpha;
    y_s2 = squaredExponentialKernel(xs,xs,sigmaF,l) - v.'*v;
end

