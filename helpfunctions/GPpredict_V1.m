function [y_mu,y_s2] = GPpredict_V1(K_y,x,y,xs,sigmaF,l)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    n = length(x);
    ks = zeros(n,1);
    for i = 1 : n
         ks(i) = kernel(x(i,:),xs,sigmaF,l); 
    end

%    K_y = (K+sigmaN*eye(n))^-1;
    y_mu = (ks.' * K_y * y);
    y_s2 = (kernel(xs,xs,sigmaF,l) - ks.'* K_y * ks);
end

