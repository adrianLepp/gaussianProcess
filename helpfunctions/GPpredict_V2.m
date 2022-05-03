function [y_mu,y_s2] = GPpredict_V2(K,x,y,xs,sigmaF,sigmaN,l)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    n = length(x);
    ks = zeros(n,1);
    for i = 1 : n
         ks(i) = kernel(x(i,:),xs,sigmaF,l); 
    end     
    L = chol(K+sigmaN^2*eye(n));
    alpha = L.'\(L\y);
    
    v = L\ks;

    y_mu = ks.' * alpha;
    y_s2 = kernel(xs,xs,sigmaF,l) - v.'*v;
end

