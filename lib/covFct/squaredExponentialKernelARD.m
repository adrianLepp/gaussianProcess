function k = squaredExponentialKernelARD(x1,x2,hyperparameter)
%squared exponential / Gaussian kernel
%   sigmaF l    Hyperparameter
%   x1 x2       Inputs (vectors allowed)
    
% 
%      x1(1,4) = x1(1,4) * 2;
%      x2(1,4) = x2(1,4) * 2;
%      
%      x1(1,1) = x1(1,1) * 2;
%      x2(1,1) = x2(1,1) * 2;
    
    
    
    k = hyperparameter.sigmaF^2 * exp(-1/2 *  (x1-x2)*((hyperparameter.w^2)^-1)*(x1-x2).');
    
   % k = k * abs(1 - (x1(1,4)-x2(1,4))/2);
end

