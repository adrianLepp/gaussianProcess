function dx1 = meanDx1R(x,hyperparameter)
%meanDx1 Mean Fct for dx1 of threeTankSystem
%   calculates dx = h(x)^T * beta
%   x       : 1x3 input vector
%   h       : 2x1 basis fct
%   beta    : 2x1 parameters [u; c13]
    beta = hyperparameter.beta;

    xR = [0.5322    0.0931    0.3134 5.064000000000000e-05];
    x = x - xR;
    h = zeros(2,1);
    h(2,1) = 307.4504 *(x(1,3)-x(1,1));
    h(1,1) = 64.9351 * x(1,4);
    
    dx1 = h.' * beta;
end

