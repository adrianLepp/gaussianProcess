function dx2 = meanDx2R(x,hyperparameter)
%meanDx1 Mean Fct for dx1 of threeTankSystem
%   calculates dx = h(x)^T * beta
%   x       : 1x3 input vector
%   h       : 2x1 basis fct
%   beta    : 2x1 parameters: [c32; c2R]
    beta = hyperparameter.beta;

    xR = [0.5322    0.0931    0.3134 5.064000000000000e-05];
    x = x - xR;
    h = zeros(2,1);
    h(1,1) = 306.4020 *(x(1,3)-x(1,2));
    h(2,1) = - 471.3288* x(1,2);
    
    
    dx2 = h.' * beta;
end