function dx1 = meanDx1(x,hyperparameter)
%meanDx1 Mean Fct for dx1 of threeTankSystem
%   calculates dx = h(x)^T * beta
%   x       : 1x3 input vector
%   h       : 2x1 basis fct
%   beta    : 2x1 parameters [u; c13]
    A = 0.0154;
    g = 9.81;
    beta = hyperparameter.beta;

    h = zeros(2,1);
    h(1,1) = 1 / A * x(1,4);
    h(2,1) = - 1 / A * sign(x(1,1)-x(1,3)) * sqrt(2*g*abs(x(1,1)-x(1,3)));
    
    dx1 = h.' * beta;
end

