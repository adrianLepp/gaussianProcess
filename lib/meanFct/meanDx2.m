function dx2 = meanDx2(x,hyperparameter)
%meanDx1 Mean Fct for dx1 of threeTankSystem
%   calculates dx = h(x)^T * beta
%   x       : 1x3 input vector
%   h       : 2x1 basis fct
%   beta    : 2x1 parameters: [c32; c2R]
    A = 0.0154;
    g = 9.81;
    beta = hyperparameter.beta;

    h = zeros(2,1);
    h(1,1) = 1 / A * sign(x(1,3)-x(1,2))*sqrt(2*g*abs(x(1,3)-x(1,2)));
    h(2,1) = - 1 / A * sqrt(2*g*abs(x(1,2)));
    
    dx2 = h.' * beta;
end