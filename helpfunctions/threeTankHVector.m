function h = threeTankHVector(x,parameter)
%threeTankHVector
%   calculate the h-Vector of the linear model of the three Tank
%   dx = h(x) * b, where b are the system parameters
%   parameters need to be type of three tank parameters
    
    if (isstruct(parameter) && strcmp(parameter.system,'threeTank'))
        h = zeros(4,3);
        dt = parameter.dt;
        A = parameter.A;
        u = parameter.u;
        g = parameter.g;
    
        h11 = u * 1/A;
        h21 = - 1/A * sign(x(1)-x(3))*sqrt(2*g*abs(x(1)-x(3)));
        h32 = 1/A * sign(x(3)-x(2))*sqrt(2*g*abs(x(3)-x(2)));
        h42 = - 1/A * sqrt(2*g*abs(x(2)));
        h23 = 1/A * sign(x(1)-x(3))*sqrt(2*g*abs(x(1)-x(3)));
        h33 = - 1/A * sign(x(3)-x(2))*sqrt(2*g*abs(x(3)-x(2)));
        
        h(:,1) = dt.*[h11; h21; 0; 0];
        h(:,2) = dt.*[0; 0; h32; h42];  
        h(:,3) = dt.*[0; h23; h33; 0];
    else
        error('parameters are not correct.')
    end
    
end