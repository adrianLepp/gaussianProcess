%   %This is what I did ...
function h = h1(x)
   
    A = 0.0154;
    g = 9.81;
   
    %h = 1 / A * x(1,4); %1
    %h = 1 / A * sign(x(1,3)-x(1,2))*sqrt(2*g*abs(x(1,3)-x(1,2))); %2
    h = 1 / A *sign(x(1,1)-x(1,3))*sqrt(2*g*abs(x(1,1)-x(1,3))); %3
end

