classdef WaterTank < DynamicSys
   properties
      parameter  
   end
   methods
      function xKK = stateTransition(obj,xK,dt)
          xKK = rungeKuttaForthOrder(xK,@solveThreeTankLinear,obj.parameter,dt);
          for i = 1 : obj.n
              if xKK(1,i) < 0 
                  xKK(1,i) = 0;
              end
          end
          xKK = xKK + randn(1,3) * sqrt(obj.sigmaX);
      end
      
      function yKK = measurement(obj,xKK)
          yKK = xKK + randn(1,3) * sqrt(obj.sigmaY);
      end
      
      function obj = WaterTank(n,m,sigmaX,sigmaY,parameter)
           obj.n = n;
           obj.m = m;
           obj.sigmaX = sigmaX;
           obj.sigmaY = sigmaY;
           obj.parameter = parameter;
      end
   end
end

function dx =  solveThreeTankLinear(x, parameter,dt)
        u1 = parameter.u1;
        u2 = parameter.u2;
        c13 = parameter.c13;
        c32 = parameter.c32;
        c2R = parameter.c2R;
        A = parameter.A;
        g = parameter.g;
        dx = zeros(1,3);

        dx(1) = dt*1/A*(u1-c13*sign(x(1)-x(3))*sqrt(2*g*abs(x(1)-x(3))));
        dx(2) = dt*1/A*(u2 + c32*sign(x(3)-x(2))*sqrt(2*g*abs(x(3)-x(2)))-c2R*sqrt(2*g*abs(x(2))));
        dx(3) = dt*1/A*(c13*sign(x(1)-x(3))*sqrt(2*g*abs(x(1)-x(3)))-c32*sign(x(3)-x(2))*sqrt(2*g*abs(x(3)-x(2))));
end