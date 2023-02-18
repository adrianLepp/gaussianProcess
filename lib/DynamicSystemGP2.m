classdef DynamicSystemGP2
% Does not implement absMethodOne
% defined as abstract in AbsClass
    properties
      transitionGP(:,1) GaussianProcess
      n
      sigmaX
      sigmaY
      m
      measurement
    end
    methods
      function obj = DynamicSystemGP2(xD,uD,dxD,hyperparameter,kernel,mean,measurement,sigmaX,sigmaY,m)
           obj.n = size(dxD,2);
           obj.m = m;
           obj.sigmaX = sigmaX;
           obj.sigmaY = sigmaY;
           obj.measurement = measurement;

           for i = 1 : obj.n
               obj.transitionGP(i) = GaussianProcess([xD,uD],dxD(:,i),hyperparameter{i},kernel{i},mean{i});
           end
           
      end
      
      function [xKK,dxKK, sigma] = stateTransition(obj,xK,u,dt)
          dxKK = zeros(1,obj.n);
          sigma = zeros(obj.n,obj.n);
          for i = 1 : obj.n
              [dxKK(i), sigma(i,i)] = obj.transitionGP(i).predict([xK,u]);
          end
          
          xKK = xK + dxKK;
          
      end
   end
end