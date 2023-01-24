classdef DynamicSystemGP < DynamicSys
% Does not implement absMethodOne
% defined as abstract in AbsClass
    properties
      transitionGP(:,1) GaussianProcess
      measurementGP(:,1) GaussianProcess
    end
    methods
      function obj = DynamicSystemGP(xD,dxD,yD,hyperparameter,kernel,mean,sigmaX,sigmaY)
           obj.n = size(xD,2);
           obj.m = size(yD,2);
           obj.sigmaX = sigmaX;
           obj.sigmaY = sigmaY;

           for i = 1 : obj.n
               obj.transitionGP(i) = GaussianProcess(xD,dxD(:,i),hyperparameter{i},kernel{i},mean{i});
           end
           
           for i = 1 : obj.m
               obj.measurementGP(i) = GaussianProcess(xD,yD(:,i),hyperparameter{i},kernel{i},mean{i});
           end
      end
      
      function [xKK,dxKK, sigma] = stateTransition(obj,xK,dt)
          dxKK = zeros(1,obj.n);
          sigma = zeros(1,obj.n);
          for i = 1 : obj.n
              [dxKK(i), sigma(i)] = obj.transitionGP(i).predict(xK);
          end
          
          xKK = xK + dxKK;
          
      end
      
      function yKK = measurement(obj,xKK)
          yKK = zeros(1,obj.m);
          for i = 1 : obj.m
              yKK(i) = obj.measurementGP(i).predict(xKK);
          end
      end
   end
end