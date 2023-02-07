classdef DynamicSystemGP < DynamicSys
% Does not implement absMethodOne
% defined as abstract in AbsClass
    properties
      transitionGP(:,1) GaussianProcess
      measurementGP(:,1) GaussianProcess
    end
    methods
      function obj = DynamicSystemGP(xD,uD,dxD,yD,transHyperparameter,transKernel,transMean,measHyperparameter,measKernel,measMean,sigmaX,sigmaY)
           obj.n = size(dxD,2);
           obj.m = size(yD,2);
           obj.sigmaX = sigmaX;
           obj.sigmaY = sigmaY;

           for i = 1 : obj.n
               obj.transitionGP(i) = GaussianProcess([xD,uD],dxD(:,i),transHyperparameter{i},transKernel{i},transMean{i});
           end
           
           for i = 1 : obj.m
               obj.measurementGP(i) = GaussianProcess([xD,uD],yD(:,i),measHyperparameter{i},measKernel{i},measMean{i});
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
      
      function [yKK,sigma] = measurement(obj,xKK,u)
          yKK = zeros(1,obj.m);
           sigma = zeros(obj.m,obj.m);
          for i = 1 : obj.m
              [yKK(i),sigma(i,i)] = obj.measurementGP(i).predict([xKK,u]);
          end
      end
   end
end