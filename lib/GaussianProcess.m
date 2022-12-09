classdef GaussianProcess
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n
        d
        hyperParameter
        kernel
        K
        KyInv
        xD
        yD
        xStd
        xMean
        yStd
        yMean
        
    end
    
    methods
        function obj = GaussianProcess(xD,yD,hyperParameter,kernel)
            %GaussianProcess Construct an instance of this class
            %   Detailed explanation goes here
            obj.n = size(xD,1);
            obj.hyperParameter = hyperParameter;
            obj.kernel = kernel;
            
            obj.xStd = std(xD);
            obj.xMean = mean(xD);
            obj.yStd = std(yD);
            obj.yMean = mean(yD);
            
            obj.xD = (xD - obj.xMean)./ obj.xStd;
            obj.yD = (yD - obj.yMean)./ obj.yStd;
            
            obj.K = CovMatrix(obj.xD, obj.hyperParameter.sigmaF, obj.hyperParameter.l,kernel);
            obj.KyInv = (obj.K + eye(obj.n) * hyperParameter.sigmaN )^-1 ;
           
        end
        
        function [yS,std] = predict(obj,xS)
            %predict Summary of this method goes here
            %   Detailed explanation goes here
            xS = (xS - obj.xMean) ./ obj.xStd ;
            
            ks = zeros(obj.n,1);
            for i = 1 : obj.n
                ks(i) = obj.kernel(obj.xD(i,:),xS,obj.hyperParameter.sigmaF,obj.hyperParameter.l); 
            end

            yS = (ks.' * obj.KyInv * obj.yD);
            std = (obj.kernel(xS,xS,obj.hyperParameter.sigmaF,obj.hyperParameter.l) - ks.'* obj.KyInv * ks);
            
            yS = yS * obj.yStd + obj.yMean;
            std = std * obj.yStd;
        end
    end
end

function K = CovMatrix(x,sigmaF,l,kernel)
%calculate covMatrix for squaredExponentialKernel
%   x           input vector of test inputs
%   sigmaF l    Hyperparameter
    n = length(x);
    K = zeros(n);
    for i = 1 : n
        for j = 1 : n
            K(i,j) = kernel(x(i,:),x(j,:),sigmaF,l);
        end
    end
end

