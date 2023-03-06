classdef UnscentedKalmanFilter2 < FilterClass
    %UNTITLED2 Summary of obj class goes here
    %   Detailed explanation goes here
    
    properties
        %system  DynamicSystemGP
        s
        n
        m
        gamma
        weightsM
        weightsC
    end
    
    methods
        function obj = UnscentedKalmanFilter2(system)
            arguments
                system DynamicSystemGP2
            end
            %UNTITLED2 Construct an instance of obj class
            %   Detailed explanation goes here
            obj.system = system;
            obj.s = (2 * obj.system.n) + 1;
            obj.n = obj.system.n;
            obj.m = obj.system.m;
            
            alpha = 1 ; % my choise in combination with kappa for sigma-points in 68% confidence interval
            kappa = 1000; % my choise in combination with alpha for sigma-points in 68% confidence interval
            beta =0; % optimal for gaussian distribution
            
            lambda = alpha^2 * (obj.n + kappa) - obj.n;
            obj.weightsM = zeros(obj.s,1);
            obj.weightsC = zeros(obj.s,1);
            obj.gamma = sqrt(obj.n + lambda);
            obj.weightsM(1) = lambda / (obj.n + lambda);
            obj.weightsC(1) = obj.weightsM(1) + (1- alpha^2 + beta);
            for i = 2 : obj.s
                obj.weightsM(i) = 1 / (2 * (obj.n + lambda));
                obj.weightsC(i) = obj.weightsM(i);
            end
        end
        
       
        function [mu,Sigma] = prediction(obj,muPost,SigmaPost,yMeas,u,dt)
            
          
            xSigmaPost = obj.sigmaPoints(muPost,SigmaPost,obj.gamma); %   s x n
            
            %x Prio sigma Points
            xSigmaPrio = zeros(obj.s, obj.n); %                           s x n
            [xSigmaPrio(1,:),~,Q] = obj.system.stateTransition(xSigmaPost(1,:),u,dt);
            for i = 2 : obj.s
                [xSigmaPrio(i,:),~,~] = obj.system.stateTransition(xSigmaPost(i,:),u,dt);
            end
            
            %x Prio estimate and variance
            muPrio = zeros(1,obj.n);
            for i = 1 : obj.s
                muPrio = muPrio + obj.weightsM(i) * xSigmaPrio(i,:);
            end
            sigmaPrio = Q;
            for i = 1 : obj.s
                sigmaPrio = sigmaPrio + obj.weightsC(i) * (xSigmaPrio(i,:)-muPrio).' *(xSigmaPrio(i,:)-muPrio);
            end
            xSigmaPrio = obj.sigmaPoints(muPrio,sigmaPrio,obj.gamma);
            
            % y sigma Points
            ySigma = zeros(obj.s, obj.m);
            ySigma(1,:) = obj.system.measurement(xSigmaPrio(1,:),obj.system.sigmaY,obj.system.m);
            R = obj.system.sigmaY;
            for i = 2 : obj.s
                ySigma(i,:) = obj.system.measurement(xSigmaPrio(i,:),obj.system.sigmaY,obj.system.m);
            end
            
            %y Calc with variance
            yCalc = zeros(1,obj.m);
            for i = 1 : obj.s
                yCalc = yCalc + obj.weightsM(i) * ySigma(i,:);
            end
            S = R;
            for i = 1 : obj.s
                S = S + obj.weightsC(i) * (ySigma(i,:)-yCalc).' *(ySigma(i,:)-yCalc);
            end
            
            %Kalman step
            sigmaXY = zeros(obj.n,obj.m); % nx1 * 1xm 
            for i = 1 : obj.s
               sigmaXY = sigmaXY + obj.weightsC(i) * (xSigmaPrio(i,:) - muPrio).'* (ySigma(i,:)- yCalc);
            end
            
            K = sigmaXY * S^-1; %n x m
            mu = muPrio + (K * (yMeas - yCalc).').';
            Sigma = sigmaPrio - K * S * K.';
            
        end
        
        function x = sigmaPoints (obj, mu, Sigma, gamma)
            x = zeros(obj.s,obj.n);
            x(1,:) = mu;
            SigmaSqrt = chol(Sigma);
            for i = 2 : (obj.n + 1)
                
                %x(i,:) = mu + gamma * sqrt(Sigma(i-1,:));
                x(i,:) = mu + gamma * SigmaSqrt(i-1,:);
            end
            for i = (obj.n + 2) : obj.s
                %x(i,:) = mu + gamma * sqrt(Sigma(i-1-obj.n,:));
                x(i,:) = mu + gamma * SigmaSqrt(i-1-obj.n,:);
            end
        end
    end
end



