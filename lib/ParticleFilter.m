classdef ParticleFilter < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        s
        system % DynamicSys
        xPost
    end
    
    methods
        function obj = ParticleFilter(system,s)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.system = system;
            obj.s = s;
            obj.xPost = zeros(s,system.n);
            for l = 1 : s
                obj.xPost(l,:) = randn(1,system.n) * system.sigmaX;
            end
            
        end
        
        function xMean = prediction(obj,y,dt)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj ParticleFilter
                y
                dt
            end
            n = obj.system.n;
            xPrio = zeros(obj.s,n);  
            w = zeros(1,obj.s);
    
            for l = 1 : obj.s
                xPrio(l,:) = obj.system.stateTransition(obj.xPost(l,:),dt);
                yCalc = obj.system.measurement(xPrio(l,:));
        
                %w(1,l) = 2*pi ^(- n / 2) * det(obj.system.sigmaY)^(-1/2) * exp(-1/2 * (y - yCalc) * inv(obj.system.sigmaY) * (y - yCalc).');
                w(1,l) = 2*pi ^(- n / 2) * det(obj.system.sigmaY)^(-1/2) * exp(-1/2 * (y - yCalc) * inv(obj.system.sigmaY) *  (y - yCalc).');
                if w(1,l) < 1e-30
                    w(1,l) = 1e-30;
                end

            end
            summe = sum(w);
            w = w./summe;

            obj.xPost = lowVarianceSampling(xPrio,w,obj.s,obj.system.n);

            for l = 1 : obj.s
                obj.xPost(l,:) = obj.xPost(l,:) + randn(1,3) *  sqrt(obj.system.sigmaX);
            end
            
            xMean = zeros(1,n);
            for i = 1 : n
                xMean(1,i) = mean(obj.xPost(:,i));
            end         
        end
    end
end

function xPost = lowVarianceSampling(xPrio,w,s,n)
%lowVarianceSampling: Resampling Method from Probalibistic Robots by Thrun
%for PF
%   xPrio: Prio Particleset
%   w: weights
    M = s; 
    xPost  = zeros(s,n);
    r = rand(1)*M^-1;
    c = w(1);
    i = 1;
    j = 1;

    for m = 1 : M
        U = r + (m-1)*M^-1;
        while U > c
            i = i + 1;
            c = c + w(i);
        end
        xPost (j,:) = xPrio(i,:);
        j = j + 1;
    end
end

