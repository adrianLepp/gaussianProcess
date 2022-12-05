function [x_post_kk] = particlefilter(xPost,S,y,solveSystem,systemParameter)
%Partikelfilter: 1 step of particlefilter
%   x_post_k    Partikelset k-1
%   S           Partikelanzahl
%   Q_x         Prozessrauschen
%   R_y         Messrauschen
%   z           observation
%   PredMeth    Wie wird der prediction step durchgeührt
%   solveSystem solve the system equations (eg solveThreeTank() ,...)
  
    xPrio = zeros(size(xPost));
    yTheor = zeros(size(xPost)); % change it to dimension of y
    w = zeros(1,length(xPost));
    
    for m = 1 : S
        %% a priori Partikel
        %x_post_k(:,m) = x_post_k(:,m) + sqrt(Q_x) * [randn; randn; randn];
        
        
        [xPrio(:,m),yTheor(:,m)] = solveSystem(xPost(:,m),systemParameter);
        
        %% Gewichte bestimmen
        P_y = 1/((det(2*pi*systemParameter.sigmaY))^(0.5)) * exp(-0.5*(y - y_theor).' * inv(systemParameter.sigmaY) * (y - y_theor));
        w(m) = P_y; 
        if w(m) < 1e-30
            w(m) = 1e-30;
        end
        summe = sum(w);
        w = w./summe;
    end
    %% a posteriori Partikel ziehen
    %summe2 = cumsum(w);
    %for m = 1 : S
    %    x_post_kk(:,m) = x_prio_kk(:,find(rand <= summe2,1));
    %end
    x_post_kk = lowVarianceSampling(xPrio,w);
    
    for m = 1 : S
        x_post_kk(:,m) = x_post_kk(:,m) + sqrt(systemParameter.sigmaX) * [randn; randn; randn];
    end

end

