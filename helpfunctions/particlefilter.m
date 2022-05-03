function [x_post_kk] = particlefilter(x_post_k,S,Q_x,R_y,z)
%Partikelfilter: 1 step of particlefilter
%   x_post_k    Partikelset k-1
%   S           Partikelanzahl
%   Q_x         Prozessrauschen
%   R_y         Messrauschen
%   z           observation
%   PredMeth    Wie wird der prediction step durchgeührt
  
    x_prio_kk = zeros(size(x_post_k));
    x_post_kk = zeros(size(x_post_k));
    w = zeros(1,length(x_post_k));
    
    for m = 1 : S
        %% a priori Partikel
        %x_post_k(:,m) = x_post_k(:,m) + sqrt(Q_x) * [randn; randn; randn];
        
        
        x_prio_kk(:,m) = solveDreitank(x_post_k(:,m));
        
        %% Gewichte bestimmen
        z_theor = diag([1 1 1]) * x_prio_kk(:,m); %y_theor
        P_y = 1/((det(2*pi*R_y))^(0.5)) * exp(-0.5*(z - z_theor).' * inv(R_y) * (z - z_theor));
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
    x_post_kk = lowVarianceSampling(x_prio_kk,w);
    
    x_post_kk(:,m) = x_post_kk(:,m) + sqrt(Q_x) * [randn; randn; randn];

end

