function [sigmaF, sigmaN, W] = minimizeLogLik(x,y,sigmaF,sigmaN,W,eta,steps)
% minimizeLogLik Try to optimize hyperparameters. not working
%   Gradient descent wird angewandt. Evtl einzelne Rechenschritte sind
%   übernehmbar aber Gesamt konzept müsste neu durchdacht werden
    n = length(x);
    delKW = zeros (n);
    delKsigmaF = zeros(n);
    delKsigmaN = zeros (n);
    
    k = 1;
    while true
       k=k+1;
       K = CovMatrix(x,sigmaF,W); 
       LogL = logLikelihood_V1(y,K,sigmaN)
       
       delKsigmaN = 2 * sigmaN * eye(n);
        for i = 1 : n
            for j = 1: n
            d =(x(i,:)-x(j,:))*((x(i,:)-x(j,:)).');
            delKsigmaF(i,j) = 2* sigmaF * exp(-0.5*W*d);
            delKW(i,j) = -0.5*d * exp(-0.5*W*d);
            end
        end
        a = K^(-1)*y
        delSigmaN = 0.5 * trace((a*a.' - K^(-1))*delKsigmaN); % Hier Optimierung, da nur diagonale Werte Ã¼bernommen
        delSigmaF = 0.5 * trace((a*a.' - K^(-1))*delKsigmaF);
        delW = 0.5 * trace((a*a.' - K^(-1))*delKW);
        
        sigmaF_neu = sigmaF + eta * delSigmaF;
        sigmaN_neu = sigmaN + eta * delSigmaN;
        W_neu = W + eta * delW;

        sigmaF = sigmaF_neu;
        sigmaN = sigmaN_neu;
        W = W_neu;
        
        if k == steps
            break;        
        end     
    end %while
end% Function

