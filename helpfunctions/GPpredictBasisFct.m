function [dxEst] = GPpredictBasisFct(K_y,x,dx,xs,xP,sigmaF,l,betaEst,H,hVector,parameter)
%GPpredictBasisFct Prediction for one time step of GP with basic fct and unknown parameters 
%   hVector is a fct which calculates h for the specific system   
    hs = hVector(xP,parameter);

    n = length(x);
    ks = zeros(n,1);
    for i = 1 : n
         ks(i) = kernel(x(i,:),xs,sigmaF,l); 
    end

    dxEst = zeros(3,1);
    for i = 1 : 3
        dxEst(i) = hs(:,i).' * betaEst(:,i) + ks.' * K_y * (dx(:,i)-H(:,:,i).' * betaEst(:,i));
    end
end

