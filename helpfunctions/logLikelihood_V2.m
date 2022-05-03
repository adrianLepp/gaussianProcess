function logL = logLikelihood_V2(y,K,sigmaN)
    %Log Likelihood Version 2 with Cholesky
    n = length(y);
    L = chol(K+sigmaN ^2 * eye(n));
    alpha = L.'\(L\y);
    LL = zeros(n,1);
    for i = 1 : n
        LL(i)=log(L(i,i)); 
    end
    logL = -1/2*y.'*alpha-sum(LL)-n/2*log(2*pi);
end




