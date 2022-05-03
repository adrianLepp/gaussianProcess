function logL = logLikelihood_V1(y,K,sigmaN)
    %Log Likelihood Version 1
    n = length(y);
    Ky = K + sigmaN^2 * eye(n);
    logL = -1/2 * y.' * Ky^-1 * y - 1/2 * log(det(Ky)) - 1/2*n* log(2*pi);
end
