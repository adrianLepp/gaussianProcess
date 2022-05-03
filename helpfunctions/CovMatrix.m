function K = CovMatrix(x,sigmaF,l)
    n = length(x);
    K = zeros(n);
    for i = 1 : n
        for j = 1 : n
            K(i,j) = kernel(x(i,:),x(j,:),sigmaF,l);
        end
    end
end

