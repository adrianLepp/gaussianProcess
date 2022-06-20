function K = CovMatrix(x,sigmaF,l)
%calculate covMatrix for squaredExponentialKernel
%   x           input vector of test inputs
%   sigmaF l    Hyperparameter
    n = length(x);
    K = zeros(n);
    for i = 1 : n
        for j = 1 : n
            K(i,j) = squaredExponentialKernel(x(i,:),x(j,:),sigmaF,l);
        end
    end
end

