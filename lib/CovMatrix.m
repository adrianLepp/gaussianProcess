function K = CovMatrix(x,hyperparameter,kernel)
%calculate covMatrix for squaredExponentialKernel
%   x           input vector of test inputs
%   sigmaF l    Hyperparameter
    n = length(x);
    K = zeros(n);
    for i = 1 : n
        for j = 1 : n
            K(i,j) = kernel(x(i,:),x(j,:),hyperparameter);
        end
    end
end
