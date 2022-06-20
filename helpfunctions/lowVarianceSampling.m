function [X_post] = lowVarianceSampling(x,w)
%lowVarianceSampling: Resampling Method from Probalibistic Robots by Thrun
%for PF
%   x: Prio Particleset
%   w: weights
    M = length(x); 
    X_post = zeros(size(x));
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
        X_post(:,j) = x(:,i);
        j = j + 1;
    end
end

