function [w_neu,Emin] = scaledConjugateGradient(E,dE,w,x,y,sigma,lambda)
%scaledConjugateGradient: Minimize error function E through optimization of
%weights w
%   w       Weights
%   alpha   step size
%   r       steepest descent - dE
%   p       orthogonal projection of r
%   k       step
%   lambda  scale parameter

% 1. Initialisierung
r = - dE(w,x,y);
p = r;
k = 1;

N = length(w);

lambdaB = 0;

success = true;
abbruch = false;

while abbruch == false
    pL = norm(p);
    w_neu = w;
    p_neu = p;
    r_neu = r;
    % 2. Calculate second order Information?
    if success == true
        %disp('2.')
        sigma_k = sigma /pL;
        s = (dE(w+sigma_k*p,x,y)-dE(w,x,y))/sigma_k;
        delta = p.'*s;
    end
    % 3. Scale delta
    delta = delta + (lambda - lambdaB)*pL^2;
    % 4. Make Hessian matrix positive Definite?
    if delta <= 0
        %disp('4.')
        lambdaB = 2*(lambda - delta/pL^2);
        delta = - delta + lambda* pL^2;
        lambda = lambdaB;
    end
    % 5. Calculate step size
    mu = p.'*r;
    alpha = mu/delta;
    % 6. Calculate comparision parameter
    Delta = 2* delta* (E(w,x,y)-E(w+alpha*p,x,y))/mu^2
    
    % 7. Error reduction?
    if Delta >= 0
        %disp('7.')
        w_neu = w + alpha*p;
        r_neu = -dE(w_neu,x,y);
        lambdaB = 0;
        success = true;
         if mod(k,N) == 0 % restart algorithm
             disp('restart')
             p_neu =r_neu;
         else
            beta = (norm(r_neu)^2-r_neu.'*r)/mu;
            p_neu = r_neu + beta *p;
         end
        if Delta >= 0.75 %reduce scale parameter
            disp('reduce scale')
            lambda =  1/4*lambda;
        end
    else
        lambdaB = lambda;
    end
    % 8. increase scale Parameter?
    if Delta < 0.25
        disp('increase scale')
        lambda = lambda + (delta*(1-Delta)/pL^2);
    end
    % 9. steepest descend direction r != 0?
    if r ~= 0 % Problem: r wird nicht gleich null, daher nie ein Abbruch. Toleranz einfügen
        disp('9.')
        k = k+1;
        w = w_neu
        r = r_neu
        p = p_neu;
        Emin = E(w_neu,x,y)
    else
        abbruch = true;
    end
    if k > N
        abbruch = true;
    end
end


end

