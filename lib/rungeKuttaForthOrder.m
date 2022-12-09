function xKK = rungeKuttaForthOrder(xK,systemSolver,parameter,Dt)
%rungeKuttaForthOrder Summary of this function goes here
%   xK
%   systemSolver(xK,parameter):xKK

    steps = 2;
    dt = Dt/steps; % solve in 2 steps

    % Runge-Kutta 4. order
    for tau = dt : dt : Dt
        xtemp = xK;
        dx1 = systemSolver(xtemp, parameter,dt);
        xtemp = xK + dx1 / 2;
        dx2 = systemSolver(xtemp, parameter,dt);
        xtemp = xK + dx2 / 2;
        dx3 = systemSolver(xtemp, parameter,dt);
        xtemp = xK + dx3 / 2;
        dx4 = systemSolver(xtemp, parameter,dt);

        xK = xK + (dx1 + 2 * dx2 + 2 * dx3 + dx4) / 6;
    end
    xKK = xK;
end