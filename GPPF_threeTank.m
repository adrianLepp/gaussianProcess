% Gaussian Process Particle Filter
% with explizit basis functions g() and residual GP f()
%Author Adrian Lepp
% Last change: 24.03.2022

% new Version with new and correct solveThreeTank function


%% clear
close all
clear 
clc

%% load parameters
    load 'threeTankData.mat'
    load 'dreiTank.mat'
    dt = parameter.dt;
    n = length(xTrain);

%% init
S = 10; % change number of particles
T = 100; % Simulationsdauer (Sekunden)

t = T/dt;

% initial state
x0 = zeros(3,1);
dx0 = zeros (3,1);
u0 = uTrain(1);

% noise
sigmaX_default = diag([1e-09 1e-09 1e-09]);    % Systemrauschen
s2 = zeros(3,S);
s2Post = zeros(3,S);
sigmaY = diag([5e-07 5e-07 5e-07]);    % Messrauschen
sigmaY_default = 5e-07;
V_mf = 1e-12;                       % Varianz des Fehlermittelwerts


xPost = x0 .*ones(3,S) + sqrt(parameter.sigmaX)* randn(3,S);
dx = zeros(3,S);
dxPost = zeros(3,S);
yEst = zeros(3,1);
xPrio = zeros(3,S);
wPost = zeros(1,S);

%output values
xOut = zeros(3,t);
yOut = zeros(3,t);
xEst = zeros(3,t);
dxEst = zeros(3,t);
sigmaXout = zeros(3,t);

%% GP init: set hyperparameters and calculate the covariance Matrix for test inputs

    theta = zeros(3,3,2); % number of hyperparameters, dimension of x, number of GP's
    %sigmaF
    theta(1,:,1) = 100;
    theta(1,:,2) = 0.1;
    %l
    theta(2,:,1) = 2;
    theta(2,:,2) = 5;
    %sigmaN
    theta(3,:,1) = 0.35;
    
    K_ux = zeros(n,n,3);
    K_x= K_ux;
    K_y = K_ux;
    K_dx = K_ux;
    
    %GP 1: Prediction
    for i = 1 :3
        K_ux(:,:,i) = CovMatrix([uTrain,xTrain],theta(1,i,1),theta(2,i,1)); %K
        K_dx(:,:,i) = (K_ux(:,:,i) + theta(3,i,1)*eye(n))^-1;               %(K + I *sigmaN)^-1
        logLikelihood_V1(dxTrain(:,i),K_ux(:,:,i),theta(3,i,1))
    end
    
    %GP 2: Observation
    K_x = zeros(n,n,3);
    for i = 1 :3
        K_x(:,:,i) = CovMatrix(xTrain,theta(1,i,2),theta(2,i,2));   %K
        K_y(:,:,i) = (K_x(:,:,i) + theta(3,i,2)*eye(n))^-1;         %(K + I *sigmaN)^-1
        logLikelihood_V1(yTrain(:,i),K_x(:,:,i),theta(3,i,2))
    end   
    
%% PF

for k = 1 : t
    %Simulation reales System
    if k == 1
        [xOut(:,k),y] = solveThreeTank(x0,parameter);
    else
        [xOut(:,k),y] = solveThreeTank(xOut(:,k-1),parameter);
    end
    %% Partikelfilter loesen
    
    for l = 1 : S
        %% a priori Partikel
        
        %GP für Systemgleichung / prediction model
        for i = 1 : 3
          [dx(i,l),s2(i,l)] = GPpredict_V1(K_dx(:,:,i),[uTrain,xTrain],dxTrain(:,i),[u0,xPost(:,l).'],theta(1,i,1),theta(2,i,1));  
        end
        
        % Der Gp bestimmt nur dx, daher Addition mit x_post_k-1
        xPrio(:,l) = xPost(:,l) + dx(:,l) + sqrt(parameter.sigmaX) * [randn; randn; randn];
        
        %% Gewichte bestimmen
        % GP Für Ausgangsgleichung  / observation model
        for i = 1 :3
            [yEst(i),sigmaY(i,i)] = GPpredict_V1(K_y(:,:,i),xTrain,yTrain(:,i),xPrio(:,l).',theta(1,i,2),theta(2,i,2));   
        end

        wPost(l) = 1/((det(2*pi*parameter.sigmaY))^(0.5)) * exp(-0.5*(yOut(:,k) - yEst).' * inv(parameter.sigmaY) * (yOut(:,k) - yEst));
        
        if wPost(l) < 1e-30
            wPost(l) = 1e-30;
        end
        summe = sum(wPost);
        wPost = wPost./summe;
    end
    %% a posteriori Partikel ziehen 
    resampledParticles = lowVarianceSampling([xPrio; dx; s2],wPost);
    xPost = resampledParticles(1:3,:);
    dxPost = resampledParticles(4:6,:);
    s2Post = resampledParticles(7:9,:);
    
    xEst(:,k) = [mean(xPost(1,:)); mean(xPost(2,:)); mean(xPost(3,:))];
    dxEst(:,k) = [mean(dxPost(1,:)); mean(dxPost(2,:)); mean(dxPost(3,:))];
    sigmaXout(:,k) = [mean(s2Post(1,:)); mean(s2Post(2,:)); mean(s2Post(3,:))];
    
end
 time = linspace(dt,T,t);
 
 figure(1)
 plot(time, xOut(1,:),time, xOut(2,:),time, xOut(3,:), time, xEst(1,:), time, xEst(2,:), time, xEst(3,:));
 legend('x_1','x_2','x_3','x_1 est','x_2 est','x_3 est');

%Die Varianz ist sehr hoch hier: Das liegt an der Wahl der Hyperparameter.
%Diese sind momentan schwer zu bestimmen.
figure(2)
f1 = [dxEst(1,:)+2*sqrt(sigmaXout(1,:)), flip(dxEst(1,:)-2*sqrt(sigmaXout(1,:)),2)];
f2 = [dxEst(2,:)+2*sqrt(sigmaXout(2,:)), flip(dxEst(2,:)-2*sqrt(sigmaXout(2,:)),2)];
f3 = [dxEst(3,:)+2*sqrt(sigmaXout(3,:)), flip(dxEst(3,:)-2*sqrt(sigmaXout(3,:)),2)];
fill([time, flip(time,2)], f1, [7 7 7]/8)
hold on; 
fill([time, flip(time,2)], f2, [7 7 7]/8)
fill([time, flip(time,2)], f3, [7 7 7]/8)

plot(time, dxEst(1,:), 'k', time, dxEst(2,:), 'b', time, dxEst(3,:), 'r');
legend('Varianz dx1','Varianz dx2','Varianz dx3', 'dx1 geschätzt', 'dx2 geschätzt', 'dx3 geschätzt');
xlabel('Zeit t /s')
ylabel('Fuellstand /m')
hold off;

 
 