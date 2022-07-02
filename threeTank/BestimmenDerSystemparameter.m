% Bestimmen der Systemparameter ZRG Pa1
%Adrian Lepp
% 25.11.20

close all
clear all
clc

%% Bestimmen der Sensorkennlinie h(U)= b - m * U 

load ('Leerlaufen_400_0.mat');
%plot(Zeit, Messung);

U_tank1 = Messung(:,1);
U_tank2 = Messung(:,2);
U_tank3 = Messung(:,3);

m_tank1 = 0.4/(U_tank1(1078)-U_tank1(1));
m_tank2 = 0.4/(U_tank2(1078)-U_tank2(1));
m_tank3 = 0.4/(U_tank3(1078)-U_tank3(1));

b_tank1 = m_tank1*U_tank1(1078);
b_tank2 = m_tank2*U_tank2(1078);
b_tank3 = m_tank3*U_tank3(1078);

h1 = b_tank1 - m_tank1 * U_tank1;
h2 = b_tank2 - m_tank2 * U_tank2;
h3 = b_tank3 - m_tank3 * U_tank3;

figure 
plot(U_tank1,h1);
plot(U_tank2,h2);
subplot(h3);
title('Sensorkennlinie')
xlabel('U /V')
ylabel('h /m')
legend('Tank1')

%% Bestimmen der Pumpenleistung Qz1 in m^3/s

A=0.0154; %m^2
g = 9.81;

load('Pumpenleistung.mat');
plot(Zeit, Messung);

t_delta = Zeit(500)-Zeit(200);
h_delta = (b_tank1 - m_tank1 * Messung(500,1))-(b_tank1 - m_tank1 * Messung(200,1));

Qz1 = h_delta*A/t_delta; %m^3/s

%% Bestimmen von c_A2


load('Ventilquerschnitt.mat');% Wasser läuft von Tank 1 in 3 (c_ und von Tank 2 in Reservoir


%plot(Zeit, Messung);

U_1zu3 = Messung(:,1);
U_von2 = Messung(:,2);
U_3von1 = Messung(:,3);

h1 = b_tank1 - m_tank1 * U_1zu3;
h2 = b_tank2 - m_tank2 * U_von2;
h3 = b_tank3 - m_tank3 * U_3von1;


for n = 2:1400
    dh2(n) = (h2(n+1)-h2(n-1))/0.2;
    ca2(n) = -(dh2(n)*A)/(sqrt(2*g*h2(n)));
end
 ca2mean = mean(ca2(100:800));
 
 for n = 2:1400
    dh3(n) = (h3(n+1)-h3(n-1))/0.2;
    c13(n) = (dh3(n)*A)/(sqrt(2*g*(h1(n)-h3(n))));
end
 c13mean = mean(c13(100:800));
 
figure
load('initial_400_Tank_1.mat');

h1 = b_tank1 - m_tank1 * Messung(:,1);
h2 = b_tank2 - m_tank2 * Messung(:,2);
h3 = b_tank3 - m_tank3 * Messung(:,3);


%plot(Zeit,h1)
