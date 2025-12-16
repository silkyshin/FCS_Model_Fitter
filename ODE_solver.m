addpath('path to ODE folder')
% Initial concentrations (uM)
% model_1
% M0 = [50; 0; 0; 0; 0; 0];
%
% model_2
% M0 = [50; 0; 0; 0];
%
% model_3
% M0 = [50; 0; 0; 0];
%
% model_4
% M0 = [50; 0; 0];
%
% Assumes all monomers at t=0

% Arbitrary rate constants
%
% model_1
% k = [0.1, 0.01, 0.05, 0.005, 0.02, 0.002, 0.01, 0.001, 0.005, 0.0005];
%
% model_2
% k = [0.1, 0.01, 0.05, 0.005, 0.02, 0.002];
%
% model_3
% k = [0.1, 0.01, 0.05, 0.005, 0.02, 0.002];
%
% model_4
% k = [0.1, 0.01, 0.05, 0.005];
%

% time vector (convert to seconds!!!)
% it goes from 0 seconds to the end of a given FCS timepoint experiment
tspan = [0 3600];

% Solve ODEs
[t, M_t] = ode45(@(t,M) odemodel(t,M,k), tspan, M0);

% Plot results
% Alter the number of series accordingly
figure; hold on;
plot(t, M_t(:,1),'r','DisplayName','M1');
plot(t, M_t(:,2),'g','DisplayName','M2');
plot(t, M_t(:,3),'b','DisplayName','M3');
plot(t, M_t(:,4),'c','DisplayName','M4');
plot(t, M_t(:,5),'m','DisplayName','M5');
plot(t, M_t(:,6),'k','DisplayName','M6');
xlabel('Time (s)');
ylabel('Concentration (µM)');
legend;
title('Oligomerization');
