% Load the csv
data = readmatrix('ODE_solution.csv');
t = data(:,1);
M_t = data(:,2:end);

% Fraction of fluorescent protein
f = 100e-9/50e-6; %eg 100nM labeled in 50uM. You could just normalize and make this go away

% Number of monomers per species
n = [1 2 3 4 5 6];

% Predicted number of fluorescent monomers per species
P = M_t .* n .* f

% Predicted G0(t)
G0_pred = sum(P.^2,2) ./ (sum(P,2).^2);

% I should play with G0

% Predicted diffusion time
tau_n = [1 2^v 3^v 4^v 5^v 6^v];
tauD_pred = sum(P .* tau_n, 2) ./ sum(P,2);

% Compare predicted to experimental data
figure;
subplot(2,1,1); plot(t, G0_pred,'-'); hold on; plot(t_exp, G0_exp,'o');
xlabel('Time (s)'); ylabel('G0'); legend('Model','Data');

subplot(2,1,2); plot(t, tauD_pred,'-'); hold on; plot(t_exp, tauD_exp,'o');
xlabel('Time (s)'); ylabel('\tau_D (s)'); legend('Model','Data');
