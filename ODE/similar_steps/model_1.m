function dMdt = model_1(t,M,k)
% ODEs of stepwise oligomerization 1-2-3-4-5-6
% reduced parameters, assumes k_n is identical
% Inputs: 
%  t - time
%  M - vector of concentrations M = [M1; M2; M3; M4; M5; M6]
%  k - vector of rate constants
%    k = [k1, K1, k_assoc, K23, K34, K45, K56]
%    k1 - monomer dimer forward rate
%    K1 = k1/k_1 (equilibrium for the first step)
%    k_assoc - assume same forward rate for remainder of steps
%    K23, K34, K45, K56 - equilibrium constants for remainder of steps
% Outputs:
%  dMdt - derivative of each oligomers

% Unpack concentrations
M1 = M(1);
M2 = M(2);
M3 = M(3);
M4 = M(4);
M5 = M(5);
M6 = M(6);

% Unpack rate constants
k1 = k(1); K1 = k(2); k_1 = k1 / K1;
k_assoc = k(3);
k_2 = k_assoc / k(4); % K23
k_3 = k_assoc / k(5); % K34
k_4 = k_assoc / k(6); % K45
k_5 = k_assoc / k(7); % K56



% Stepwise ODEs

% Monomer
dM1 = -2*k1*M1^2 + 2*k_1*M2 - k_assoc*M1*M2 + k_2*M3 - k_assoc*M1*M3 + k_3*M4 - k_assoc*M1*M4 + k_4*M5 - k_assoc*M1*M5 + k_5*M6;

% Dimer
dM2 = k1*M1^2 - k_1*M2 - k_assoc*M1*M2 + k_2*M3;

% Trimer
dM3 = k2*M1*M2 - k_2*M3 - k_assoc*M1*M3 + k_3*M4;

% Tetramer
dM4 = k3*M1*M3 - k_3*M4 - k_assoc*M1*M4 + k_4*M5;

% Pentamer
dM5 = k4*M1*M4 - k_4*M5 - k_assoc*M1*M5 + k_5*M6;

% Hexamer
dM6 = k_assoc*M1*M5 - k_5*M6;

% Return derivatives
dMdt = [dM1; dM2; dM3; dM4; dM5; dM6];

end
