function dMdt = model_1(t,M,k)
% ODEs of stepwise oligomerization 1-2-3-4-5-6
% reduced parameters
% Inputs: 
%  t - time
%  M - vector of concentrations M = [M1; M2; M3; M4; M5; M6]
%  k - vector of rate constants
%    k = [k1,k_1,k2,k_2,k3,k_3,k4,k_4,k5,k_5]
%    where K_n = k_n / k_-n
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
k2 = k(3); K2 = k(4); k_2 = k2 / K2;
k3 = k(5); K3 = k(6); k_2 = k3 / K3;
k4 = k(7); K4 = k(8); k_2 = k4 /K4;
k5 = k(9); K5 = k(10); k_2 = k5 /K5;

% Stepwise ODEs

% Monomer
dM1 = -2*k1*M1^2 + 2*k_1*M2 - k2*M1*M2 + k_2*M3 - k3*M1*M3 + k_3*M4 - k4*M1*M4 + k_4*M5 - k5*M1*M5 + k_5*M6;

% Dimer
dM2 = k1*M1^2 - k_1*M2 - k2*M1*M2 + k_2*M3;

% Trimer
dM3 = k2*M1*M2 - k_2*M3 - k3*M1*M3 + k_3*M4;

% Tetramer
dM4 = k3*M1*M3 - k_3*M4 - k4*M1*M4 + k_4*M5;

% Pentamer
dM5 = k4*M1*M4 - k_4*M5 - k5*M1*M5 + k_5*M6;

% Hexamer
dM6 = k5*M1*M5 - k_5*M6;

% Return derivatives
dMdt = [dM1; dM2; dM3; dM4; dM5; dM6];

end
