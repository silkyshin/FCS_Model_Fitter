function dMdt = model_4(t, M, k)
% ODE system for trimer based nucleation mode: 1 <-> 3 <-> 6
% reduced parameters
%
% Inputs:
%   t - time (required by ode45)
%   M - vector of concentrations: M = [M1; M3; M6]
%   k - vector of rate constants:
%       k = [k1, k_1, k2, k_2]
%
% Outputs:
%   dMdt - derivative of each species

% Unpack concentrations
M1 = M(1);
M3 = M(3);
M6 = M(6);

% Unpack rate constants
k1  = k(1); K1 = k(2); k_1 = k1 / K1;
k2  = k(3); K2 = k(4); k_2 = k2 / K2;

% Stepwise ODEs

% Monomer
dM1 = -3*k1*M1^3 + 3*k_1*M3;

% Trimer
dM3 = 3*k1*M1^3 - k_1*M3 - 2*k2*M3^2 + 2*k_2*M6;

% Hexamer
dM6 = 2*k2*M3^2 - k_2*M6;

% Return derivatives
dMdt = [dM1; dM3; dM6];

end
