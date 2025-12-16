function dMdt = ode_model2(t,M,k)
% ODEs of dimer-higher order oligomer model 1-2-4-8
% Inputs: 
%  t - time
%  M - vector of concentrations M = [M1; M2; M4; M8]
%  k - vector of rate constants
%    k = [k1, k_1, k2, k_2, k3, k_3]
% Outputs:
%  dMdt - derivative of each oligomers

% Unpack concentrations
M1 = M(1);
M2 = M(2);
M4 = M(3);
M8 = M(4);

% Unpack rate constants
k1 = k(1); k_1 = k(2);
k2 = k(3); k_2 = k(4);
k3 = k(5); k_3 = k(6);

% Stepwise ODEs

% Monomer
dM1 = -2*k1*M1^2 + 2*k_1*M2;

% Dimer
dM2 = k1*M1^2 - k_1*M2 - 2*k2*M2^2 + 2*k_2*M4;

% Tetramer
dM4 = k2*M2^2 - k_2*M4 - 2*k3*M4^2 + 2*k_3*M8;

% Octamer
dM8 = k3*M4^2 - k_3*M8;

% Return derivatives
dMdt = [dM1; dM2; dM4; dM8];

end
