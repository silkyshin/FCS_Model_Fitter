function dMdt = model_2(t,M,k)
% ODEs of even order oligomerization model 1-2-4-6
% reduced parameters, assumes k_n is identical
% Inputs: 
%  t - time
%  M - vector of concentrations M = [M1; M2; M4; M6]
%  k - vector of rate constants
%    k = [k1, k_1, k_assoc, K24, K46]
% Outputs:
%  dMdt - derivative of each oligomers

% Unpack concentrations
M1 = M(1);
M2 = M(2);
M4 = M(3);
M6 = M(4);

% Unpack rate constants
k1 = k(1); K1 = k(2); k_1 = k1 / K1;
k_assoc = k(3);
k_2 = k_assoc / k(4);
k_3 = k_assoc / k(5);

% Stepwise ODEs

% Monomer
dM1 = -2*k1*M1^2 + 2*k_1*M2;

% Dimer
dM2 = k1*M1^2 - k_1*M2 - 2*k_assoc*M2^2 - k_assoc*M2*M4 + 2*k_2*M4 + k_3*M6;

% Tetramer
dM4 = k_assoc*M2^2 - k_2*M4 - k_assoc*M2*M4 + k_3*M6;

% Hexamer
dM6 = k_assoc*M2*M4 - k_3*M6;

% Return derivatives
dMdt = [dM1; dM2; dM4; dM6];

end
