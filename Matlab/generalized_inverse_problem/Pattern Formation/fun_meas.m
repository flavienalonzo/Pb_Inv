function [Measure] = fun_meas(Psi)
%	Psi est sous forme de vecteur, Measure aussi
%   Detailed explanation goes here
global N Nm P index_1time
trick = ones(Nm,P).*(0:N:N*P-1);
index_tricked = index_1time + trick;
index = index_tricked(:);
Measure = [Psi(index);Psi(N+index)];
end

