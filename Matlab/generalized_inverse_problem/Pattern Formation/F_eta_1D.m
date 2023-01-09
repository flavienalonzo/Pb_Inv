function [outputArg1] = F_eta_1D(Coefs,A_0_ex,B_0_ex,A_exact,B_exact)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global Big_matrix Big_matrix_B N
Big_matrix_B=sparse(Coefs(1)*Big_matrix+eye(N));
[A_eta,B_eta,L1_eta,L2_eta] = calc_AB_eta_1D(Coefs,A_0_ex,B_0_ex,A_exact,B_exact);
[outputArg1] = calc_F_1D(A_eta,B_eta,L1_eta,L2_eta);
end
