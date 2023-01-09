function [outputArg1,outputArg2] = calc_relev(L1_t0,L2_t0)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global N M Ca
Corr_A_t0 = Ca*L1_t0(:);
Corr_B_t0 = Ca*L2_t0(:);
outputArg1 = reshape(Corr_A_t0,[N,M]);
outputArg2 = reshape(Corr_B_t0,[N,M]);
end

