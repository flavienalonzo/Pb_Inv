function [Corr_A_t0,Corr_B_t0] = calc_relev_1D(L1_t0,L2_t0)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global Ca
Corr_A_t0 = Ca*L1_t0;
Corr_B_t0 = Ca*L2_t0;
end

