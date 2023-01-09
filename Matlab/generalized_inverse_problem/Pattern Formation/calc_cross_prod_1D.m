function [outputArg1,outputArg2] = calc_cross_prod_1D(L1,L2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global N P delta_t Cq
outputArg1 = zeros(N,P);
outputArg2 = zeros(N,P);
for i=1:N
    sous_L1 = L1(i,:);sous_L2 = L2(i,:);
    outputArg1(i,:) = delta_t*Cq*sous_L1(:);
    outputArg2(i,:) = delta_t*Cq*sous_L2(:);
end
end

