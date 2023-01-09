function [outputArg1,outputArg2] = calc_cross_prod(L1,L2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global N M P delta_t
Cq = (1e-10)^2*delta_t*min((1:P),(1:P)');
outputArg1 = zeros(N,M,P);
outputArg2 = zeros(N,M,P);
for i=1:N
    for j=1:M
        sous_L1 = L1(i,j,:);sous_L2 = L2(i,j,:);
        outputArg1(i,j,:) = delta_t*Cq*sous_L1(:);
        outputArg2(i,j,:) = delta_t*Cq*sous_L2(:);
    end
end
end

