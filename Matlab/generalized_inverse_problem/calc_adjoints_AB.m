function [outputArg1,outputArg2] = calc_adjoints_AB(Coefs,Tout_A,Tout_B,Tout_Corr1,Tout_Corr2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global N M P delta_mea
outputArg1 = zeros(N,M,P);
outputArg2 = zeros(N,M,P);
iter = P-1;
L1_next = zeros(N,M); L2_next = zeros(N,M);
while (iter>=1)
    outputArg1(:,:,iter+1) = L1_next;
    outputArg2(:,:,iter+1) = L2_next;
    A = Tout_A(:,:,iter);B = Tout_B(:,:,iter);
    if (mod(iter+1,delta_mea)==0)
        Corr_1=Tout_Corr1(:,:,iter+1);Corr_2=Tout_Corr2(:,:,iter+1);
    else
        Corr_1 = zeros(N,M); Corr_2 = zeros(N,M);
    end
    [L1_prev,L2_prev] = calc_adjoint_time(Coefs,L1_next,L2_next,A,B,Corr_1,Corr_2);
    L1_next = L1_prev; L2_next = L2_prev;
    iter = iter - 1;
end
outputArg1(:,:,iter+1) = L1_next;
outputArg2(:,:,iter+1) = L2_next;
end

