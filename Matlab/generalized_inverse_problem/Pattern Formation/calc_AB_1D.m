function [outputArg1,outputArg2] = calc_AB_1D(Coefs,A_0,B_0,Cross1,Cross2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global N P delta_t tf %sigma_q
outputArg1 = zeros(N,P);
outputArg2 = zeros(N,P);
iter = 0;
A_prev = A_0;B_prev = B_0;
while ((iter+1)*delta_t<=tf)
    outputArg1(:,iter+1) = A_prev ;%+sigma_q*sqrt(delta_t)*randn(N,1);
    outputArg2(:,iter+1) = B_prev ;%+sigma_q*sqrt(delta_t)*randn(N,1);
    [A_next,B_next] = calc_AB_time_1D(Coefs,A_prev,B_prev,Cross1(:,iter+2),Cross2(:,iter+2));
    A_prev = A_next; B_prev = B_next;
    iter = iter + 1;
end
outputArg1(:,iter+1) = A_prev;
outputArg2(:,iter+1) = B_prev;
end

