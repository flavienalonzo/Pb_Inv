function [outputArg1,outputArg2] = calc_AB_long(Coefs,A_0,B_0,delta_f)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global N M P delta_t tf
outputArg1 = zeros(N,M,P);
outputArg2 = zeros(N,M,P);
iter = 0;
A_prev = A_0;B_prev = B_0;
while ((iter+1)*delta_t<=delta_f*tf)
    if (mod(iter,delta_f)==0)
    outputArg1(:,:,iter+1) = A_prev;
    outputArg2(:,:,iter+1) = B_prev;
    end
    [A_next,B_next] = calc_AB_time(Coefs,A_prev,B_prev,zeros(N,M),zeros(N,M));
    A_prev = A_next; B_prev = B_next;
    iter = iter + 1;
end
outputArg1(:,:,P) = A_prev;
outputArg2(:,:,P) = B_prev;
end