function [outputArg1,outputArg2] = calc_AB(Coefs,A_0,B_0,Cross1,Cross2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global N M P delta_t tf Big_matrix_B Big_matrix
outputArg1 = zeros(N,M,P);
outputArg2 = zeros(N,M,P);
Big_matrix_B=sparse(Coefs(1)*Big_matrix+eye(N*M));
iter = 0;
A_prev = A_0;B_prev = B_0;
while ((iter+1)*delta_t<=tf)
    outputArg1(:,:,iter+1) = A_prev;
    outputArg2(:,:,iter+1) = B_prev;
    [A_next,B_next] = calc_AB_time(Coefs,A_prev,B_prev,Cross1(:,:,iter+1),Cross2(:,:,iter+1));
    A_prev = A_next; B_prev = B_next;
    iter = iter + 1;
end
outputArg1(:,:,iter+1) = A_prev;
outputArg2(:,:,iter+1) = B_prev;
end

