function [outputArg1] = calc_F(A_eta,B_eta,L1_eta,L2_eta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global Big_matrix Coef_a Coef_rho Coef_K Coef_alp Coef_b h delta_t N M P
A = A_eta(:);B = B_eta(:); L1 = L1_eta(:); L2 = L2_eta(:);
f1 = @(u,v) Coef_a-u-Coef_rho.*u.*v./(1+u+Coef_K.*u.^2);
f2 = @(u,v) Coef_alp*(Coef_b-v)-Coef_rho.*u.*v./(1+u+Coef_K.*u.^2);
Vec21 = -f1(A,B);
Vec22 = -f2(A,B);
Vec12 = zeros(N*M*P,1);
for i=1:P
    sous_B = B_eta(:,:,i);
    Vec12((i-1)*N*M+1:i*N*M,1) = 1e10*sparse(Big_matrix)*sous_B(:);
end
outputArg1 = h*delta_t*[sum(Vec12.*L2);sum(1e10*Vec21.*L1)+sum(1e10*Vec22.*L2)];
end

