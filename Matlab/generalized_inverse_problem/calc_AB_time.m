function [outputArg1,outputArg2] = calc_AB_time(Coefs,A_prev,B_prev,Cross1,Cross2)
global delta_t N M Coef_a Coef_b Coef_alp Coef_rho Coef_K Big_matrix_A Big_matrix_B
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
d=Coefs(1);gamma=Coefs(2);
u_prev = A_prev(:); c_prev = B_prev(:);
f = @(u,v) gamma*[Coef_a-u-Coef_rho*u.*v./(1+u+Coef_K*u.^2);Coef_alp*(Coef_b-v)-Coef_rho*u.*v./(1+u+Coef_K*u.^2)];
k1 = f(u_prev,c_prev); k1u = k1(1:N*M,:); k1c = k1(N*M+1:2*N*M,:);
k2 = f(u_prev+delta_t/2*k1u,c_prev+delta_t/2*k1c); k2u = k2(1:N*M,:); k2c = k2(N*M+1:2*N*M,:);
k3 = f(u_prev+delta_t/2*k2u,c_prev+delta_t/2*k2c); k3u = k3(1:N*M,:); k3c = k3(N*M+1:2*N*M,:);
k4 = f(u_prev+delta_t*k3u,c_prev+delta_t*k3c); k4u = k4(1:N*M,:); k4c = k4(N*M+1:2*N*M,:);
A_inter = u_prev + delta_t/6*(k1u+2*k2u+2*k3u+k4u);
B_inter = c_prev + delta_t/6*(k1c+2*k2c+2*k3c+k4c);
A_next = Big_matrix_A\(A_inter+delta_t*Cross1(:));
B_next = Big_matrix_B\(B_inter+delta_t*Cross2(:));
outputArg1 = reshape(A_next,[N,M]);
outputArg2 = reshape(B_next,[N,M]);
end

