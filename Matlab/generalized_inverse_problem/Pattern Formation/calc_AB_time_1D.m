function [A_next,B_next] = calc_AB_time_1D(Coefs,A_prev,B_prev,Cross1,Cross2)
global delta_t N Coef_a Coef_b Coef_alp Coef_rho Coef_K Big_matrix_A Big_matrix_B Problem
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
d=Coefs(1);gamma=Coefs(2);
switch Problem
    case 'Simple'
        f = @(u,v) gamma*[Coef_a-u+u.^2.*v,Coef_b-u.^2.*v];
    case 'Thomas'
        f = @(u,v) gamma*[Coef_a-u-Coef_rho*u.*v./(1+u+Coef_K*u.^2),Coef_alp*(Coef_b-v)-Coef_rho*u.*v./(1+u+Coef_K*u.^2)];
end
k1 = f(A_prev,B_prev); k1u = k1(:,1); k1c = k1(:,2);
k2 = f(A_prev+delta_t/2*k1u,B_prev+delta_t/2*k1c); k2u = k2(:,1); k2c = k2(:,2);
k3 = f(A_prev+delta_t/2*k2u,B_prev+delta_t/2*k2c); k3u = k3(:,1); k3c = k3(:,2);
k4 = f(A_prev+delta_t*k3u,B_prev+delta_t*k3c); k4u = k4(:,1); k4c = k4(:,2);
A_inter = A_prev + delta_t/6*(k1u+2*k2u+2*k3u+k4u);
B_inter = B_prev + delta_t/6*(k1c+2*k2c+2*k3c+k4c);
A_next = Big_matrix_A\(A_inter+delta_t*Cross1);
B_next = Big_matrix_B\(B_inter+delta_t*Cross2);
end

