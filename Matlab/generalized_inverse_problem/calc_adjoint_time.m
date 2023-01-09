function [outputArg1,outputArg2] = calc_adjoint_time(Coefs,Lambda_1_next,Lambda_2_next,A,B,Corr_1,Corr_2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global delta_t N M Coef_alp Coef_rho Coef_K Big_matrix_A Big_matrix_B
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
d=Coefs(1);gamma=Coefs(2);
u_prev = A(:); c_prev = B(:);
C1_prev = Corr_1(:);C2_prev = Corr_2(:);
l1_next = Lambda_1_next(:); l2_next = Lambda_2_next(:);
h_u = @(u,v) (Coef_rho*v.*(1+u+Coef_K.*u.^2)-Coef_rho*u.*v.*(1+2*Coef_K.*u))./(1+u+Coef_K.*u.^2).^2;
h_v = @(u,v) Coef_rho.*u./(1+u+Coef_K.*u.^2);
f_u = @(u,v) -1-h_u(u,v);
f_v = @(u,v) -h_v(u,v);
g_u = @(u,v) -h_u(u,v);
g_v = @(u,v) -Coef_alp-h_v(u,v);
Mat11 = gamma*delta_t*diag(-f_u(u_prev,c_prev));
Mat12 = gamma*delta_t*diag(-f_v(u_prev,c_prev));
Mat21 = gamma*delta_t*diag(-g_u(u_prev,c_prev));
Mat22 = gamma*delta_t*diag(-g_v(u_prev,c_prev));

Giga_matrix = sparse([-Big_matrix_A-Mat11,-Mat12;-Mat21,-Big_matrix_B-Mat22]);
Lambda = Giga_matrix\[-l1_next-delta_t*C1_prev;-l2_next-delta_t*C2_prev];

l1_prev = Lambda(1:N*M,1);l2_prev = Lambda(N*M+1:2*M*N,1);
outputArg1 = reshape(l1_prev,[N,M]);
outputArg2 = reshape(l2_prev,[N,M]);
end

