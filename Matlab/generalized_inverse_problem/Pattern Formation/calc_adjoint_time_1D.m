function [l1_prev,l2_prev] = calc_adjoint_time_1D(Coefs,Lambda_1_next,Lambda_2_next,A,B,Corr_1,Corr_2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global delta_t N Coef_alp Coef_rho Coef_K Big_matrix_A Big_matrix_B Problem index_1time
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
d=Coefs(1);gamma=Coefs(2);
C1_prev = zeros(N,1); C2_prev = zeros(N,1);
C1_prev(index_1time) = Corr_1;C2_prev(index_1time) = Corr_2;
l1_next = Lambda_1_next; l2_next = Lambda_2_next;
switch Problem
    case 'Simple'
        f_u = @(u,v) -1+2.*u.*v;
        f_v = @(u,v) u.^2;
        g_u = @(u,v) -2.*u.*v;
        g_v = @(u,v) -u.^2;
    case 'Thomas'
        h_u = @(u,v) (Coef_rho*v.*(1+u+Coef_K.*u.^2)-Coef_rho*u.*v.*(1+2*Coef_K.*u))./(1+u+Coef_K.*u.^2).^2;
        h_v = @(u,v) Coef_rho.*u./(1+u+Coef_K.*u.^2);
        f_u = @(u,v) -1-h_u(u,v);
        f_v = @(u,v) -h_v(u,v);
        g_u = @(u,v) -h_u(u,v);
        g_v = @(u,v) -Coef_alp-h_v(u,v);
end
Mat11 = gamma*delta_t*diag(-f_u(A,B));
Mat12 = gamma*delta_t*diag(-f_v(A,B));
Mat21 = gamma*delta_t*diag(-g_u(A,B));
Mat22 = gamma*delta_t*diag(-g_v(A,B));

Giga_matrix = sparse([-Big_matrix_A-Mat11,-Mat12;-Mat21,-Big_matrix_B-Mat22]);
Lambda = Giga_matrix\[-l1_next-delta_t*C1_prev;-l2_next-delta_t*C2_prev];

l1_prev = Lambda(1:N,1);l2_prev = Lambda(N+1:2*N,1);
end

