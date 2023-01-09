function [U,C] = calc_Psi(Coefs, CrossProduct1,CrossProduct2,U_prev,C_prev)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global h N diff_u diff_c chi delta_t;
rho = exp(Coefs(1));delta = exp(Coefs(2));alpha = exp(Coefs(3));beta = exp(Coefs(4));gamma = exp(Coefs(5));
A = diff_u*delta_t*(2*eye(N)-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1))/h^2;
B = diff_c*delta_t*(2*eye(N)-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1))/h^2;
C = zeros(N,N);
D = chi*delta_t*(2*eye(N)-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1))/h^2;

Grande_matrice = [A,D;C,B];
f = @(x,y) [rho*x.*(1-x)-delta*x,alpha*x-beta*y-gamma*x.*y];
k1 = f(U_prev,C_prev); k1u = k1(:,1); k1c = k1(:,2);
k2 = f(U_prev+delta_t/2*k1u,C_prev+delta_t/2*k1c); k2u = k2(:,1); k2c = k2(:,2);
k3 = f(U_prev+delta_t/2*k2u,C_prev+delta_t/2*k2c); k3u = k3(:,1); k3c = k3(:,2);
k4 = f(U_prev+delta_t*k3u,C_prev+delta_t*k3c); k4u = k4(:,1); k4c = k4(:,2);
U_inter = U_prev + delta_t/6*(k1u+2*k2u+2*k3u+k4u);
C_inter = C_prev + delta_t/6*(k1c+2*k2c+2*k3c+k4c);
Mat_spar = sparse(Grande_matrice+diag(ones(2*N,1)));
Sec_membre = [U_inter+delta_t*CrossProduct1;C_inter+delta_t*CrossProduct2];
Z = Mat_spar\Sec_membre;
U = Z(1:N);C =Z(N+1:2*N);

end

