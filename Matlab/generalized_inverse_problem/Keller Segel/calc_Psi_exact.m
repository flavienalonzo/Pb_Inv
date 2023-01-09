function [U,C] = calc_Psi_exact(Coefs, CrossProduct1,CrossProduct2,U_prev,C_prev)
%	du/dt - d/dx.(D1 d/dx u) + d/dx.(D1 chi u d/dx c) = g1(u,c,theta)
%	+(KqqL)1
%   dc/dt - d/dx.(D2 d/dx c) = g2(u,c,theta) + (KqqL)2

global h N diff_u diff_c chi delta_t;
%   Schéma Runge-Kutta d'ordre 4
rho = exp(Coefs(1));delta = exp(Coefs(2));alpha = exp(Coefs(3));beta = exp(Coefs(4));gamma = exp(Coefs(5));
f = @(x,y) [rho*x.*(1-x)-delta*x,alpha*x-beta*y-gamma*x.*y];
k1 = f(U_prev,C_prev); k1u = k1(:,1); k1c = k1(:,2);
k2 = f(U_prev+delta_t/2*k1u,C_prev+delta_t/2*k1c); k2u = k2(:,1); k2c = k2(:,2);
k3 = f(U_prev+delta_t/2*k2u,C_prev+delta_t/2*k2c); k3u = k3(:,1); k3c = k3(:,2);
k4 = f(U_prev+delta_t*k3u,C_prev+delta_t*k3c); k4u = k4(:,1); k4c = k4(:,2);

U_inter = U_prev + delta_t/6*(k1u+2*k2u+2*k3u+k4u);
C_inter = C_prev + delta_t/6*(k1c+2*k2c+2*k3c+k4c);

%   Résolution implicite de l'équation en c
B = diff_c*delta_t*(2*eye(N)-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1))/h^2;
Grande_mat_c = sparse(B+diag(ones(N,1)));
C = Grande_mat_c\(C_inter+delta_t*CrossProduct2);
%   Résolution implicite de l'équation en u
A = diff_u*delta_t*(2*eye(N)-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1))/h^2;

Vec_plus = (eye(N)-diag(ones(N-1,1),1))*C;Vec_plus(N)=0;
Vec_moins = (eye(N)-diag(ones(N-1,1),-1))*C;Vec_moins(1)=0;
Vec_0 = zeros(N,1);
C_plus_pos = max([Vec_plus,Vec_0],[],2);
C_plus_neg = max([-Vec_plus,Vec_0],[],2);
C_moins_pos = max([Vec_moins,Vec_0],[],2);
C_moins_neg = max([-Vec_moins,Vec_0],[],2);
D = diff_u*chi*delta_t*(-diag(C_plus_neg)-diag(C_moins_neg)+diag(C_plus_pos(1:N-1),1)+diag(C_moins_pos(2:N),-1))/h^2;

Grande_matrice_u = sparse(A-D+diag(ones(N,1)));
U = Grande_matrice_u\(U_inter+delta_t*CrossProduct1);

end


