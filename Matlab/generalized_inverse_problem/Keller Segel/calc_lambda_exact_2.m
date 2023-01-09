function [L1,L2] = calc_lambda_exact_2(Coefs, U, C,L1_next, L2_next, Corr1, Corr2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global h N diff_u diff_c chi delta_t;
%   Schéma runge-Kutta d'ordre 4
rho = exp(Coefs(1));delta = exp(Coefs(2));alpha = exp(Coefs(3));beta = exp(Coefs(4));gamma = exp(Coefs(5));
E = delta_t*(2*rho*diag(U)+(-rho+delta)*eye(N));
F = delta_t*(gamma*diag(C)-alpha*eye(N));
G = zeros(N,N);
H = delta_t*(gamma*diag(U)+beta*eye(N));
f = @(x,y) fun_help_lambda(E,F,G,H,x,y);
k1 = f(L1_next,L2_next); k1u = k1(1,:); k1c = k1(2,:);
k2 = f(L1_next+delta_t/2*k1u,L2_next+delta_t/2*k1c); k2u = k2(1,:); k2c = k2(2,:);
k3 = f(L1_next+delta_t/2*k2u,L2_next+delta_t/2*k2c); k3u = k3(1,:); k3c = k3(2,:);
k4 = f(L1_next+delta_t*k3u,L2_next+delta_t*k3c); k4u = k4(1,:); k4c = k4(2,:);
L1_inter = L1_next - delta_t/6*(k1u+2*k2u+2*k3u+k4u);
L2_inter = L2_next - delta_t/6*(k1c+2*k2c+2*k3c+k4c);

%   Résolution implicite de l'équation en L1
A = diff_u*delta_t*(2*eye(N)-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1))/h^2;
M_vit = diag(ones(N-1,1),1)-diag(ones(N-1,1),-1);M_vit(1,1) = -1;M_vit(N,N) = 1;
Vit = 0.5*M_vit*C/h;
Trans = diff_u*chi*delta_t*0.5*diag(Vit)*M_vit/h;
Grande_matrice_L1 = sparse(-A+Trans-eye(N));
L1 = Grande_matrice_L1\(-L1_inter-delta_t*Corr1);

%   Résolution implicite de l'équation en L2
B = diff_c*delta_t*(2*eye(N)-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1))/h^2;

Vec_plus = (eye(N)-diag(ones(N-1,1),1))*L1;Vec_plus(N)=0;
Vec_moins = (eye(N)-diag(ones(N-1,1),-1))*L1;Vec_moins(1)=0;
Vec_0 = zeros(N,1);
L1_plus_pos = max([Vec_plus,Vec_0],[],2);
L1_plus_neg = max([-Vec_plus,Vec_0],[],2);
L1_moins_pos = max([Vec_moins,Vec_0],[],2);
L1_moins_neg = max([-Vec_moins,Vec_0],[],2);
D = diff_u*chi*delta_t*(-diag(L1_plus_neg)-diag(L1_moins_neg)+diag(L1_plus_pos(1:N-1),1)+diag(L1_moins_pos(2:N),-1))/h^2;
Conv = D*U;
Grande_matrice_L2 = sparse(-B-eye(N));
L2 = Grande_matrice_L2\(-L2_inter-Conv-delta_t*Corr2);

end

