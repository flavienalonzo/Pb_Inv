function [L1,L2] = calc_lambda(Coefs, U, C,L1_next, L2_next, Corr1, Corr2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global h N diff_u diff_c chi delta_t;
rho = exp(Coefs(1));delta = exp(Coefs(2));alpha = exp(Coefs(3));beta = exp(Coefs(4));gamma = exp(Coefs(5));
A = diff_u*delta_t*(2*eye(N)-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1))/h^2;
B = diff_c*delta_t*(2*eye(N)-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1))/h^2;
D = chi*delta_t*(2*eye(N)-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1))/h^2;

E = delta_t*(2*rho*diag(U)+delta*diag(ones(N,1)));
F = delta_t*(gamma*diag(C)-alpha*diag(ones(N,1)));
G = zeros(N,N);
H = delta_t*(gamma*diag(U)+beta*diag(ones(N,1)));

f = @(x,y) [E,F;G,H]*[x;y];
k1 = f(L1_next,L2_next); k1u = k1(1:N,1); k1c = k1(N+1:2*N,1);
k2 = f(L1_next+delta_t/2*k1u,L2_next+delta_t/2*k1c); k2u = k2(1:N,1); k2c = k2(N+1:2*N,1);
k3 = f(L1_next+delta_t/2*k2u,L2_next+delta_t/2*k2c); k3u = k3(1:N,1); k3c = k3(N+1:2*N,1);
k4 = f(L1_next+delta_t*k3u,L2_next+delta_t*k3c); k4u = k4(1:N,1); k4c = k4(N+1:2*N,1);
L1_inter = L2_next - delta_t/6*(k1u+2*k2u+2*k3u+k4u);
L2_inter = L2_next - delta_t/6*(k1c+2*k2c+2*k3c+k4c);

Grande_matrice = [A,G;D,B];
Mat_spar = sparse(-Grande_matrice-diag(ones(2*N,1)));
Sec_membre = [-L1_inter+delta_t*Corr1;-L2_inter+delta_t*Corr2];

% Sec_membre = [U_inter+delta_t*CrossProduct1;C_inter+delta_t*CrossProduct2];
Z = Mat_spar\Sec_membre;
L1 = Z(1:N); L2 =Z(N+1:2*N);
end

