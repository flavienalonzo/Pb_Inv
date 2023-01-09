function [c,ceq] = contraintes_nonlin(rho,delta,alpha,beta,gamma)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global diff_u delta_t N h diff_c
ceq=[];
c = zeros(4,1);
%Matrice de Psi1
scale = ones(N,1);scale(1)=2;scale(N)=2;
A = diff_u*delta_t*(2*eye(N)-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1))/h^2;
A(1,1)=diff_u*delta_t/h^2;A(N,N)=diff_u*delta_t/h^2;
c_11 = rho*delta_t./scale;
c_12 = delta*delta_t./scale;
A_spar_Psi1 = sparse(A-diag(c_11)+diag(c_12)+diag(1./scale));
c(1) = -rcond(full(A_spar_Psi1))+1e-10;
%Matrice de Psi2
B = diff_c*delta_t*(2*eye(N)-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1))/h^2;
B(1,1)=diff_c*delta_t/h^2;B(N,N)=diff_c*delta_t/h^2;
c_21 = beta*delta_t./scale;
c_22 = gamma*delta_t./scale;
A_spar_Psi2 = sparse(B+diag(c_22)+diag(c_21)+diag(1./scale));
c(2) = -rcond(full(A_spar_Psi2))+1e-10;
%Matrice de Lambda1
c_31 = delta_t*delta./scale;
A_spar_Lambda1 = sparse(-A-diag(1./scale)-diag(c_31));
c(3) = -rcond(full(A_spar_Lambda1))+1e-10;
%Matrice de Lambda2
c_41 = delta_t*(beta)./scale;
A_spar_Lambda2 = sparse(-B-diag(1./scale)-diag(c_41));
c(4) = -rcond(full(A_spar_Lambda2))+1e-10;

end