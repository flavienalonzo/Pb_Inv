function [Vect] = Int_model_error_PF(Lambda1,Lambda2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global N M P delta_t
Cq = (1e-7)^2*delta_t*min((1:P),(1:P)');
L1_col=zeros(N*M,P);L2_col=zeros(N*M,P);Vect = zeros(2*N*M,P);
for i=1:P
    L1i = Lambda1(:,:,i);L2i = Lambda2(:,:,i);
    L1i_col = L1i(:);L2i_col = L2i(:);
    L1_col(:,i)=L1i_col;
    L2_col(:,i)=L2i_col;
end
for k=1:N*M
    Vect(k,i) = Vect(k,i)+Cq(i,:)*L1_col(k,:)';
    Vect(N*M+k,i) = Vect(N*M+k,i)+ Cq(i,:)*L2_col(k,:)';
end

end

