function [e1,e2,e3] = calc_errors(A_eta,B_eta,L1_eta,L2_eta,A_exact,B_exact)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global P N M Ca delta_mea W_eps
Vect = Int_model_error_PF(L1_eta,L2_eta); % donne le premier terme de l'erreur
L1_col=zeros(N*M,P);L2_col=zeros(N*M,P);
for i=1:P
    L1i = L1_eta(:,:,i);L2i = L2_eta(:,:,i);
    L1i_col = L1i(:);L2i_col = L2i(:);
    L1_col(:,i)=L1i_col;
    L2_col(:,i)=L2i_col;
end
Big_lambda = [L1_col;L2_col]; % Lambda sous forme de vecteur
e1 = fun_int_approx(Vect,Big_lambda,2);

Err_cond_ini = [Ca*L1_col(:,1);Ca*L2_col(:,1)]; % Premier terme de l'erreur e2
e2 = fun_int_approx(Err_cond_ini,Big_lambda(:,1),1);%h*Err_cond_ini'*W_a*Err_cond_ini;
    
index_mea = (delta_mea+1:delta_mea:P);
A_exact_mea = A_exact(:,:,index_mea);B_exact_mea = B_exact(:,:,index_mea);
A_mea = A_eta(:,:,index_mea);B_mea = B_eta(:,:,index_mea);
Err_measure = ([A_exact_mea(:);B_exact_mea(:)]-[A_mea(:);B_mea(:)]);
e3 = trace(Err_measure'*W_eps*Err_measure);
end

