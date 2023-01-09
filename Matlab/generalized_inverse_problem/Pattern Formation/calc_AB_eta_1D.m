function [A_eta,B_eta,L1_eta,L2_eta] = calc_AB_eta_1D(Coefs,A0,B0,A_exact,B_exact)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
global N Nm P Big_matrix_B Big_matrix W_eps eps
%Initialisation
number_max_iter = 15;
Big_matrix_B=sparse(Coefs(1)*Big_matrix+eye(N));
[Tout_A,Tout_B] = calc_AB_1D(Coefs,A0,B0,zeros(N,P),zeros(N,P));
Corr = W_eps.*(fun_meas([A_exact(:);B_exact(:)])-fun_meas([Tout_A(:);Tout_B(:)]));  
Corr1 = reshape(Corr(1:Nm*P),[Nm,P]);Corr2 = reshape(Corr(Nm*P+1:2*Nm*P),[Nm,P]);
[Tout_L1,Tout_L2] = calc_adjoints_AB_1D(Coefs,Tout_A,Tout_B,Corr1,Corr2);
[Cross_prod1,Cross_prod2] = calc_cross_prod_1D(Tout_L1,Tout_L2);
[Relev_A,Relev_B] = calc_relev_1D(Tout_L1(:,1),Tout_L2(:,1));
iteration = 1;cond=1;theta=1;
while((~cond || iteration==1) && iteration<=number_max_iter)
    iteration = iteration+1;
    [Tout_A_next,Tout_B_next] = calc_AB_1D(Coefs,A0+Relev_A,B0+Relev_B,Cross_prod1,Cross_prod2);
    Tout_A_next = theta*Tout_A_next + (1-theta)*Tout_A;
    Tout_B_next = theta*Tout_B_next + (1-theta)*Tout_B;
    cond = norm(Tout_A_next(:)-Tout_A(:))<=1e-2*norm(Tout_A(:));
    cond = cond && (norm(Tout_B_next(:)-Tout_B(:))<=1e-2*norm(Tout_B(:)));
    Tout_A=Tout_A_next ; Tout_B = Tout_B_next;
    Corr = W_eps.*(fun_meas([A_exact(:);B_exact(:)])-fun_meas([Tout_A(:);Tout_B(:)]));  
    Corr1 = reshape(Corr(1:Nm*P),[Nm,P]);Corr2 = reshape(Corr(Nm*P+1:2*Nm*P),[Nm,P]);
    [Tout_L1_next,Tout_L2_next] = calc_adjoints_AB_1D(Coefs,Tout_A,Tout_B,Corr1,Corr2);
    Tout_L1_next = theta*Tout_L1_next + (1-theta)*Tout_L1;
    Tout_L2_next = theta*Tout_L2_next + (1-theta)*Tout_L2;
    cond = cond && (norm(Tout_L1_next(:)-Tout_L1(:))<=1e-2*norm(Tout_L1(:)));
    cond = cond && (norm(Tout_L2_next(:)-Tout_L2(:))<=1e-2*norm(Tout_L2(:)));
    Tout_L1=Tout_L1_next ; Tout_L2 =Tout_L2_next;
    [Cross_prod1,Cross_prod2] = calc_cross_prod_1D(Tout_L1,Tout_L2);
    [Relev_A,Relev_B] = calc_relev_1D(Tout_L1(:,1),Tout_L2(:,1));
    if (iteration==number_max_iter)
        theta = 9/10*theta;
        iteration = 0;number_max_iter = round(5/4*number_max_iter);
        [Tout_A,Tout_B] = calc_AB_1D(Coefs,A0,B0,zeros(N,P),zeros(N,P));
        Corr = W_eps.*(fun_meas([A_exact(:);B_exact(:)])-fun_meas([Tout_A(:);Tout_B(:)])); 
        Corr1 = reshape(Corr(1:Nm*P),[Nm,P]);Corr2 = reshape(Corr(Nm*P+1:2*Nm*P),[Nm,P]);
        [Tout_L1,Tout_L2] = calc_adjoints_AB_1D(Coefs,Tout_A,Tout_B,Corr1,Corr2);
        [Cross_prod1,Cross_prod2] = calc_cross_prod_1D(Tout_L1,Tout_L2);
        [Relev_A,Relev_B] = calc_relev_1D(Tout_L1(:,1),Tout_L2(:,1));
        disp('Nombre d iterations maximales atteint');disp(theta);
    end
end

A_eta = Tout_A;B_eta = Tout_B;
L1_eta = Tout_L1; L2_eta = Tout_L2;
end


