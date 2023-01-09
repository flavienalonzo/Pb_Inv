function [result] = Resolution(rho,delta,alpha,beta,gamma,u_0,c_0,int_result)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global t0 tf delta_t N P h delta_mea Measures W_eps C_a;
Coefs = [rho,delta,alpha,beta,gamma];
%Coefs_exact = log([0.2,0.1,0.1,0.03,0.08]);

L1_N = zeros(N,1);
L2_N = zeros(N,1);

index = 1;
Psi1 = zeros(N,P);Psi2 = zeros(N,P);
u_prev = u_0;
c_prev = c_0;
Psi1(:,1) = u_0;Psi2(:,1) = c_0;
CrossProduct = zeros(N,1);
while (index*delta_t<=tf)
%     u_t = calc_Psi_1(Coefs, CrossProduct,u_prev,c_prev);
%     c_t = calc_Psi_2(Coefs, CrossProduct,c_prev,u_t);
    [u_t,c_t] = calc_Psi_exact(Coefs, CrossProduct,CrossProduct,u_prev,c_prev);
    Psi1(:,index+1)=u_t;
    Psi2(:,index+1)=c_t;
    u_prev = u_t;
    c_prev = c_t;
    index = index + 1;
end

index = P-1;
Lambda1 = zeros(N,P);Lambda2 = zeros(N,P);
L1_next = L1_N; L2_next = L2_N;
Lambda1(:,P) = L1_N;Lambda2(:,P) = L2_N;
%Corr2 = zeros(N,1);
while(index*delta_t>=t0)
    if (mod(index+1,delta_mea)==0)
        Corr = W_eps*(Measures(:,index+1)-M(Psi1(:,index+1),Psi2(:,index+1)));
        Corr1 = Corr(1:N,1);
        Corr2 = Corr(N+1:2*N,1);
    else
        Corr1 = zeros(N,1);
        Corr2 = zeros(N,1);
    end
%     L2_t = calc_lambda_2(Coefs, Psi1(:,index+1), L1_next,L2_next, Corr2);
%     L1_t = calc_lambda_1(Coefs, Psi1(:,index+1), Psi2(:,index+1),L1_next, L2_t, Corr1);
    [L1_t,L2_t] = calc_lambda_exact_2(Coefs, Psi1(:,index+1), Psi2(:,index+1),L1_next, L2_next, Corr1, Corr2);
    Lambda1(:,index+1) = L1_t;
    Lambda2(:,index+1) = L2_t;
    L1_next = L1_t; L2_next = L2_t;
    index = index - 1;
end

iter = 1; itermax = 1000;
cond = false;
Psi1_prev = Psi1;Psi2_prev = Psi2;Lambda1_prev = Lambda1;Lambda2_prev = Lambda2;
%bar = waitbar(0,'Lets get started');
while (iter<=itermax && ~cond)
    %waitbar(iter/itermax,bar,'Progression');
    %Psi
    index = 1;
    Psi1 = zeros(N,P);Psi2 = zeros(N,P);
    Err_cond_ini = C_a*[Lambda1(:,1);Lambda2(:,1)];
    u_prev = u_0 + (Err_cond_ini(1:N,1)); % à compléter
    c_prev = c_0 + (Err_cond_ini(N+1:2*N,1)); % à compléter
    Psi1(:,1) = u_prev;Psi2(:,1) = c_prev;
    CrossProduct = Int_model_error(Lambda1,Lambda2);
    while (index*delta_t<=tf)
%         u_t = calc_Psi_1(Coefs, CrossProduct(1:N,index),u_prev,c_prev);
%         c_t = calc_Psi_2(Coefs, CrossProduct(N+1:2*N,index),c_prev,u_t);
        [u_t,c_t] = calc_Psi_exact(Coefs, CrossProduct(1:N,index+1),CrossProduct(N+1:2*N,index+1),u_prev,c_prev);
        Psi1(:,index+1)=u_t;
        Psi2(:,index+1)=c_t;
        u_prev = u_t;
        c_prev = c_t;
        index = index + 1;
    end
    %Lambda
    index = P-1;
    Lambda1 = zeros(N,P);Lambda2 = zeros(N,P);
    L1_next = L1_N; L2_next = L2_N;
    Lambda1(:,P) = L1_N;Lambda2(:,P) = L2_N;
    while(index*delta_t>=t0)
        if (mod(index+1,delta_mea)==0)
            Corr = W_eps*(Measures(:,index+1)-M(Psi1(:,index+1),Psi2(:,index+1)));
            Corr1 = Corr(1:N,1);
            Corr2 = Corr(N+1:2*N,1);
        else
            Corr1 = zeros(N,1);
            Corr2 = zeros(N,1);
        end
%         L2_t = calc_lambda_2(Coefs, Psi1(:,index+1), L1_next, L2_next, Corr2);
%         L1_t = calc_lambda_1(Coefs, Psi1(:,index+1), Psi2(:,index+1),L1_next, L2_t, Corr1);
        [L1_t,L2_t] = calc_lambda_exact_2(Coefs, Psi1(:,index+1), Psi2(:,index+1),L1_next, L2_next, Corr1, Corr2);
        Lambda1(:,index+1) = L1_t;
        Lambda2(:,index+1) = L2_t;
        L1_next = L1_t; L2_next = L2_t;
        index = index - 1;
    end
    iter = iter+1;
    
    Vect = Int_model_error(Lambda1,Lambda2); % donne le premier terme de l'erreur
    Big_lambda = [Lambda1;Lambda2]; % Lambda sous forme de vecteur
    e1 = fun_int_approx(Vect,Big_lambda,2);

    Err_cond_ini = C_a*[Lambda1(:,1);Lambda2(:,1)]; % Premier terme de l'erreur e2
    e2 = fun_int_approx(Err_cond_ini,Big_lambda(:,1),1);%h*Err_cond_ini'*W_a*Err_cond_ini;
    
    index_mea = (delta_mea+1:delta_mea:P);
    Err_measure = Measures(:,index_mea) - M(Psi1(:,index_mea),Psi2(:,index_mea));
    e3 = trace(Err_measure'*W_eps*Err_measure);
    
   
    norm_prev = sqrt(norm(Psi1_prev,'fro') + norm(Psi2_prev,'fro') + norm(Lambda1_prev,'fro') + norm(Lambda2_prev,'fro'));
    norm_dif = sqrt(norm(Psi1-Psi1_prev,'fro') + norm(Psi2-Psi2_prev,'fro') + norm(Lambda1-Lambda1_prev,'fro') + norm(Lambda2-Lambda2_prev,'fro'));
    cond = ((norm_dif/norm_prev)<0.001);
    Psi1_prev = Psi1; Psi2_prev = Psi2; Lambda1_prev = Lambda1 ; Lambda2_prev = Lambda2;
    %
%     if(cond)
%         disp('cond obtained')
%         disp(iter);
%         break
%     end
end
total_error = sqrt(e1+e2+e3);
%disp(e1);disp(e2);disp(e3);
Int_approx = calc_int_param(Lambda1, Lambda2, Psi1, Psi2, Coefs);
if (int_result==1)
    result = total_error;
elseif (int_result==2) 
    result = Int_approx;
elseif (int_result==3)
        result = [total_error;Int_approx];
elseif (int_result==4)
    result = [Psi1(:,P),Psi2(:,P)];
    [X,Y]=meshgrid(linspace(0,1,N)',(0:P-1)'*delta_t);
    figure;p1=subplot(2,3,1);mesh(X,Y,Psi1');title('\Psi_1','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);ylabel('time','Rotation',-40,'FontSize', 20);p1.FontSize = 20;colorbar;
    p2=subplot(2,3,2);mesh(X,Y,Measures(1:N,:)');title('\Psi_1 exact','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);ylabel('time','Rotation',-40,'FontSize', 20);p2.FontSize = 20;colorbar;
    p3=subplot(2,3,3);mesh(X,Y,Lambda1');title('\lambda_1','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);ylabel('time','Rotation',-40,'FontSize', 20);p3.FontSize = 20;colorbar;
    p4=subplot(2,3,4);mesh(X,Y,Psi2');title('\Psi_2','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);ylabel('time','Rotation',-40,'FontSize', 20);p4.FontSize = 20;colorbar;
    p5=subplot(2,3,5);mesh(X,Y,Measures(N+1:2*N,:)');title('\Psi_2 exact','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);ylabel('time','Rotation',-40,'FontSize', 20);p5.FontSize = 20;colorbar;
    p6=subplot(2,3,6);mesh(X,Y,Lambda2');title('\lambda_2','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);ylabel('time','Rotation',-40,'FontSize', 20);p6.FontSize = 20;colorbar;
    disp(e1);disp(e2);disp(e3);
elseif (int_result==5)
    result = CrossProduct;
end
if(iter>=itermax)
    disp('Warning: the maximum number of iteration has been reached. The relative error at that step is ');disp((norm_dif/norm_prev));
end
%figure;mesh(Psi1);colorbar;
%figure;mesh(Psi2);colorbar;
%figure;mesh(Lambda1);colorbar;
%figure;mesh(Lambda2);colorbar;
end