clear all;close all;
global t0 tf delta_t N Nm P x_space h Coef_a Coef_b Coef_alp Coef_rho Coef_K Coefs_exact index_1time delta_mea Big_matrix Big_matrix_A Big_matrix_B W_eps Ca Problem Cq;
tic;
t0=0e0; tf=2.0e0; delta_t=1e-3; N=100; delta_mea = 100; h=1/(N-1); P = round((tf-t0)/delta_t+1); Nm = N;
x_space = linspace(0,1,N)';
%Problem = 'Simple';
Problem = 'Thomas';
switch Problem
    case 'Simple'
        Coef_a=0.2;Coef_b=2.0;
        Coefs_exact = [50;30];%(d,gamma)
        lb=[40;20];ub = [60;40];A=[];b=[];
        W_eps = 1e2;
        Ca = (1e-2)^2*eye(N);
        Cq = (1e-4)^2*delta_t*min((1:P),(1:P)');
        f = @(u,v) Coef_a-u+u^2*v;g = @(u,v) Coef_b-u^2*v;
        f_u = @(u,v) -1+2*u*v;  f_v = @(u,v) u^2;
        g_u = @(u,v) -2*u*v;    g_v = @(u,v) -u^2;
        x_ex = [Coef_a+Coef_b,Coef_b/(Coef_a+Coef_b)^2];
    case 'Thomas'
        Coef_a=92; Coef_b=64; Coef_alp=1.5; Coef_rho=18.5; Coef_K=0.1;      
        Coefs_exact = [20;30];%(d,gamma)
        lb=[15;25];ub = [25;35];A=[];b=[];              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        W_eps = (1e1)^2; 
        %eps = 1e0*randn(2*N*P,1);
        Ca = (1e-4)^2*eye(N);                       
        Cq = (1e-4)^2*delta_t*min((1:P),(1:P)');    
        index_1time = sort(randperm(N,Nm)');
        f = @(u,v) Coef_a-u-Coef_rho.*u.*v./(1+u+Coef_K.*u.^2);
        g = @(u,v) Coef_alp.*(Coef_b-v)-Coef_rho.*u.*v./(1+u+Coef_K.*u.^2);
        h_u = @(u,v) (Coef_rho*v*(1+u+Coef_K*u^2)-Coef_rho*u*v*(1+2*Coef_K*u))/(1+u+Coef_K*u^2)^2;
        h_v = @(u,v) Coef_rho*u/(1+u+Coef_K*u^2);
        f_u = @(u,v) -1-h_u(u,v);
        f_v = @(u,v) -h_v(u,v);
        g_u = @(u,v) -h_u(u,v);
        g_v = @(u,v) -Coef_alp-h_v(u,v);
        fun = @(x) abs(f(x(1),x(2)))+abs(g(x(1),x(2)));
        options = optimset('TolFun',1e-13,'TolX',1e-13,'MaxFunEvals',1e6,'MaxIter',1e6);
        x_ex = fminsearch(fun,[10,9],options);
end

instab_1 = f_u(x_ex(1),x_ex(2))+g_v(x_ex(1),x_ex(2))<0;
instab_2 = f_u(x_ex(1),x_ex(2))*g_v(x_ex(1),x_ex(2))-f_v(x_ex(1),x_ex(2))*g_u(x_ex(1),x_ex(2))>0;
instab_3 = Coefs_exact(1)*f_u(x_ex(1),x_ex(2))+g_v(x_ex(1),x_ex(2))>0;
instab_4 = (Coefs_exact(1)*f_u(x_ex(1),x_ex(2))+g_v(x_ex(1),x_ex(2)))^2-4*Coefs_exact(1)*(f_u(x_ex(1),x_ex(2))*g_v(x_ex(1),x_ex(2))-f_v(x_ex(1),x_ex(2))*g_u(x_ex(1),x_ex(2)))>0;
if (~(instab_1 && instab_2 && instab_3 && instab_4))
    disp('Attention le problème n_est pas instable !!!!')
end 

%% Quels sont les modes instables du problème ?
P_dc = [f_u(x_ex(1),x_ex(2))^2 2*(2*f_v(x_ex(1),x_ex(2))*g_u(x_ex(1),x_ex(2))-f_u(x_ex(1),x_ex(2))*g_v(x_ex(1),x_ex(2))) g_v(x_ex(1),x_ex(2))^2];
roots_dc = roots(P_dc);
d_crit = max(roots_dc);
d = Coefs_exact(1); gamma=Coefs_exact(2);
qnte = d*f_u(x_ex(1),x_ex(2))+g_v(x_ex(1),x_ex(2));
Mat_fg = [f_u(x_ex(1),x_ex(2)),f_v(x_ex(1),x_ex(2));g_u(x_ex(1),x_ex(2)),g_v(x_ex(1),x_ex(2))];
k1 = sqrt(gamma/(2*d)*(qnte-sqrt(qnte^2-4*d*det(Mat_fg))));
k2 = sqrt(gamma/(2*d)*(qnte+sqrt(qnte^2-4*d*det(Mat_fg))));
fun_modes = @(x) d.*x.^2-gamma*(d*f_u(x_ex(1),x_ex(2))+g_v(x_ex(1),x_ex(2))).*x+gamma^2*det(Mat_fg);
figure;subplot(2,1,1);hold on;fplot(fun_modes,[0,3/2*k2^2]);plot(k1^2,0,'r*');plot(k2^2,0,'r*');xlabel('k^2');ylabel('h(k^2)');
fun_wavenumber_1 = @(x) 0.5*real(-x*(1+d)+gamma*(f_u(x_ex(1),x_ex(2))+g_v(x_ex(1),x_ex(2)))-sqrt((x*(1+d)-gamma*(f_u(x_ex(1),x_ex(2))+g_v(x_ex(1),x_ex(2)))).^2-4*fun_modes(x)));
fun_wavenumber_2 = @(x) 0.5*real(-x*(1+d)+gamma*(f_u(x_ex(1),x_ex(2))+g_v(x_ex(1),x_ex(2)))+sqrt((x*(1+d)-gamma*(f_u(x_ex(1),x_ex(2))+g_v(x_ex(1),x_ex(2)))).^2-4*fun_modes(x)));
fun_wavenumber = @(x) max(fun_wavenumber_1(x),fun_wavenumber_2(x));
n1 = k1/pi;n2 = k2/pi;
list_int = (ceil(n1):floor(n2));
list_wave = (list_int*pi).^2;
subplot(2,1,2);hold on;fplot(fun_wavenumber,[0,3/2*k2^2]);plot(k1^2,0,'r*');plot(k2^2,0,'r*');plot(list_wave,fun_wavenumber(list_wave),'b*');xlabel('k^2');ylabel('Re \lambda');
[l_mode,n_mode] = max(fun_wavenumber(list_wave));
A_0 = x_ex(1)*ones(N,1)+0.01*cos(n_mode*pi*x_space);
B_0 = x_ex(2)*ones(N,1);%+0.01*randn(N,1);
A_0_ex = A_0;%+0.001*randn(N,1);%0.001*randn(1)*cos(n_mode*pi*x_space);
B_0_ex = B_0;%+0.001*randn(N,1);%0.001*randn(N,1);%cos(n_mode*pi*x_space);%+0.15*cos(pi*x_space);
if (lb(1)<d_crit)
    disp('Il y a des valeurs de d qui ne produisent pas d_instabilité');
end
% [Xgrid,Ygrid] = meshgrid(A_0_ex,B_0_ex);
% sdfghj=f(Xgrid,Ygrid);
% figure;meshc(Xgrid,Ygrid,sdfghj);
%% Pré-calcul des matrices pour ne pas le faire dans la boucle en temps
Big_matrix = delta_t*(2*eye(N)-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1))/h^2;
Big_matrix(1,1) = delta_t/h^2;Big_matrix(N,N) = delta_t/h^2;
Big_matrix_A=sparse(Big_matrix+eye(N));
Big_matrix_B=sparse(Coefs_exact(1)*Big_matrix+eye(N));
Spectre_A = eig(Big_matrix_A);Spectre_B = eig(Big_matrix_B);
%% Calcul de la solution avec les paramètres de référence
tic;
%sigma_q = 1e-4;
[A_exact,B_exact] = calc_AB_1D(Coefs_exact,A_0_ex,B_0_ex,zeros(N,P),zeros(N,P));
[X,Y]=meshgrid(x_space,(0:100:P-1)'*delta_t);
figure;
subplot(1,2,1);mesh(X,Y,A_exact(1:N,1:100:P)');xlabel('space');ylabel('time');title('A(t,x)');
subplot(1,2,2);mesh(X,Y,B_exact(1:N,1:100:P)');xlabel('space');ylabel('time');title('B(t,x)');
toc;
%sigma_q = 0e0;
%% Calcul en un eta
Coefs_0 = [23.23;26.75];
Coefs = Coefs_0;tic;
[A_eta,B_eta,L1_eta,L2_eta] = calc_AB_eta_1D(Coefs,A_0,B_0,A_exact,B_exact);toc;
% figure;subplot(2,2,1);mesh(A_eta);subplot(2,2,2);mesh(B_eta);subplot(2,2,3);mesh(L1_eta);subplot(2,2,4);mesh(L2_eta);
% figure;subplot(1,2,1);meshc(A_eta-A_exact);subplot(1,2,2);meshc(B_eta-B_exact);
% calc_F_1D(A_eta,B_eta,L1_eta,L2_eta)
figure('Position', get(groot, 'Screensize'));hold on;
[X,Y]=meshgrid(linspace(0,1,N)',(1:100:P)'*delta_t);
p1=subplot(2,3,1);mesh(X,Y,A_eta(:,1:100:P)');title('A(t,x)','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);set(p1,'zlim',[0,30]);ylabel('time','Rotation',-40,'FontSize', 20);p1.FontSize = 20;colorbar; caxis([0 30]);
p2=subplot(2,3,2);mesh(X,Y,A_exact(:,1:100:P)');title('A(t,x) exact','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);set(p2,'zlim',[0,30]);ylabel('time','Rotation',-40,'FontSize', 20);p2.FontSize = 20;colorbar;caxis([0 30]);
p3=subplot(2,3,3);mesh(X,Y,L1_eta(:,1:100:P)');title('\lambda_1','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);ylabel('time','Rotation',-40,'FontSize', 20);p3.FontSize = 20;colorbar;
p4=subplot(2,3,4);mesh(X,Y,B_eta(:,1:100:P)');title('B(t,x)','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);set(p4,'zlim',[7.5,12]);ylabel('time','Rotation',-40,'FontSize', 20);p4.FontSize = 20;colorbar; caxis([7.5 12]);
p5=subplot(2,3,5);mesh(X,Y,B_exact(:,1:100:P)');title('B(t,x) exact','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);set(p5,'zlim',[7.5,12]);ylabel('time','Rotation',-40,'FontSize', 20);p5.FontSize = 20;colorbar; caxis([7.5 12]);
p6=subplot(2,3,6);mesh(X,Y,L2_eta(:,1:100:P)');title('\lambda_2','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);ylabel('time','Rotation',-40,'FontSize', 20);p6.FontSize = 20;colorbar;
    
%% Calcul de F en certains eta
% N_integer=10;
% d_values = linspace(lb(1),ub(1),N_integer);gamma_values = linspace(lb(2),ub(2),N_integer);
% F_values_1 = zeros(N_integer,N_integer);F_values_2 = zeros(N_integer,N_integer);
% for i=1:N_integer
%     disp(i);
%     for j=1:N_integer
%         Coefs = [d_values(i);gamma_values(j)];
%         Big_matrix_B=sparse(Coefs(1)*Big_matrix+eye(N));
%         [A_eta,B_eta,L1_eta,L2_eta] = calc_AB_eta_1D(Coefs,A_0,B_0,A_exact,B_exact);
%         [F_val] = calc_F_1D(A_eta,B_eta,L1_eta,L2_eta);
%         F_values_1(i,j) = F_val(1);
%         F_values_2(i,j) = F_val(2);
%     end
% end
% %figure;subplot(1,2,1);meshc(d_values,gamma_values,F_values_1);subplot(1,2,2);meshc(d_values,gamma_values,F_values_2);
% figure;subplot(1,2,1);surfl(d_values,gamma_values,F_values_1);subplot(1,2,2);surfl(d_values,gamma_values,F_values_2);
% %figure;subplot(2,1,1);plot(d_values,F_values_1(1,:));subplot(2,1,2);plot(d_values,F_values_2(1,:));

%%  Iteration en eta
%Coefs_0 = lb + rand(2,1).*(ub-lb);
A=[];b=[];All_coefs = Coefs_0;
Big_matrix_B=sparse(Coefs_0(1)*Big_matrix+eye(N));
[A_eta,B_eta,L1_eta,L2_eta] = calc_AB_eta_1D(Coefs_0,A_0,B_0,A_exact,B_exact);
[F_val] = calc_F_1D(A_eta,B_eta,L1_eta,L2_eta);
F = @(Coefs) F_eta_1D(Coefs,A_0,B_0,A_exact,B_exact);
[Mat_jac] = Jac_approx(Coefs_0,F,F_val,1e-6,2);
Diff_eta = -Mat_jac\F_val;
Coefs_new = Coefs_0 + Diff_eta;iteration=1;
Coefs_new = proj_ortho(Coefs_new,lb,ub,A,b);Diff_Coefs = Coefs_new-Coefs_0;
while((norm(Diff_Coefs)>1e-2 && iteration<=100)  || iteration <=2)
    Coefs = Coefs_new;All_coefs(:,iteration+1)=Coefs_new;disp(iteration);%disp(Coefs_new);
    Big_matrix_B=sparse(Coefs(1)*Big_matrix+eye(N));
    [A_eta,B_eta,L1_eta,L2_eta] = calc_AB_eta_1D(Coefs,A_0,B_0,A_exact,B_exact);
    [F_val] = calc_F_1D(A_eta,B_eta,L1_eta,L2_eta);
    [Mat_jac] = Jac_approx(Coefs,F,F_val,1e-6,2);
    Diff_eta = -Mat_jac\F_val;
    Coefs_new = Coefs + Diff_eta;iteration=iteration+1;
    Coefs_new = proj_ortho(Coefs_new,lb,ub,A,b);Diff_Coefs = Coefs_new-Coefs;
end
%disp(Coefs_new);
%% 
figure('Position', get(groot, 'Screensize'));hold on;
[X,Y]=meshgrid(linspace(0,1,N)',(1:100:P)'*delta_t);
p1=subplot(2,3,1);mesh(X,Y,A_eta(:,1:100:P)');title('A(t,x)','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);ylabel('time','Rotation',-40,'FontSize', 20);p1.FontSize = 20;colorbar;
p2=subplot(2,3,2);mesh(X,Y,A_exact(:,1:100:P)');title('A(t,x) exact','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);ylabel('time','Rotation',-40,'FontSize', 20);p2.FontSize = 20;colorbar;
p3=subplot(2,3,3);mesh(X,Y,L1_eta(:,1:100:P)');title('\lambda_1','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);ylabel('time','Rotation',-40,'FontSize', 20);p3.FontSize = 20;colorbar;
p4=subplot(2,3,4);mesh(X,Y,B_eta(:,1:100:P)');title('B(t,x)','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);ylabel('time','Rotation',-40,'FontSize', 20);p4.FontSize = 20;colorbar;
p5=subplot(2,3,5);mesh(X,Y,B_exact(:,1:100:P)');title('B(t,x) exact','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);ylabel('time','Rotation',-40,'FontSize', 20);p5.FontSize = 20;colorbar;
p6=subplot(2,3,6);mesh(X,Y,L2_eta(:,1:100:P)');title('\lambda_2','FontSize', 40);xlabel('space','Rotation',30,'FontSize', 20);ylabel('time','Rotation',-40,'FontSize', 20);p6.FontSize = 20;colorbar;
    

%% 
figure('Position', get(groot, 'Screensize'));hold on;
a1=subplot(1,2,1);xlim([0,iteration-1]);ylim([lb(1)-1,ub(1)+1]);grid on;
plot((0:iteration-1),All_coefs(1,:)','*',[0,iteration-1],[Coefs_exact(1),Coefs_exact(1)],'-',[0,iteration-1],[lb(1),lb(1)],'--',[0,iteration-1],[ub(1),ub(1)],'--','Color',[0 0.4470 0.7410],'Linewidth',3,'MarkerSize',15);
a1.FontSize = 20;a1.YGrid = 'on';
xlabel('Number of iteration','FontSize', 30);
title('Estimation of d','FontSize', 40);
legend('(d^{k})_{k\geq 0}','d_{exact}','lower bound','upper bound');
a2=subplot(1,2,2);xlim([0,iteration-1]);ylim([lb(2)-1,ub(2)+1]);grid on;
plot((0:iteration-1),All_coefs(2,:)','*',[0,iteration-1],[Coefs_exact(2),Coefs_exact(2)],'-',[0,iteration-1],[lb(2),lb(2)],'--',[0,iteration-1],[ub(2),ub(2)],'--','Color',[0.8500 0.3250 0.0980],'Linewidth',3,'MarkerSize',15);
a2.FontSize = 20;a2.YGrid = 'on';
xlabel('Number of iteration','FontSize', 30);
title('Estimation of \mu','FontSize', 40);
legend('(\mu^{k})_{k\geq 0}','\mu_{exact}','lower bound','upper bound');
%legend('d','\gamma');

%% Enregistrer les résultats
save('Psi1_estime.mat','A_eta');
save('Psi1_exact.mat','A_exact');
save('Psi2_estime.mat','B_eta');
save('Psi2_exact.mat','B_exact');
save('eta_estime.mat','All_coefs');
save('eta_exact.mat','Coefs_exact');
save('lambda1_estime.mat','L1_eta');
save('lambda2_estime.mat','L2_eta');

%% Analyse globale de sensibilité
%Par convention :   1 l'algo converge vers la bonne valeur
%                   0 l'algo ne converge pas en 100 itérations
%                  -1 l'algo converge mais pas vers la bonne valeur
% N_coefs = 100;
% All_eta = zeros(N_coefs,3);
% sz = zeros(N_coefs,1);
% for i=1:N_coefs
%     disp(i);
%     Coefs_0=lb + rand(2,1).*(ub-lb);
%     All_eta(i,1:2) = Coefs_0;
%     [Coefs_new,cvg] = calc_iter_eta_1D(Coefs_0,lb,ub,A,b,A_0,B_0,A_exact,B_exact);
%     if (cvg==1)%(norm(Coefs_new-Coefs_exact)<=1e-1)
%         All_eta(i,3) = 1;
%     elseif (cvg==0)
%         All_eta(i,3) = 0;
%     elseif (cvg == -1)
%         All_eta(i,3) = -1;
%     end
%     sz(i) = norm(Coefs_new-Coefs_exact);
% end
% %figure;scatter(All_eta(:,1),All_eta(:,2),[],All_eta(:,3),'filled');
% figure;scatter(All_eta(:,1),All_eta(:,2),[],sz,'filled');
