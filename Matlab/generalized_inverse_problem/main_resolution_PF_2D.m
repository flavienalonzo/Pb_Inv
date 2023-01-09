clear all;close all;
global t0 tf delta_t N M P x_space y_space h Coef_a Coef_b Coef_alp Coef_rho Coef_K Coefs_exact delta_mea Big_matrix Big_matrix_A Big_matrix_B W_eps Ca;
tic;
t0=0e0; tf=1.0e1; delta_t=1e-3; Lx = 1; Ly = 1/2; N=10; delta_mea = 100; h=1/(N-1); M=round(Ly/h+1); P = round((tf-t0)/delta_t+1);
x_space = linspace(0,Lx,N)'; y_space = linspace(0,Ly,M)'; 
Coef_a=92; Coef_b=64; Coef_alp=1.5; Coef_rho=18.5; Coef_K=0.1;
Coefs_exact = [20;50];%(d,gamma)
lb=[19;49];ub = [21;51];
W_eps = 1e-6;
Ca = (1e-3)^2*eye(N*M);
f = @(u,v) Coef_a-u-Coef_rho*u*v/(1+u+Coef_K*u^2);
g = @(u,v) Coef_alp*(Coef_b-v)-Coef_rho*u*v/(1+u+Coef_K*u^2);
h_u = @(u,v) (Coef_rho*v*(1+u+Coef_K*u^2)-Coef_rho*u*v*(1+2*Coef_K*u))/(1+u+Coef_K*u^2)^2;
h_v = @(u,v) Coef_rho*u/(1+u+Coef_K*u^2);
f_u = @(u,v) -1-h_u(u,v);
f_v = @(u,v) -h_v(u,v);
g_u = @(u,v) -h_u(u,v);
g_v = @(u,v) -Coef_alp-h_v(u,v);

fun = @(x) abs(f(x(1),x(2)))+abs(g(x(1),x(2)));
options = optimset('TolFun',1e-3,'TolX',1e-2);
x = fminsearch(fun,[10,9],options);
[Y_space,X_space]=meshgrid(y_space,x_space);
A_0 = x(1)*ones(N,M);%+0.0001*randn(N,M);
B_0 = x(2)*ones(N,M);%+0.0001*randn(N,M);

% disp(f_u(x(1),x(2))+g_v(x(1),x(2))<0);
% disp(f_u(x(1),x(2))*g_v(x(1),x(2))-f_v(x(1),x(2))*g_u(x(1),x(2))>0);
% disp(Coefs_exact(1)*f_u(x(1),x(2))+g_v(x(1),x(2))>0);
% disp((Coefs_exact(1)*f_u(x(1),x(2))+g_v(x(1),x(2)))^2-4*Coefs_exact(1)*(f_u(x(1),x(2))*g_v(x(1),x(2))-f_v(x(1),x(2))*g_u(x(1),x(2)))>0);

options = optimset('TolFun',1e-10,'TolX',1e-10);

x_ex = fminsearch(fun,x,options);
%A_0_ex = x_ex(1)*ones(N,M);
%B_0_ex = x_ex(2)*ones(N,M);

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
[liste_couples,list_wave] = find_nm_couples(k1,k2,Lx,Ly);
subplot(2,1,2);hold on;fplot(fun_wavenumber,[0,3/2*k2^2]);plot(k1^2,0,'r*');plot(k2^2,0,'r*');plot(list_wave,fun_wavenumber(list_wave),'b*');xlabel('k^2');ylabel('Re \lambda');
[l_mode,modes] = max(fun_wavenumber(list_wave));
nm_modes=liste_couples((randi(size(list_wave))),:);
%nm_modes=liste_couples(modes,:);
A_0_ex = x_ex(1)*ones(N,M)+0.01*cos(nm_modes(1)*pi*x_space/Lx).*cos(nm_modes(2)*pi*y_space'/Ly);
B_0_ex = x_ex(2)*ones(N,M)+0.001*randn(N,M);%cos(n_mode*pi*x_space);%+0.15*cos(pi*x_space);
if (lb(1)<d_crit)
    disp('Il y a des valeurs de d qui ne produisent pas d_instabilité');
end


%% Pré-calcul des matrices pour ne pas le faire dans la boucle en temps
Big_matrix = delta_t*(4*eye(M*N)-diag(ones(N*M-1,1),1)-diag(ones(N*M-1,1),-1)-diag(ones(N*(M-1),1),N)-diag(ones(N*(M-1),1),-N))/h^2;
index_0 = [(N:N:(M-1)*N),(N+1:N:(M-1)*N+1);(N+1:N:(M-1)*N+1),(N:N:(M-1)*N)];
for i=1:size(index_0,2)
    Big_matrix(index_0(1,i),index_0(2,i)) = 0;
end
index_1 = [(1:N:(M-1)*N+1),(N:N:M*N);(N:N:M*N),(1:N:(M-1)*N+1)];
for i=1:size(index_1,2)
    Big_matrix(index_1(1,i),index_1(2,i)) = -delta_t/h^2;
end
index_2 = [(1:N),((M-1)*N+1:M*N);(1:N),((M-1)*N+1:M*N)];%(1:N:(M-1)*N+1),(N:N:M*N)
for i=1:size(index_2,2)
    Big_matrix(index_2(1,i),index_2(2,i)) = Big_matrix(index_2(1,i),index_2(2,i))-delta_t/h^2;
end
Big_matrix_A=sparse(Big_matrix+eye(N*M));
Big_matrix_B=sparse(Coefs_exact(1)*Big_matrix+eye(N*M));
Spectre_A = eig(Big_matrix_A);Spectre_B = eig(Big_matrix_B);
%% Calcul de la solution avec les paramètres de référence
tic;
[A_exact,B_exact] = calc_AB(Coefs_exact,A_0_ex,B_0_ex,zeros(N,M,P),zeros(N,M,P));%figure;subplot(1,5,1);imagesc(A_exact(:,:,P)>mean(A_exact(:,:,P),'all'));subplot(1,5,[2,3]);surf(A_exact(:,:,P));subplot(1,5,[4,5]);surf(B_exact(:,:,P));
toc;
tic;
[A_exact_non,B_exact_non] = calc_AB(Coefs_exact+1,A_0_ex,B_0_ex,zeros(N,M,P),zeros(N,M,P));%figure;subplot(1,5,1);imagesc(A_exact_non(:,:,P)>mean(A_exact_non(:,:,P),'all'));subplot(1,5,[2,3]);surf(A_exact_non(:,:,P));subplot(1,5,[4,5]);surf(B_exact_non(:,:,P));
toc;
%figure;surf(A_exact(:,:,P)-A_exact_non(:,:,P));
% tic;
% delta_f = 2000;
% [A_long,B_long] = calc_AB_long(Coefs_exact,A_0,B_0,delta_f);
% figure;subplot(1,5,1);imagesc(A_long(:,:,P)>mean(A_long(:,:,P),'all'));subplot(1,5,[2,3]);surf(A_long(:,:,P));subplot(1,5,[4,5]);surf(B_long(:,:,P));
% toc;
%% calcul de (A,B,L1,L2) pour une valeur de eta fixée
Coefs = [18,52];
Big_matrix_B=sparse(Coefs(1)*Big_matrix+eye(N*M));
[A_eta,B_eta,L1_eta,L2_eta] = calc_AB_eta(Coefs,A_0_ex,B_0_ex,A_exact,B_exact);
figure;subplot(1,5,1);imagesc(A_eta(:,:,P)>mean(A_eta(:,:,P),'all'));subplot(1,5,[2,3]);surf(A_eta(:,:,P));subplot(1,5,[4,5]);surf(B_eta(:,:,P));
%% Calcul de F(A,B,L1,L2)
[F_val] = calc_F(A_eta,B_eta,L1_eta,L2_eta);

%% Calcul d'une approximation de la matrice jacobienne de F en Coefs

% F = @(Coefs) F_eta(Coefs,A_0_ex,B_0_ex,A_exact,B_exact);
% [Mat_jac] = Jac_approx(Coefs,F,F_val,1e-6,2);
% Diff_eta = -Mat_jac\F_val;
% Coefs_new = Coefs + Diff_eta;

%% Iteration en eta
% Coefs_0 = [19;51];A=[];b=[];All_coefs = Coefs_0;
% Big_matrix_B=sparse(Coefs_0(1)*Big_matrix+eye(N*M));
% [A_eta,B_eta,L1_eta,L2_eta] = calc_AB_eta(Coefs_0,A_0_ex,B_0_ex,A_exact,B_exact);
% [F_val] = calc_F(A_eta,B_eta,L1_eta,L2_eta);
% F = @(Coefs) F_eta(Coefs,A_0_ex,B_0_ex,A_exact,B_exact);
% [Mat_jac] = Jac_approx(Coefs_0,F,F_val,1e-6,2);
% Diff_eta = -Mat_jac\F_val;
% Coefs_new = Coefs_0 + Diff_eta;iteration=1;
% Coefs_new = proj_ortho(Coefs_new,lb,ub,A,b);Diff_Coefs = Coefs_new-Coefs_0;
% while(norm(Diff_Coefs)>1e-2  || iteration <=2)
%     Coefs = Coefs_new;All_coefs(:,iteration+1)=Coefs_new;disp(iteration);disp(Coefs_new);
%     Big_matrix_B=sparse(Coefs(1)*Big_matrix+eye(N*M));
%     [A_eta,B_eta,L1_eta,L2_eta] = calc_AB_eta(Coefs,A_0_ex,B_0_ex,A_exact,B_exact);
%     [F_val] = calc_F(A_eta,B_eta,L1_eta,L2_eta);
%     [Mat_jac] = Jac_approx(Coefs,F,F_val,1e-8,2);
%     Diff_eta = -Mat_jac\F_val;
%     Coefs_new = Coefs + Diff_eta;iteration=iteration+1;
%     Coefs_new = proj_ortho(Coefs_new,lb,ub,A,b);Diff_Coefs = Coefs_new-Coefs;
% end
% disp(Coefs_new);
% figure;plot(All_coefs');legend('d','\gamma');

%% Résolution par algo génétique
% Coefs_0 = [13.0;5.0];A_ineq=[];b_ineq=[];
% tic;
% fun2 = @(x) sqrt(norm(F_eta(x,A_0_ex,B_0_ex,A_exact,B_exact),2));
% Aeq=[];beq=[]; 
% nonlcon=[]; 
% nvars=2;
% options = optimoptions('ga','PlotFcn', @gaplotbestf,'Display','iter','InitialPopulationMatrix',[Coefs_0(1),Coefs_0(2)]);
% [Coefs_ga,fval,exitflag,output,population,scores] = ga(fun2,nvars,A_ineq,b_ineq,Aeq,beq,lb,ub,nonlcon,options);
% toc;
% 
% [A_ga,B_ga] = calc_AB(Coefs_ga,A_0_ex,B_0_ex,zeros(N,M,P),zeros(N,M,P));figure;subplot(1,5,1);imagesc(A_ga(:,:,P)>mean(A_ga(:,:,P),'all'));subplot(1,5,[2,3]);surf(A_ga(:,:,P));subplot(1,5,[4,5]);surf(B_ga(:,:,P));


%% Calcul de F en certains eta
% N_integer=101;lb=[10-1e-7;14.995];ub = [10+1e-7;15.005];
% d_values = linspace(lb(1),ub(1),N_integer);gamma_values = linspace(lb(2),ub(2),N_integer);
% F_values_1 = zeros(N_integer,N_integer);F_values_2 = zeros(N_integer,N_integer);
% for i=1:1%N_integer
%     %disp(i);
%     for j=1:N_integer
%         disp(j);
%         Coefs = [d_values(j),Coefs_exact(2)];%[d_values(i);gamma_values(j)];
%         Big_matrix_B=sparse(Coefs(1)*Big_matrix+eye(N*M));
%         [A_eta,B_eta,L1_eta,L2_eta] = calc_AB_eta(Coefs,A_0_ex,B_0_ex,A_exact,B_exact);
%         [F_val] = calc_F(A_eta,B_eta,L1_eta,L2_eta);
%         F_values_1(i,j) = F_val(1);
%         F_values_2(i,j) = F_val(2);
%     end
% end
% %figure;subplot(1,2,1);meshc(d_values,gamma_values,F_values_1);subplot(1,2,2);meshc(d_values,gamma_values,F_values_2);
% %figure;subplot(1,2,1);surfl(d_values,gamma_values,F_values_1);subplot(1,2,2);surfl(d_values,gamma_values,F_values_2);
% figure;subplot(2,1,1);plot(d_values,F_values_1(1,:));subplot(2,1,2);plot(d_values,F_values_2(1,:));

%% Erreurs
% Coefs = [11,16];
% [A_eta,B_eta,L1_eta,L2_eta] = calc_AB_eta(Coefs,A_0_ex,B_0_ex,A_exact,B_exact);figure;subplot(2,5,[1,6]);imagesc(A_eta(:,:,P)>mean(A_eta(:,:,P),'all'));subplot(2,5,[2,3]);surf(A_eta(:,:,P));subplot(2,5,[4,5]);surf(B_eta(:,:,P));subplot(2,5,[7,8]);surf(L1_eta(:,:,1));subplot(2,5,[9,10]);surf(L2_eta(:,:,1));
% [e1,e2,e3] = calc_errors(A_eta,B_eta,L1_eta,L2_eta,A_exact,B_exact);
% disp(e1);disp(e2);disp(e3);
%% 
% Coefs = [11,16];
% [A_eta,B_eta,L1_eta,L2_eta] = calc_AB_eta(Coefs,A_0_ex,B_0_ex,A_exact,B_exact);figure;subplot(2,5,[1,6]);imagesc(A_eta(:,:,P)>mean(A_eta(:,:,P),'all'));subplot(2,5,[2,3]);surf(A_eta(:,:,P));subplot(2,5,[4,5]);surf(B_eta(:,:,P));subplot(2,5,[7,8]);surf(L1_eta(:,:,1));subplot(2,5,[9,10]);surf(L2_eta(:,:,1));
% disp(calc_F(A_eta,B_eta,L1_eta,L2_eta));
%% 
% tic;
% Corr = W_eps.*([A_exact(:);B_exact(:)]-[Tout_A(:);Tout_B(:)]);
% Corr1 = reshape(Corr(1:N*M*P),[N,M,P]);Corr2 = reshape(Corr(N*M*P+1:2*N*M*P),[N,M,P]);
% [Tout_L1,Tout_L2] = calc_adjoints_AB(Coefs,Tout_A,Tout_B,Corr1,Corr2);
% toc;
%% 
% tic;
% [Cross_prod1,Cross_prod2] = calc_cross_prod(Tout_L1,Tout_L2);
% [Relev_A,Relev_B] = calc_relev(Tout_L1(:,:,1),Tout_L2(:,:,1));
% toc;
%% 
% figure;data=A_exact;
% s = surf(data(:,:,1));
% for ii=1:10:P
%     s.ZData = data(:,:,ii);
% 
%     pause(0.05)
% end
