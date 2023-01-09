clearvars -except u_f c_f Coefs_new
%% 
global t0 tf delta_t N P x h diff_u diff_c chi Coefs_exact delta_mea Measures Cq C_eps W_eps C_a W_a;

t0=0e0; tf=7.1e0; delta_t=1e-2;   N=100; delta_mea = 100; h=1/(N-1);  P = round((tf-t0)/delta_t+1);
x = linspace(0,1,N)';
diff_u=5e-5; diff_c = 1e-3; chi=1e1;

Cq = (1e-2)^2*delta_t*min((1:P),(1:P)');
C_eps = sparse(diag(1e-1*ones(2*N,1)));    W_eps = inv(C_eps);
C_a = sparse(diag(1e-2*ones(2*N,1)));      W_a = inv(C_a);

lb=[log10(7e-2);log10(6e-2);log10(7e-2);log10(9e-3);log10(5e-2)]/log10(exp(1));%Contraintes sur theta : minoration
ub=[log10(6e-1);log10(5e-1);log10(4e-1);log10(6e-2);log10(2e-1)]/log10(exp(1));%Contraintes sur theta : majoration

choix = input('Voulez-vous commencer une nouvelle simulation (0) ou continuez la précédente (1) ? ');
if (choix==0)
    close all;
    a = 1/2;b = 7/10;c = (a+b)/2;
    u_0 = zeros(N,1);
    c_0 = zeros(N,1)+2e-1;
    for i=1:N
        if (abs(x(i)-c)<(b-a)/2)
            y = (x(i)-c)/((b-a)/2);
            u_0(i) = exp(1-1/(1-y^2));
        end
        c_0(i) = 0.5*u_0(i);%1e0-0.5*u_0(i);
    end
    Coefs_0 = log10([0.2864;0.0647;0.3075;0.0529;0.1281])/log10(exp(1));%(ub-lb).*rand(5,1)+lb;%log10([1e-1;5e-2;7e-2;6e-2;1e-1])/log10(exp(1));
    A_ineq=[-1 1 0 0 0; 0 0 -1 1 0];b_ineq=[-1e-8;-1e-8];
    while(~all(A_ineq*Coefs_0<=b_ineq))
        Coefs_0 = (ub-lb).*rand(5,1)+lb;
    end
elseif (choix==1)
     u_0 = u_f;
     c_0 = c_f;
     Coefs_0 = Coefs_new;
else
    disp('Mauvais choix')
end

Coefs_exact = log10([0.2;0.1;0.1;0.03;0.08])/log10(exp(1));
%Coefs_0 = Coefs_exact;
Measures = measurements(u_0,c_0,Coefs_exact);

disp(exp(Coefs_0'));
%Resolution(Coefs_0(1),Coefs_0(2),Coefs_0(3),Coefs_0(4),Coefs_0(5),u_0,c_0,4);

%% 
% lb=[log10(2e-1);log10(1e-1);log10(1e-1);log10(3e-2);log10(3e-2)]/log10(exp(1));%Contraintes sur theta : minoration
% ub=[log10(2e-1);log10(1e-1);log10(1e-1);log10(3e-2);log10(4e-1)]/log10(exp(1));%Contraintes sur theta : majoration
% Coefs_0 = (ub-lb).*rand(5,1)+lb;%log10([1e-1;5e-2;7e-2;6e-2;1e-1])/log10(exp(1));
% A_ineq=[-1 1 0 0 0; 0 0 -1 1 0];b_ineq=[-1e-8;-1e-8];
% while(~all(A_ineq*Coefs_0<=b_ineq))
%     Coefs_0 = (ub-lb).*rand(5,1)+lb;
% end
% itermax = 20;
% save_file = zeros(itermax,11);
% couleur = zeros(itermax,3);
% for iter=1:itermax
%     donnee = zeros(1,10);
%     Coefs_0 = (ub-lb).*rand(5,1)+lb;%log10([1e-1;5e-2;7e-2;6e-2;1e-1])/log10(exp(1));
%     A_ineq=[-1 1 0 0 0; 0 0 -1 1 0];b_ineq=[-1e-8;-1e-8];
%     while(~all(A_ineq*Coefs_0<=b_ineq))
%         Coefs_0 = (ub-lb).*rand(5,1)+lb;
%     end
%     Coefs = Coefs_0;
%     donnee(1,1:5) = Coefs_0;
%     [error_total] = [1e5;1e5;1e3;1e5;1e4].*Resolution(Coefs(1),Coefs(2),Coefs(3),Coefs(4),Coefs(5),u_0,c_0,2);
%     donnee(1,6:10) = error_total;%
%     donnee(1,11) = norm(error_total);
%     couleur(iter,1:3) = [0,0,1];
%     save_file(iter,:) = donnee;
% end
% figure;%bubblechart(exp(save_file(:,2)),exp(save_file(:,4)),save_file(:,6));
% subplot(2,3,1);bubblechart(exp(save_file(:,5)),save_file(:,6),save_file(:,11));
% subplot(2,3,2);bubblechart(exp(save_file(:,5)),save_file(:,7),save_file(:,11));
% subplot(2,3,3);bubblechart(exp(save_file(:,5)),save_file(:,8),save_file(:,11));
% subplot(2,3,4);bubblechart(exp(save_file(:,5)),save_file(:,9),save_file(:,11));
% subplot(2,3,5);bubblechart(exp(save_file(:,5)),save_file(:,10),save_file(:,11));


%% Résolution par utilisation d'algo génétique

% tic;
% fun2 = @(x) norm([1e5;1e5;1e3;1e5;1e4].*Resolution(x(1),x(2),x(3),x(4),x(5),u_0,c_0,2),2);
% Aeq=[];beq=[]; 
% nonlcon=[]; 
% nvars=5;
% options = optimoptions('ga','PlotFcn', @gaplotbestf,'Display','iter','InitialPopulationMatrix',[Coefs_0(1),Coefs_0(2),Coefs_0(3),Coefs_0(4),Coefs_0(5)]);%,'FunctionTolerance',1e-10);%,'InitialPopulation',Coefs_0','PopulationSize',3
% [Coefs,fval,exitflag,output,population,scores] = ga(fun2,nvars,A_ineq,b_ineq,Aeq,beq,lb,ub,nonlcon,options);
% Resolution(Coefs(1),Coefs(2),Coefs(3),Coefs(4),Coefs(5),u_0,c_0,4);
% Coefs_new = Coefs;
% Psi_f = Resolution(Coefs(1),Coefs(2),Coefs(3),Coefs(4),Coefs(5),u_0,c_0,4);
% u_f = Psi_f(:,1); c_f = Psi_f(:,2);
% toc;

%% Résolution par méthode de Newton-Raphson
tol_NR = 1e-3; iter=0;itermax=2;
tic;
Coefs = Coefs_0;
Coefs_all = zeros(5,itermax+1);
error_total = zeros(2,itermax);
%[Int_approx] = [1e5;1e5;1e3;1e5;1e4].*Resolution(Coefs(1),Coefs(2),Coefs(3),Coefs(4),Coefs(5),u_0,c_0,2);
[result] = Resolution(Coefs(1),Coefs(2),Coefs(3),Coefs(4),Coefs(5),u_0,c_0,3);
Int_approx = result(2:end);%[1e5;1e5;1e3;1e5;1e4].*result(2:end);
Coefs_all(:,1) = Coefs_0;
Resolution(Coefs_0(1),Coefs_0(2),Coefs_0(3),Coefs_0(4),Coefs_0(5),u_0,c_0,4);
fun_int_ap = @(x) Resolution(x(1),x(2),x(3),x(4),x(5),u_0,c_0,2);%[1e5;1e5;1e3;1e5;1e4].*Resolution(x(1),x(2),x(3),x(4),x(5),u_0,c_0,2);
%[M] = Jac_approx(Coefs,fun_int_ap,Int_approx,1e-3,5);
%Difference_theta = linsolve(M,-Int_approx);
Coefs_new = Coefs ;%+ tol_NR+1e0;%proj_ortho(Coefs + Difference_theta,lb,ub,A_ineq,b_ineq);
while ((norm(Coefs_new-Coefs)>tol_NR && iter <itermax)||(iter==0))
    if(iter>0)
        Coefs_all(:,iter+1) = Coefs_new;
        error_total(1,iter)=result(1);error_total(2,iter)=norm(Int_approx);disp(iter);
    end
    iter = iter+1;Coefs = Coefs_new;
    [result] = Resolution(Coefs(1),Coefs(2),Coefs(3),Coefs(4),Coefs(5),u_0,c_0,3);
    Int_approx = result(2:end);%[1e5;1e5;1e3;1e5;1e4].*result(2:end);
    % Méthode de Newton-Raphson    
    %On cherche les zéros d'une fonction
    Point_connu = fun_int_ap(Coefs);
    [Matrix] = Jac_approx(Coefs,fun_int_ap,Point_connu,1e-8,5);
    Difference_theta = linsolve(Matrix,-Point_connu);
    Coefs_new = Coefs + Difference_theta;
    Coefs_new = proj_ortho(Coefs_new,lb,ub,A_ineq,b_ineq);
end
Coefs_all(:,iter+1) = Coefs_new;error_total(1,iter)=result(1);
error_total(2,iter)=norm(Int_approx);
%disp(iter);disp(Coefs_new);
Psi_f = Resolution(Coefs(1),Coefs(2),Coefs(3),Coefs(4),Coefs(5),u_0,c_0,4);
u_f = Psi_f(:,1); c_f = Psi_f(:,2);
%CrossProduct = Resolution(Coefs(1),Coefs(2),Coefs(3),Coefs(4),Coefs(5),u_0,c_0,5); 
toc;
%% 

colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
figure('Position', get(groot, 'Screensize'));
slg1=subplot(2,3,1);%hold on;
semilogy((0:iter),exp(Coefs_all(1,1:iter+1)'),'*',[0,iter],[exp(Coefs_exact(1)),exp(Coefs_exact(1))],'-',[0,iter],[exp(lb(1)),exp(lb(1))],'--',[0,iter],[exp(ub(1)),exp(ub(1))],'--','Color',colors(1,:),'Linewidth',3,'MarkerSize',10);
slg1.FontSize = 20;slg1.YGrid = 'on';
title('Estimation of \rho','FontSize', 40);
legend('(\rho^k)_{k\geq 0}','\rho_{exact}','lower bound','upper bound');

slg2=subplot(2,3,2); %hold on;
semilogy((0:iter),exp(Coefs_all(2,1:iter+1)'),'*',[0,iter],[exp(Coefs_exact(2)),exp(Coefs_exact(2))],'-',[0,iter],[exp(lb(2)),exp(lb(2))],'--',[0,iter],[exp(ub(2)),exp(ub(2))],'--','Color',colors(2,:),'Linewidth',3,'MarkerSize',10);
slg2.FontSize = 20;slg2.YGrid = 'on';
title('Estimation of \delta','FontSize', 40);
legend('(\delta^k)_{k\geq 0}','\delta_{exact}','lower bound','upper bound');

% subplot(2,3,3);% hold on;
% semilogy((0:iter-1),error_total(1,1:iter),(0:iter-1),error_total(2,1:iter),'DisplayName','error_total','LineWidth',3);
% title('Errors');
% legend('e_{total}','||F(\theta,\Psi,\lambda)||');grid on;

slg4=subplot(2,3,4); %hold on;
semilogy((0:iter),exp(Coefs_all(3,1:iter+1)'),'*',[0,iter],[exp(Coefs_exact(3)),exp(Coefs_exact(3))],'-',[0,iter],[exp(lb(3)),exp(lb(3))],'--',[0,iter],[exp(ub(3)),exp(ub(3))],'--','Color',colors(3,:),'Linewidth',3,'MarkerSize',10);%
slg4.FontSize = 20;slg4.YGrid = 'on';
title('Estimation of \alpha','FontSize', 40);
legend('(\alpha^k)_{k\geq 0}','\alpha_{exact}','lower bound','upper bound');

slg5=subplot(2,3,5); %hold on;
semilogy((0:iter),exp(Coefs_all(4,1:iter+1)'),'*',[0,iter],[exp(Coefs_exact(4)),exp(Coefs_exact(4))],'-',[0,iter],[exp(lb(4)),exp(lb(4))],'--',[0,iter],[exp(ub(4)),exp(ub(4))],'--','Color',colors(4,:),'Linewidth',3,'MarkerSize',10);
slg5.FontSize = 20;slg5.YGrid = 'on';
title('Estimation of \beta','FontSize', 40);
legend('(\beta^k)_{k\geq 0}','\beta_{exact}','lower bound','upper bound');

slg6=subplot(2,3,6);% hold on;
semilogy((0:iter),exp(Coefs_all(5,1:iter+1)'),'*',[0,iter],[exp(Coefs_exact(5)),exp(Coefs_exact(5))],'-',[0,iter],[exp(lb(5)),exp(lb(5))],'--',[0,iter],[exp(ub(5)),exp(ub(5))],'--','Color',colors(5,:),'Linewidth',3,'MarkerSize',10);
slg6.FontSize = 20;slg6.YGrid = 'on';
title('Estimation of \gamma','FontSize', 40);
legend('(\gamma^k)_{k\geq 0}','\gamma_{exact}','lower bound','upper bound');

