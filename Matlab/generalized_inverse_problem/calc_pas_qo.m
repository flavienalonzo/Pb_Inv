function [tau,iter] = calc_pas_qo(c1,c2,xn,hn,F)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    tau = 1e1; T = 1e20;iter = 1;
    F_xn = F(xn);
    [GradF_xn] = Grad_approx(xn,F,F_xn,1e-8,5);
    y = xn + tau*hn;
    F_y = F(y);
    [GradF_y] = Grad_approx(y,F,F_y,1e-8,5); 
    Y = xn + T*hn;
    F_Y = F(Y);
    if ~(hn'*GradF_xn<0)
        disp('La direction hn n_est pas correcte');
    elseif (F_Y<=F_xn+c1*T*hn'*GradF_xn)
        disp('T n_est pas assez grand');
    elseif ~(F_y<=F_xn+c1*tau*hn'*GradF_xn)
        disp('tau n_est pas assez petit');
    end
    cond_Wolf = (hn'*GradF_y>=c2*hn'*GradF_xn);
    while ~cond_Wolf
        tau_n = exp((log(tau)+log(T))/2);
        cond_Armijo = (F_y<=F_xn+c1*tau_n*hn'*GradF_xn);
        if cond_Armijo
            tau = tau_n;
        else
            T = tau_n;
        end
        disp([tau,T]);
        y = xn + tau*hn;
        F_y = F(y);
        [GradF_y] = Grad_approx(y,F,F_y,1e-8,5);     
        cond_Wolf = (hn'*GradF_y>=c2*hn'*GradF_xn);
        iter = iter+1;
    end
end

