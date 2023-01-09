function [M] = Jac_approx(x,F,F_already_known,epsilon,n)
%Calculate an approximation of the Jacobian matrix of a function that is
%time-consuming to get the value of one point.
% F : R^n -> R^n with Jac(f)(x) = M(x) in R^n*n
%M(x)_ij = ei^T*(F(x+h*ej)-F(x))/h
%We need to know the point x, the function F, the value F(x) and the scale h
M = zeros(n,n);%pas=randn(1)*epsilon;
vec_zero = zeros(n,1);
for i=1:n
    e_i = vec_zero;
    e_i(i,1) = 1;
    x_i = x+epsilon*e_i;%pas*e_i;
    F_i = F(x_i);
    for j=1:n
        e_j = vec_zero;
        e_j(j,1) = 1;
        M(j,i) = e_j'*(F_i-F_already_known)/epsilon;%pas;
    end
end

end

%Test 
% f =@(x)[x(1)^2*x(2);5*x(1)+sin(x(2))];disp(Jac_approx([0.5;1],f,f([0.5;1]),0.001,2));disp([1 0.25;5 cos(1)]);