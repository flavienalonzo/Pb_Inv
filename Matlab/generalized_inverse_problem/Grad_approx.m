function [GradF] = Grad_approx(x,F,F_already_known,epsilon,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
GradF = zeros(n,1);
for i=1:n
    e_i = zeros(n,1);
    e_i(i,1) = 1;
    GradF(i,1) = (F(x+epsilon*e_i) - F_already_known)/epsilon;
end
end

