function [error] = fun_int_approx(X,Y,Nbr_var)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global delta_t h
if (Nbr_var==1)
    error = h*sum(X.*Y);
elseif (Nbr_var==2)
    Z = X.*Y;
    Z_int_t = delta_t*sum(Z,2);
    error = h*sum(Z_int_t);
elseif (Nbr_var==3)
    Zint= X(:).*Y(:);
    error = delta_t*h*sum(Zint);
end
end

