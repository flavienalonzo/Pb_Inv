function [f,g] = fun_ortho(x,Coefs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
f = 0.5*norm(x-Coefs)^2;
g = x-Coefs;
end

