function [Coefs_proj] = proj_ortho(Coefs,lb,ub,A,b)
%On calcule le projeté orthogonale d'un vecteur Coefs 
%sachant que l'espace K est {x : lb<x<ub, Ax<=b}
Aeq=[];beq=[];
fun = @(x)fun_ortho(x,Coefs);
options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',true);
[Coefs_proj,fval,exitflag,output]= fmincon(fun,Coefs,A,b,Aeq,beq,lb,ub,[],options); 
if (exitflag~=1)
    disp('Problem with proj_ortho');
end
end