function [Coefs_new,cvg] = calc_iter_eta_1D(Coefs_0,lb,ub,A,b,A_0_ex,B_0_ex,A_exact,B_exact)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global N Big_matrix Big_matrix_B
All_coefs = Coefs_0;
Big_matrix_B=sparse(Coefs_0(1)*Big_matrix+eye(N));
[A_eta,B_eta,L1_eta,L2_eta] = calc_AB_eta_1D(Coefs_0,A_0_ex,B_0_ex,A_exact,B_exact);
[F_val] = calc_F_1D(A_eta,B_eta,L1_eta,L2_eta);
F = @(Coefs) F_eta_1D(Coefs,A_0_ex,B_0_ex,A_exact,B_exact);
[Mat_jac] = Jac_approx(Coefs_0,F,F_val,1e-6,2);
Diff_eta = -Mat_jac\F_val;
Coefs_new = Coefs_0 + Diff_eta;iteration=1;
Coefs_new = proj_ortho(Coefs_new,lb,ub,A,b);Diff_Coefs = Coefs_new-Coefs_0;
cvg = 1;
while((norm(Diff_Coefs)>1e-2 && iteration<=100)  || iteration <=2)
    Coefs = Coefs_new;All_coefs(:,iteration+1)=Coefs_new;
    Big_matrix_B=sparse(Coefs(1)*Big_matrix+eye(N));
    [A_eta,B_eta,L1_eta,L2_eta] = calc_AB_eta_1D(Coefs,A_0_ex,B_0_ex,A_exact,B_exact);
    [F_val] = calc_F_1D(A_eta,B_eta,L1_eta,L2_eta);
    [Mat_jac] = Jac_approx(Coefs,F,F_val,1e-6,2);
    Diff_eta = -Mat_jac\F_val;
    Coefs_new = Coefs + Diff_eta;iteration=iteration+1;
    Coefs_new = proj_ortho(Coefs_new,lb,ub,A,b);Diff_Coefs = Coefs_new-Coefs;
end
if (~(iteration<=100))
    cvg = 0;
end
if (~(Coefs_new==Coefs + Diff_eta))
    cvg = -1;
end
end

