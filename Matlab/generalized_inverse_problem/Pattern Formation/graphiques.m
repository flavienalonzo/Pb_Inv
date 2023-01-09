clear all;
% load('eta_estime.mat');             %All_coefs
% load('eta_exact.mat');              %Coefs_exact
% load('lambda1_estime.mat');         %L1_eta
% load('lambda2_estime.mat');         %L2_eta
% load('Psi1_estime.mat');            %A_eta
% load('Psi2_estime.mat');            %B_eta
% load('Psi1_exact.mat');             %A_exact
% load('Psi2_exact.mat');             %B_exact
% 
t0=0e0; tf=2.0e0; delta_t=1e-3; N=100; delta_mea = 100; h=1/(N-1); P = round((tf-t0)/delta_t+1);
% %% (A-A_{ex})_{L^{\infty}(t_0,t_f,L^2(\Omega))}
% 
% figure;plot((t0:delta_t:tf),sqrt(h*sum((A_eta-A_exact).^2,1)));title('||A-A_{ex}||_{L^2(\Omega)}');
% ecart_A = max(sqrt(h*sum((A_eta-A_exact).^2,1)));
% disp('ecart_A');
% disp(ecart_A);
% 
% %% (B-B_{ex})_{L^{\infty}(t_0,t_f,L^2(\Omega))}
% 
% figure;plot((t0:delta_t:tf),sqrt(h*sum((B_eta-B_exact).^2,1)));title('||B-B_{ex}||_{L^2(\Omega)}');
% ecart_B = max(sqrt(h*sum((B_eta-B_exact).^2,1)));
% disp('ecart_B');
% disp(ecart_B);
% 
% %% ||L1||_{L^{\infty}(t_0,t_f,L^2(\Omega))}
% 
% figure;plot((t0:delta_t:tf),log(sqrt(h*sum((L1_eta).^2,1))));title('log(||\lambda_1||_{L^2(\Omega)})');
% norm_L1 = max(sqrt(h*sum((L1_eta).^2,1)));
% disp('norm_L1');
% disp(norm_L1);
% 
% %% ||L2||_{L^{\infty}(t_0,t_f,L^2(\Omega))}
% 
% figure;plot((t0:delta_t:tf),log(sqrt(h*sum((L2_eta).^2,1))));title('log(||\lambda_2||_{L^2(\Omega)})');
% norm_L2 = max(sqrt(h*sum((L2_eta).^2,1)));
% disp('norm_L2');
% disp(norm_L2);
% 
% %% N_iter et ||theta-theta_ex||_{\infty}
% 
% N_iter = size(All_coefs,2);
% disp('N_iter');
% disp(N_iter);
% ecart_theta = norm(All_coefs(:,N_iter)-Coefs_exact,Inf);
% disp('ecart_theta');
% disp(ecart_theta);
% 

%% Chargement de toutes les données
Coefs_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/001/eta_exact.mat').Coefs_exact;

data_001_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/001/eta_estime.mat').All_coefs;
data_001_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/001/lambda1_estime.mat').L1_eta;
data_001_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/001/lambda2_estime.mat').L2_eta;
data_001_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/001/Psi1_estime.mat').A_eta;
data_001_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/001/Psi2_estime.mat').B_eta;
data_001_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/001/Psi1_exact.mat').A_exact;
data_001_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/001/Psi2_exact.mat').B_exact;

data_002_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/002/eta_estime.mat').All_coefs;
data_002_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/002/lambda1_estime.mat').L1_eta;
data_002_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/002/lambda2_estime.mat').L2_eta;
data_002_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/002/Psi1_estime.mat').A_eta;
data_002_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/002/Psi2_estime.mat').B_eta;
data_002_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/002/Psi1_exact.mat').A_exact;
data_002_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/002/Psi2_exact.mat').B_exact;

data_003_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/003/eta_estime.mat').All_coefs;
data_003_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/003/lambda1_estime.mat').L1_eta;
data_003_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/003/lambda2_estime.mat').L2_eta;
data_003_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/003/Psi1_estime.mat').A_eta;
data_003_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/003/Psi2_estime.mat').B_eta;
data_003_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/003/Psi1_exact.mat').A_exact;
data_003_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/003/Psi2_exact.mat').B_exact;

data_004_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/004/eta_estime.mat').All_coefs;
data_004_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/004/lambda1_estime.mat').L1_eta;
data_004_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/004/lambda2_estime.mat').L2_eta;
data_004_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/004/Psi1_estime.mat').A_eta;
data_004_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/004/Psi2_estime.mat').B_eta;
data_004_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/004/Psi1_exact.mat').A_exact;
data_004_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/004/Psi2_exact.mat').B_exact;

data_005_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/005/eta_estime.mat').All_coefs;
data_005_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/005/lambda1_estime.mat').L1_eta;
data_005_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/005/lambda2_estime.mat').L2_eta;
data_005_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/005/Psi1_estime.mat').A_eta;
data_005_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/005/Psi2_estime.mat').B_eta;
data_005_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/005/Psi1_exact.mat').A_exact;
data_005_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/005/Psi2_exact.mat').B_exact;

data_006_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/006/eta_estime.mat').All_coefs;
data_006_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/006/lambda1_estime.mat').L1_eta;
data_006_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/006/lambda2_estime.mat').L2_eta;
data_006_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/006/Psi1_estime.mat').A_eta;
data_006_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/006/Psi2_estime.mat').B_eta;
data_006_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/006/Psi1_exact.mat').A_exact;
data_006_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/006/Psi2_exact.mat').B_exact;

data_007_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/007/eta_estime.mat').All_coefs;
data_007_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/007/lambda1_estime.mat').L1_eta;
data_007_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/007/lambda2_estime.mat').L2_eta;
data_007_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/007/Psi1_estime.mat').A_eta;
data_007_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/007/Psi2_estime.mat').B_eta;
data_007_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/007/Psi1_exact.mat').A_exact;
data_007_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/007/Psi2_exact.mat').B_exact;

data_008_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/008/eta_estime.mat').All_coefs;
data_008_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/008/lambda1_estime.mat').L1_eta;
data_008_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/008/lambda2_estime.mat').L2_eta;
data_008_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/008/Psi1_estime.mat').A_eta;
data_008_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/008/Psi2_estime.mat').B_eta;
data_008_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/008/Psi1_exact.mat').A_exact;
data_008_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/008/Psi2_exact.mat').B_exact;

data_009_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/009/eta_estime.mat').All_coefs;
data_009_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/009/lambda1_estime.mat').L1_eta;
data_009_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/009/lambda2_estime.mat').L2_eta;
data_009_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/009/Psi1_estime.mat').A_eta;
data_009_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/009/Psi2_estime.mat').B_eta;
data_009_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/009/Psi1_exact.mat').A_exact;
data_009_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/009/Psi2_exact.mat').B_exact;

data_010_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/010/eta_estime.mat').All_coefs;
data_010_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/010/lambda1_estime.mat').L1_eta;
data_010_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/010/lambda2_estime.mat').L2_eta;
data_010_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/010/Psi1_estime.mat').A_eta;
data_010_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/010/Psi2_estime.mat').B_eta;
data_010_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/010/Psi1_exact.mat').A_exact;
data_010_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/010/Psi2_exact.mat').B_exact;

data_011_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/011/eta_estime.mat').All_coefs;
data_011_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/011/lambda1_estime.mat').L1_eta;
data_011_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/011/lambda2_estime.mat').L2_eta;
data_011_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/011/Psi1_estime.mat').A_eta;
data_011_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/011/Psi2_estime.mat').B_eta;
data_011_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/011/Psi1_exact.mat').A_exact;
data_011_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/011/Psi2_exact.mat').B_exact;

data_012_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/012/eta_estime.mat').All_coefs;
data_012_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/012/lambda1_estime.mat').L1_eta;
data_012_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/012/lambda2_estime.mat').L2_eta;
data_012_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/012/Psi1_estime.mat').A_eta;
data_012_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/012/Psi2_estime.mat').B_eta;
data_012_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/012/Psi1_exact.mat').A_exact;
data_012_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/012/Psi2_exact.mat').B_exact;

data_013_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/013/eta_estime.mat').All_coefs;
data_013_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/013/lambda1_estime.mat').L1_eta;
data_013_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/013/lambda2_estime.mat').L2_eta;
data_013_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/013/Psi1_estime.mat').A_eta;
data_013_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/013/Psi2_estime.mat').B_eta;
data_013_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/013/Psi1_exact.mat').A_exact;
data_013_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/013/Psi2_exact.mat').B_exact;

data_014_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/014/eta_estime.mat').All_coefs;
data_014_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/014/lambda1_estime.mat').L1_eta;
data_014_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/014/lambda2_estime.mat').L2_eta;
data_014_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/014/Psi1_estime.mat').A_eta;
data_014_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/014/Psi2_estime.mat').B_eta;
data_014_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/014/Psi1_exact.mat').A_exact;
data_014_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/014/Psi2_exact.mat').B_exact;

data_015_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/015/eta_estime.mat').All_coefs;
data_015_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/015/lambda1_estime.mat').L1_eta;
data_015_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/015/lambda2_estime.mat').L2_eta;
data_015_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/015/Psi1_estime.mat').A_eta;
data_015_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/015/Psi2_estime.mat').B_eta;
data_015_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/015/Psi1_exact.mat').A_exact;
data_015_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/015/Psi2_exact.mat').B_exact;

data_016_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/016/eta_estime.mat').All_coefs;
data_016_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/016/lambda1_estime.mat').L1_eta;
data_016_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/016/lambda2_estime.mat').L2_eta;
data_016_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/016/Psi1_estime.mat').A_eta;
data_016_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/016/Psi2_estime.mat').B_eta;
data_016_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/016/Psi1_exact.mat').A_exact;
data_016_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/016/Psi2_exact.mat').B_exact;

data_017_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/017/eta_estime.mat').All_coefs;
data_017_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/017/lambda1_estime.mat').L1_eta;
data_017_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/017/lambda2_estime.mat').L2_eta;
data_017_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/017/Psi1_estime.mat').A_eta;
data_017_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/017/Psi2_estime.mat').B_eta;
data_017_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/017/Psi1_exact.mat').A_exact;
data_017_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/017/Psi2_exact.mat').B_exact;

data_018_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/018/eta_estime.mat').All_coefs;
data_018_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/018/lambda1_estime.mat').L1_eta;
data_018_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/018/lambda2_estime.mat').L2_eta;
data_018_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/018/Psi1_estime.mat').A_eta;
data_018_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/018/Psi2_estime.mat').B_eta;
data_018_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/018/Psi1_exact.mat').A_exact;
data_018_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/018/Psi2_exact.mat').B_exact;

data_019_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/019/eta_estime.mat').All_coefs;
data_019_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/019/lambda1_estime.mat').L1_eta;
data_019_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/019/lambda2_estime.mat').L2_eta;
data_019_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/019/Psi1_estime.mat').A_eta;
data_019_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/019/Psi2_estime.mat').B_eta;
data_019_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/019/Psi1_exact.mat').A_exact;
data_019_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/019/Psi2_exact.mat').B_exact;

data_020_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/020/eta_estime.mat').All_coefs;
data_020_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/020/lambda1_estime.mat').L1_eta;
data_020_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/020/lambda2_estime.mat').L2_eta;
data_020_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/020/Psi1_estime.mat').A_eta;
data_020_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/020/Psi2_estime.mat').B_eta;
data_020_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/020/Psi1_exact.mat').A_exact;
data_020_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/020/Psi2_exact.mat').B_exact;

data_021_eta = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/021/eta_estime.mat').All_coefs;
data_021_L1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/021/lambda1_estime.mat').L1_eta;
data_021_L2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/021/lambda2_estime.mat').L2_eta;
data_021_Psi1 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/021/Psi1_estime.mat').A_eta;
data_021_Psi2 = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/021/Psi2_estime.mat').B_eta;
data_021_Psi1_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/021/Psi1_exact.mat').A_exact;
data_021_Psi2_exact = load('/Users/Flavien/Desktop/generalized_inverse_problem/Pattern Formation/021/Psi2_exact.mat').B_exact;

%% ||A-A_ex||L2(Omega)
ecart_A_001 = sqrt(h*sum((data_001_Psi1-data_001_Psi1_exact).^2,1));
ecart_A_002 = sqrt(h*sum((data_002_Psi1-data_002_Psi1_exact).^2,1));
ecart_A_003 = sqrt(h*sum((data_003_Psi1-data_003_Psi1_exact).^2,1));
ecart_A_004 = sqrt(h*sum((data_004_Psi1-data_004_Psi1_exact).^2,1));
ecart_A_005 = sqrt(h*sum((data_005_Psi1-data_005_Psi1_exact).^2,1));
ecart_A_006 = sqrt(h*sum((data_006_Psi1-data_006_Psi1_exact).^2,1));
ecart_A_007 = sqrt(h*sum((data_007_Psi1-data_007_Psi1_exact).^2,1));
ecart_A_008 = sqrt(h*sum((data_008_Psi1-data_008_Psi1_exact).^2,1));
ecart_A_009 = sqrt(h*sum((data_009_Psi1-data_009_Psi1_exact).^2,1));
ecart_A_010 = sqrt(h*sum((data_010_Psi1-data_010_Psi1_exact).^2,1));
ecart_A_011 = sqrt(h*sum((data_011_Psi1-data_011_Psi1_exact).^2,1));
ecart_A_012 = sqrt(h*sum((data_012_Psi1-data_012_Psi1_exact).^2,1));
ecart_A_013 = sqrt(h*sum((data_013_Psi1-data_013_Psi1_exact).^2,1));
ecart_A_014 = sqrt(h*sum((data_014_Psi1-data_014_Psi1_exact).^2,1));
ecart_A_015 = sqrt(h*sum((data_015_Psi1-data_015_Psi1_exact).^2,1));
ecart_A_016 = sqrt(h*sum((data_016_Psi1-data_016_Psi1_exact).^2,1));
ecart_A_017 = sqrt(h*sum((data_017_Psi1-data_017_Psi1_exact).^2,1));
ecart_A_018 = sqrt(h*sum((data_018_Psi1-data_018_Psi1_exact).^2,1));
ecart_A_019 = sqrt(h*sum((data_019_Psi1-data_019_Psi1_exact).^2,1));
ecart_A_020 = sqrt(h*sum((data_020_Psi1-data_020_Psi1_exact).^2,1));
ecart_A_021 = sqrt(h*sum((data_021_Psi1-data_021_Psi1_exact).^2,1));

%% ||B-B_ex||L2(Omega)
ecart_B_001 = sqrt(h*sum((data_001_Psi2-data_001_Psi2_exact).^2,1));
ecart_B_002 = sqrt(h*sum((data_002_Psi2-data_002_Psi2_exact).^2,1));
ecart_B_003 = sqrt(h*sum((data_003_Psi2-data_003_Psi2_exact).^2,1));
ecart_B_004 = sqrt(h*sum((data_004_Psi2-data_004_Psi2_exact).^2,1));
ecart_B_005 = sqrt(h*sum((data_005_Psi2-data_005_Psi2_exact).^2,1));
ecart_B_006 = sqrt(h*sum((data_006_Psi2-data_006_Psi2_exact).^2,1));
ecart_B_007 = sqrt(h*sum((data_007_Psi2-data_007_Psi2_exact).^2,1));
ecart_B_008 = sqrt(h*sum((data_008_Psi2-data_008_Psi2_exact).^2,1));
ecart_B_009 = sqrt(h*sum((data_009_Psi2-data_009_Psi2_exact).^2,1));
ecart_B_010 = sqrt(h*sum((data_010_Psi2-data_010_Psi2_exact).^2,1));
ecart_B_011 = sqrt(h*sum((data_011_Psi2-data_011_Psi2_exact).^2,1));
ecart_B_012 = sqrt(h*sum((data_012_Psi2-data_012_Psi2_exact).^2,1));
ecart_B_013 = sqrt(h*sum((data_013_Psi2-data_013_Psi2_exact).^2,1));
ecart_B_014 = sqrt(h*sum((data_014_Psi2-data_014_Psi2_exact).^2,1));
ecart_B_015 = sqrt(h*sum((data_015_Psi2-data_015_Psi2_exact).^2,1));
ecart_B_016 = sqrt(h*sum((data_016_Psi2-data_016_Psi2_exact).^2,1));
ecart_B_017 = sqrt(h*sum((data_017_Psi2-data_017_Psi2_exact).^2,1));
ecart_B_018 = sqrt(h*sum((data_018_Psi2-data_018_Psi2_exact).^2,1));
ecart_B_019 = sqrt(h*sum((data_019_Psi2-data_019_Psi2_exact).^2,1));
ecart_B_020 = sqrt(h*sum((data_020_Psi2-data_020_Psi2_exact).^2,1));
ecart_B_021 = sqrt(h*sum((data_021_Psi2-data_021_Psi2_exact).^2,1));

%% ||L1||_L2(Omega)
norm_L1_001 = sqrt(h*sum((data_001_L1).^2,1));
norm_L1_002 = sqrt(h*sum((data_002_L1).^2,1));
norm_L1_003 = sqrt(h*sum((data_003_L1).^2,1));
norm_L1_004 = sqrt(h*sum((data_004_L1).^2,1));
norm_L1_005 = sqrt(h*sum((data_005_L1).^2,1));
norm_L1_006 = sqrt(h*sum((data_006_L1).^2,1));
norm_L1_007 = sqrt(h*sum((data_007_L1).^2,1));
norm_L1_008 = sqrt(h*sum((data_008_L1).^2,1));
norm_L1_009 = sqrt(h*sum((data_009_L1).^2,1));
norm_L1_010 = sqrt(h*sum((data_010_L1).^2,1));
norm_L1_011 = sqrt(h*sum((data_011_L1).^2,1));
norm_L1_012 = sqrt(h*sum((data_012_L1).^2,1));
norm_L1_013 = sqrt(h*sum((data_013_L1).^2,1));
norm_L1_014 = sqrt(h*sum((data_014_L1).^2,1));
norm_L1_015 = sqrt(h*sum((data_015_L1).^2,1));
norm_L1_016 = sqrt(h*sum((data_016_L1).^2,1));
norm_L1_017 = sqrt(h*sum((data_017_L1).^2,1));
norm_L1_018 = sqrt(h*sum((data_018_L1).^2,1));
norm_L1_019 = sqrt(h*sum((data_019_L1).^2,1));
norm_L1_020 = sqrt(h*sum((data_020_L1).^2,1));
norm_L1_021 = sqrt(h*sum((data_021_L1).^2,1));

%% ||L2||_L2(Omega)
norm_L2_001 = sqrt(h*sum((data_001_L2).^2,1));
norm_L2_002 = sqrt(h*sum((data_002_L2).^2,1));
norm_L2_003 = sqrt(h*sum((data_003_L2).^2,1));
norm_L2_004 = sqrt(h*sum((data_004_L2).^2,1));
norm_L2_005 = sqrt(h*sum((data_005_L2).^2,1));
norm_L2_006 = sqrt(h*sum((data_006_L2).^2,1));
norm_L2_007 = sqrt(h*sum((data_007_L2).^2,1));
norm_L2_008 = sqrt(h*sum((data_008_L2).^2,1));
norm_L2_009 = sqrt(h*sum((data_009_L2).^2,1));
norm_L2_010 = sqrt(h*sum((data_010_L2).^2,1));
norm_L2_011 = sqrt(h*sum((data_011_L2).^2,1));
norm_L2_012 = sqrt(h*sum((data_012_L2).^2,1));
norm_L2_013 = sqrt(h*sum((data_013_L2).^2,1));
norm_L2_014 = sqrt(h*sum((data_014_L2).^2,1));
norm_L2_015 = sqrt(h*sum((data_015_L2).^2,1));
norm_L2_016 = sqrt(h*sum((data_016_L2).^2,1));
norm_L2_017 = sqrt(h*sum((data_017_L2).^2,1));
norm_L2_018 = sqrt(h*sum((data_018_L2).^2,1));
norm_L2_019 = sqrt(h*sum((data_019_L2).^2,1));
norm_L2_020 = sqrt(h*sum((data_020_L2).^2,1));
norm_L2_021 = sqrt(h*sum((data_021_L2).^2,1));

%% ||eta-eta_ex||_infty
N_iter_001 = size(data_001_eta,2);ecart_theta_001 = norm(data_001_eta(:,N_iter_001)-Coefs_exact,Inf);
N_iter_002 = size(data_002_eta,2);ecart_theta_002 = norm(data_002_eta(:,N_iter_002)-Coefs_exact,Inf);
N_iter_003 = size(data_003_eta,2);ecart_theta_003 = norm(data_003_eta(:,N_iter_003)-Coefs_exact,Inf);
N_iter_004 = size(data_004_eta,2);ecart_theta_004 = norm(data_004_eta(:,N_iter_004)-Coefs_exact,Inf);
N_iter_005 = size(data_005_eta,2);ecart_theta_005 = norm(data_005_eta(:,N_iter_005)-Coefs_exact,Inf);
N_iter_006 = size(data_006_eta,2);ecart_theta_006 = norm(data_006_eta(:,N_iter_006)-Coefs_exact,Inf);
N_iter_007 = size(data_007_eta,2);ecart_theta_007 = norm(data_007_eta(:,N_iter_007)-Coefs_exact,Inf);
N_iter_008 = size(data_008_eta,2);ecart_theta_008 = norm(data_008_eta(:,N_iter_008)-Coefs_exact,Inf);
N_iter_009 = size(data_009_eta,2);ecart_theta_009 = norm(data_009_eta(:,N_iter_009)-Coefs_exact,Inf);
N_iter_010 = size(data_010_eta,2);ecart_theta_010 = norm(data_010_eta(:,N_iter_010)-Coefs_exact,Inf);
N_iter_011 = size(data_011_eta,2);ecart_theta_011 = norm(data_011_eta(:,N_iter_011)-Coefs_exact,Inf);
N_iter_012 = size(data_012_eta,2);ecart_theta_012 = norm(data_012_eta(:,N_iter_012)-Coefs_exact,Inf);
N_iter_013 = size(data_013_eta,2);ecart_theta_013 = norm(data_013_eta(:,N_iter_013)-Coefs_exact,Inf);
N_iter_014 = size(data_014_eta,2);ecart_theta_014 = norm(data_014_eta(:,N_iter_014)-Coefs_exact,Inf);
N_iter_015 = size(data_015_eta,2);ecart_theta_015 = norm(data_015_eta(:,N_iter_015)-Coefs_exact,Inf);
N_iter_016 = size(data_016_eta,2);ecart_theta_016 = norm(data_016_eta(:,N_iter_016)-Coefs_exact,Inf);
N_iter_017 = size(data_017_eta,2);ecart_theta_017 = norm(data_017_eta(:,N_iter_017)-Coefs_exact,Inf);
N_iter_018 = size(data_018_eta,2);ecart_theta_018 = norm(data_018_eta(:,N_iter_018)-Coefs_exact,Inf);
N_iter_019 = size(data_019_eta,2);ecart_theta_019 = norm(data_019_eta(:,N_iter_019)-Coefs_exact,Inf);
N_iter_020 = size(data_020_eta,2);ecart_theta_020 = norm(data_020_eta(:,N_iter_020)-Coefs_exact,Inf);
N_iter_021 = size(data_021_eta,2);ecart_theta_021 = norm(data_021_eta(:,N_iter_021)-Coefs_exact,Inf);

%% Paramètres des simulations
Param_001 = [1e-1,1e-4,1e-4,100,100,0,0,0,15,25,25,35];%[sigma_eps,sigma_a,sigma_q,Dm/Dt,m,q_tx,a_x,eps,lb(1),lb(2),ub(1),ub(2)]
Param_002 = [1e-2,1e-4,1e-4,100,100,0,0,0,15,25,25,35];
Param_003 = [1e-0,1e-4,1e-4,100,100,0,0,0,15,25,25,35];
Param_004 = [1e-1,1e-3,1e-4,100,100,0,0,0,15,25,25,35];
Param_005 = [1e-1,1e-5,1e-4,100,100,0,0,0,15,25,25,35];
Param_006 = [1e-1,1e-4,1e-3,100,100,0,0,0,15,25,25,35];
Param_007 = [1e-1,1e-4,1e-5,100,100,0,0,0,15,25,25,35];
Param_008 = [1e-1,1e-4,1e-4,200,100,0,0,0,15,25,25,35];
Param_009 = [1e-1,1e-4,1e-4,50,100,0,0,0,15,25,25,35];
Param_010 = [1e-1,1e-4,1e-4,100,50,0,0,0,15,25,25,35];
Param_011 = [1e-1,1e-4,1e-4,100,10,0,0,0,15,25,25,35];
Param_012 = [1e-1,1e-4,1e-4,100,1,0,0,0,15,25,25,35];
Param_013 = [1e-1,1e-4,1e-4,100,100,1e-4,0,0,15,25,25,35];
Param_014 = [1e-1,1e-4,1e-4,100,100,18.7,0,0,15,25,25,35];
Param_015 = [1e-1,1e-4,1e-4,100,100,18.55,0,0,15,25,25,35];
Param_016 = [1e-1,1e-4,1e-4,100,100,0,1e-4,0,15,25,25,35];
Param_017 = [1e-1,1e-4,1e-4,100,100,0,1e-3,0,15,25,25,35];
Param_018 = [1e-1,1e-4,1e-4,100,100,0,0,1e-1,15,25,25,35];
Param_019 = [1e-1,1e-4,1e-4,100,100,0,0,1e0,15,25,25,35];
Param_020 = [1e-1,1e-4,1e-4,100,100,0,0,0,10,20,30,40];
Param_021 = [1e-1,1e-4,1e-4,100,100,0,0,0,0,0,Inf,Inf];

%% sigma eps
figure;
X1 = categorical({'||A-A_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X1bis = categorical({'||A-A_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X2 = categorical({'||B-B_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X2bis = categorical({'||B-B_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X3 = categorical({'||\lambda_1||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X4 = categorical({'||\lambda_2||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X5 = categorical({'||\vartheta-\vartheta_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
Y = [max(ecart_A_002),max(ecart_A_001),max(ecart_A_003);
    ecart_A_002(1),ecart_A_001(1),ecart_A_003(1);
    max(ecart_B_002),max(ecart_B_001),max(ecart_B_003);
    ecart_B_002(1),ecart_B_001(1),ecart_B_003(1);
    max(norm_L1_002),max(norm_L1_001),max(norm_L1_003);
    max(norm_L2_002),max(norm_L2_001),max(norm_L2_003);
    ecart_theta_002,ecart_theta_001,ecart_theta_003];
sp1 = subplot(2,4,1);bar(X1,Y(1,:)); set(sp1,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(1,:))/2,2*max(Y(1,:))]);
sp2 = subplot(2,4,2);bar(X1bis,Y(2,:)); set(sp2,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(2,:))/2,2*max(Y(2,:))]);
sp3 = subplot(2,4,5);bar(X2,Y(3,:)); set(sp3,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(3,:))/2,2*max(Y(3,:))]);
sp4 = subplot(2,4,6);bar(X2bis,Y(4,:)); set(sp4,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(4,:))/2,2*max(Y(4,:))]);
sp5 = subplot(2,4,3);bar(X3,Y(5,:)); set(sp5,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(5,:))/2,2*max(Y(5,:))]);
sp6 = subplot(2,4,7);bar(X4,Y(6,:)); set(sp6,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(6,:))/2,2*max(Y(6,:))]);
sp7 = subplot(2,4,4);bar(X5,Y(7,:)); set(sp7,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(7,:))/2,2*max(Y(7,:))]);
st = suptitle('Influence of \sigma_{\epsilon}'); set(st,'FontSize',40);
legend('\sigma_{\epsilon} = 10^{-2}','\sigma_{\epsilon} = 10^{-1}','\sigma_{\epsilon} = 10^{0}','Location','southeastoutside','FontSize',30);

%% sigma a
figure;
X1 = categorical({'||A-A_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X1bis = categorical({'||A-A_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X2 = categorical({'||B-B_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X2bis = categorical({'||B-B_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X3 = categorical({'||\lambda_1||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X4 = categorical({'||\lambda_2||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X5 = categorical({'||\vartheta-\vartheta_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
Y = [max(ecart_A_005),max(ecart_A_001),max(ecart_A_004);
    ecart_A_005(1),ecart_A_001(1),ecart_A_004(1);
    max(ecart_B_005),max(ecart_B_001),max(ecart_B_004);
    ecart_B_005(1),ecart_B_001(1),ecart_B_004(1);
    max(norm_L1_005),max(norm_L1_001),max(norm_L1_004);
    max(norm_L2_005),max(norm_L2_001),max(norm_L2_004);
    ecart_theta_005,ecart_theta_001,ecart_theta_004];
sp1 = subplot(2,4,1);bar(X1,Y(1,:)); set(sp1,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(1,:))/2,2*max(Y(1,:))]);
sp2 = subplot(2,4,2);bar(X1bis,Y(2,:)); set(sp2,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(2,:))/2,2*max(Y(2,:))]);
sp3 = subplot(2,4,5);bar(X2,Y(3,:)); set(sp3,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(3,:))/2,2*max(Y(3,:))]);
sp4 = subplot(2,4,6);bar(X2bis,Y(4,:)); set(sp4,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(4,:))/2,2*max(Y(4,:))]);
sp5 = subplot(2,4,3);bar(X3,Y(5,:)); set(sp5,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(5,:))/2,2*max(Y(5,:))]);
sp6 = subplot(2,4,7);bar(X4,Y(6,:)); set(sp6,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(6,:))/2,2*max(Y(6,:))]);
sp7 = subplot(2,4,4);bar(X5,Y(7,:)); set(sp7,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(7,:))/2,2*max(Y(7,:))]);
st = suptitle('Influence of \sigma_{a}'); set(st,'FontSize',40);
legend('\sigma_{a} = 10^{-5}','\sigma_{a} = 10^{-4}','\sigma_{a} = 10^{-3}','southeastoutside','FontSize',30);

%% sigma q
figure;
X1 = categorical({'||A-A_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X1bis = categorical({'||A-A_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X2 = categorical({'||B-B_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X2bis = categorical({'||B-B_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X3 = categorical({'||\lambda_1||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X4 = categorical({'||\lambda_2||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X5 = categorical({'||\vartheta-\vartheta_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
Y = [max(ecart_A_007),max(ecart_A_001),max(ecart_A_006);
    ecart_A_007(1),ecart_A_001(1),ecart_A_006(1);
    max(ecart_B_007),max(ecart_B_001),max(ecart_B_006);
    ecart_B_007(1),ecart_B_001(1),ecart_B_006(1);
    max(norm_L1_007),max(norm_L1_001),max(norm_L1_006);
    max(norm_L2_007),max(norm_L2_001),max(norm_L2_006);
    ecart_theta_007,ecart_theta_001,ecart_theta_006];
sp1 = subplot(2,4,1);bar(X1,Y(1,:)); set(sp1,'YScale','linear','FontSize',20,'YGrid','on');%,'ylim',[min(Y(1,:))/2,2*max(Y(1,:))]);
sp2 = subplot(2,4,2);bar(X1bis,Y(2,:)); set(sp2,'YScale','linear','FontSize',20,'YGrid','on');%,'ylim',[min(Y(2,:))/2,2*max(Y(2,:))]);
sp3 = subplot(2,4,5);bar(X2,Y(3,:)); set(sp3,'YScale','linear','FontSize',20,'YGrid','on');%,'ylim',[min(Y(3,:))/2,2*max(Y(3,:))]);
sp4 = subplot(2,4,6);bar(X2bis,Y(4,:)); set(sp4,'YScale','linear','FontSize',20,'YGrid','on');%,'ylim',[min(Y(4,:))/2,2*max(Y(4,:))]);
sp5 = subplot(2,4,3);bar(X3,Y(5,:)); set(sp5,'YScale','linear','FontSize',20,'YGrid','on');%,'ylim',[min(Y(5,:))/2,2*max(Y(5,:))]);
sp6 = subplot(2,4,7);bar(X4,Y(6,:)); set(sp6,'YScale','linear','FontSize',20,'YGrid','on');%,'ylim',[min(Y(6,:))/2,2*max(Y(6,:))]);
sp7 = subplot(2,4,4);bar(X5,Y(7,:)); set(sp7,'YScale','linear','FontSize',20,'YGrid','on');%,'ylim',[min(Y(7,:))/2,2*max(Y(7,:))]);
st = suptitle('Influence of \sigma_{q}'); set(st,'FontSize',40);
legend('\sigma_{q} = 10^{-5}','\sigma_{q} = 10^{-4}','\sigma_{q} = 10^{-3}','southeastoutside','FontSize',30);

%% \Delta m / \Delta t
figure;
X1 = categorical({'||A-A_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X1bis = categorical({'||A-A_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X2 = categorical({'||B-B_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X2bis = categorical({'||B-B_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X3 = categorical({'||\lambda_1||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X4 = categorical({'||\lambda_2||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X5 = categorical({'||\vartheta-\vartheta_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
Y = [max(ecart_A_009),max(ecart_A_001),max(ecart_A_008);
    ecart_A_009(1),ecart_A_001(1),ecart_A_008(1);
    max(ecart_B_009),max(ecart_B_001),max(ecart_B_008);
    ecart_B_009(1),ecart_B_001(1),ecart_B_008(1);
    max(norm_L1_009),max(norm_L1_001),max(norm_L1_008);
    max(norm_L2_009),max(norm_L2_001),max(norm_L2_008);
    ecart_theta_009,ecart_theta_001,ecart_theta_008];
sp1 = subplot(2,4,1);bar(X1,Y(1,:)); set(sp1,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(1,:))/2,2*max(Y(1,:))]);
sp2 = subplot(2,4,2);bar(X1bis,Y(2,:)); set(sp2,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(2,:))/2,2*max(Y(2,:))]);
sp3 = subplot(2,4,5);bar(X2,Y(3,:)); set(sp3,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(3,:))/2,2*max(Y(3,:))]);
sp4 = subplot(2,4,6);bar(X2bis,Y(4,:)); set(sp4,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(4,:))/2,2*max(Y(4,:))]);
sp5 = subplot(2,4,3);bar(X3,Y(5,:)); set(sp5,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(5,:))/2,2*max(Y(5,:))]);
sp6 = subplot(2,4,7);bar(X4,Y(6,:)); set(sp6,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(6,:))/2,2*max(Y(6,:))]);
sp7 = subplot(2,4,4);bar(X5,Y(7,:)); set(sp7,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(7,:))/2,2*max(Y(7,:))]);
st = suptitle('Influence of the time between measurements'); set(st,'FontSize',40);
legend('\Delta m = 50\Delta t','\Delta m = 100\Delta t','\Delta m = 200\Delta t','southeastoutside','FontSize',30);

%% m
figure;
X1 = categorical({'||A-A_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X1bis = categorical({'||A-A_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X2 = categorical({'||B-B_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X2bis = categorical({'||B-B_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X3 = categorical({'||\lambda_1||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X4 = categorical({'||\lambda_2||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X5 = categorical({'||\vartheta-\vartheta_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
Y = [max(ecart_A_001),max(ecart_A_010),max(ecart_A_011),max(ecart_A_012);
    ecart_A_001(1),ecart_A_010(1),ecart_A_011(1),ecart_A_012(1);
    max(ecart_B_001),max(ecart_B_010),max(ecart_B_011),max(ecart_B_012);
    ecart_B_001(1),ecart_B_010(1),ecart_B_011(1),ecart_B_012(1);
    max(norm_L1_001),max(norm_L1_010),max(norm_L1_011),max(norm_L1_012);
    max(norm_L2_001),max(norm_L2_010),max(norm_L2_011),max(norm_L2_012);
    ecart_theta_001,ecart_theta_010,ecart_theta_011,ecart_theta_012];
sp1 = subplot(2,4,1);bar(X1,Y(1,:)); set(sp1,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(1,:))/2,2*max(Y(1,:))]);
sp2 = subplot(2,4,2);bar(X1bis,Y(2,:)); set(sp2,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(2,:))/2,2*max(Y(2,:))]);
sp3 = subplot(2,4,5);bar(X2,Y(3,:)); set(sp3,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(3,:))/2,2*max(Y(3,:))]);
sp4 = subplot(2,4,6);bar(X2bis,Y(4,:)); set(sp4,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(4,:))/2,2*max(Y(4,:))]);
sp5 = subplot(2,4,3);bar(X3,Y(5,:)); set(sp5,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(5,:))/2,2*max(Y(5,:))]);
sp6 = subplot(2,4,7);bar(X4,Y(6,:)); set(sp6,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(6,:))/2,2*max(Y(6,:))]);
sp7 = subplot(2,4,4);bar(X5,Y(7,:)); set(sp7,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(7,:))/2,2*max(Y(7,:))]);
st = suptitle('Influence of m'); set(st,'FontSize',40);
legend('m = 100','m = 50','m = 10','m = 1','southeastoutside','FontSize',30);

%% q_tx
figure;
X1 = categorical({'||A-A_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X1bis = categorical({'||A-A_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X2 = categorical({'||B-B_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X2bis = categorical({'||B-B_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X3 = categorical({'||\lambda_1||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X4 = categorical({'||\lambda_2||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X5 = categorical({'||\vartheta-\vartheta_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
Y = [max(ecart_A_001),max(ecart_A_013),max(ecart_A_015),max(ecart_A_014);
    ecart_A_001(1),ecart_A_013(1),ecart_A_015(1),ecart_A_014(1);
    max(ecart_B_001),max(ecart_B_013),max(ecart_B_015),max(ecart_B_014);
    ecart_B_001(1),ecart_B_013(1),ecart_B_015(1),ecart_B_014(1);
    max(norm_L1_001),max(norm_L1_013),max(norm_L1_015),max(norm_L1_014);
    max(norm_L2_001),max(norm_L2_013),max(norm_L2_015),max(norm_L2_014);
    ecart_theta_001,ecart_theta_013,ecart_theta_015,ecart_theta_014];
sp1 = subplot(2,4,1);bar(X1,Y(1,:)); set(sp1,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(1,:))/2,2*max(Y(1,:))]);
sp2 = subplot(2,4,2);bar(X1bis,Y(2,:)); set(sp2,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(2,:))/2,2*max(Y(2,:))]);
sp3 = subplot(2,4,5);bar(X2,Y(3,:)); set(sp3,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(3,:))/2,2*max(Y(3,:))]);
sp4 = subplot(2,4,6);bar(X2bis,Y(4,:)); set(sp4,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(4,:))/2,2*max(Y(4,:))]);
sp5 = subplot(2,4,3);bar(X3,Y(5,:)); set(sp5,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(5,:))/2,2*max(Y(5,:))]);
sp6 = subplot(2,4,7);bar(X4,Y(6,:)); set(sp6,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(6,:))/2,2*max(Y(6,:))]);
sp7 = subplot(2,4,4);bar(X5,Y(7,:)); set(sp7,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(7,:))/2,2*max(Y(7,:))]);
st = suptitle('Influence of errors in the modelling'); set(st,'FontSize',40);
legend('q_{t,x}=0','q_{t,x}=10^{-4}dW_t','k_4 = 18.55','k_4 = 18.7','southeastoutside','FontSize',30);

%% a_x
figure;
X1 = categorical({'||A-A_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X1bis = categorical({'||A-A_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X2 = categorical({'||B-B_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X2bis = categorical({'||B-B_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X3 = categorical({'||\lambda_1||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X4 = categorical({'||\lambda_2||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X5 = categorical({'||\vartheta-\vartheta_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
Y = [max(ecart_A_001),max(ecart_A_016),max(ecart_A_017);
    ecart_A_001(1),ecart_A_016(1),ecart_A_017(1);
    max(ecart_B_001),max(ecart_B_016),max(ecart_B_017);
    ecart_B_001(1),ecart_B_016(1),ecart_B_017(1);
    max(norm_L1_001),max(norm_L1_016),max(norm_L1_017);
    max(norm_L2_001),max(norm_L2_016),max(norm_L2_017);
    ecart_theta_001,ecart_theta_016,ecart_theta_017];
sp1 = subplot(2,4,1);bar(X1,Y(1,:)); set(sp1,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(1,:))/2,2*max(Y(1,:))]);
sp2 = subplot(2,4,2);bar(X1bis,Y(2,:)); set(sp2,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(2,:))/2,2*max(Y(2,:))]);
sp3 = subplot(2,4,5);bar(X2,Y(3,:)); set(sp3,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(3,:))/2,2*max(Y(3,:))]);
sp4 = subplot(2,4,6);bar(X2bis,Y(4,:)); set(sp4,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(4,:))/2,2*max(Y(4,:))]);
sp5 = subplot(2,4,3);bar(X3,Y(5,:)); set(sp5,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(5,:))/2,2*max(Y(5,:))]);
sp6 = subplot(2,4,7);bar(X4,Y(6,:)); set(sp6,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(6,:))/2,2*max(Y(6,:))]);
sp7 = subplot(2,4,4);bar(X5,Y(7,:)); set(sp7,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(7,:))/2,2*max(Y(7,:))]);
st = suptitle('Influence of errors in the initial condition'); set(st,'FontSize',40);
legend('a_x = 0','a_x = 10^{-4}N(0,1)','a_x = 10^{-3}N(0,1)','southeastoutside','FontSize',30);

%% epsilon
figure;
X1 = categorical({'||A-A_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X1bis = categorical({'||A-A_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X2 = categorical({'||B-B_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X2bis = categorical({'||B-B_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X3 = categorical({'||\lambda_1||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X4 = categorical({'||\lambda_2||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X5 = categorical({'||\vartheta-\vartheta_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
Y = [max(ecart_A_001),max(ecart_A_018),max(ecart_A_019);
    ecart_A_001(1),ecart_A_018(1),ecart_A_019(1);
    max(ecart_B_001),max(ecart_B_018),max(ecart_B_019);
    ecart_B_001(1),ecart_B_018(1),ecart_B_019(1);
    max(norm_L1_001),max(norm_L1_018),max(norm_L1_019);
    max(norm_L2_001),max(norm_L2_018),max(norm_L2_019);
    ecart_theta_001,ecart_theta_018,ecart_theta_019];
sp1 = subplot(2,4,1);bar(X1,Y(1,:)); set(sp1,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(1,:))/2,2*max(Y(1,:))]);
sp2 = subplot(2,4,2);bar(X1bis,Y(2,:)); set(sp2,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(2,:))/2,2*max(Y(2,:))]);
sp3 = subplot(2,4,5);bar(X2,Y(3,:)); set(sp3,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(3,:))/2,2*max(Y(3,:))]);
sp4 = subplot(2,4,6);bar(X2bis,Y(4,:)); set(sp4,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(4,:))/2,2*max(Y(4,:))]);
sp5 = subplot(2,4,3);bar(X3,Y(5,:)); set(sp5,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(5,:))/2,2*max(Y(5,:))]);
sp6 = subplot(2,4,7);bar(X4,Y(6,:)); set(sp6,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(6,:))/2,2*max(Y(6,:))]);
sp7 = subplot(2,4,4);bar(X5,Y(7,:)); set(sp7,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(7,:))/2,2*max(Y(7,:))]);
st = suptitle('Influence of errors in the measurement process'); set(st,'FontSize',40);
legend('\epsilon = 0','\epsilon = 10^{-1}N(0,1)','\epsilon = 10^{0}N(0,1)','southeastoutside','FontSize',30);

%% lb, ub
figure;
X1 = categorical({'||A-A_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X1bis = categorical({'||A-A_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X2 = categorical({'||B-B_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X2bis = categorical({'||B-B_{ex}(t_0,.)||_{L^{2}(\Omega)}'});
X3 = categorical({'||\lambda_1||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X4 = categorical({'||\lambda_2||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
X5 = categorical({'||\vartheta-\vartheta_{ex}||_{L^{\infty}(t_0,t_f,L^{2}(\Omega))}'});
Y = [max(ecart_A_001),max(ecart_A_020),max(ecart_A_021);
    ecart_A_001(1),ecart_A_020(1),ecart_A_021(1);
    max(ecart_B_001),max(ecart_B_020),max(ecart_B_021);
    ecart_B_001(1),ecart_B_020(1),ecart_B_021(1);
    max(norm_L1_001),max(norm_L1_020),max(norm_L1_021);
    max(norm_L2_001),max(norm_L2_020),max(norm_L2_021);
    ecart_theta_001,ecart_theta_020,ecart_theta_021];
sp1 = subplot(2,4,1);bar(X1,Y(1,:)); set(sp1,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(1,:))/2,2*max(Y(1,:))]);
sp2 = subplot(2,4,2);bar(X1bis,Y(2,:)); set(sp2,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(2,:))/2,2*max(Y(2,:))]);
sp3 = subplot(2,4,5);bar(X2,Y(3,:)); set(sp3,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(3,:))/2,2*max(Y(3,:))]);
sp4 = subplot(2,4,6);bar(X2bis,Y(4,:)); set(sp4,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(4,:))/2,2*max(Y(4,:))]);
sp5 = subplot(2,4,3);bar(X3,Y(5,:)); set(sp5,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(5,:))/2,2*max(Y(5,:))]);
sp6 = subplot(2,4,7);bar(X4,Y(6,:)); set(sp6,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(6,:))/2,2*max(Y(6,:))]);
sp7 = subplot(2,4,4);bar(X5,Y(7,:)); set(sp7,'YScale','log','FontSize',20,'YGrid','on','ylim',[min(Y(7,:))/2,2*max(Y(7,:))]);
st = suptitle('Influence of the range between l_b and u_b'); set(st,'FontSize',40);
legend('l_b = [15,25], u_b = [25,35]','l_b = [10,20], u_b = [30,40]','l_b = [0,0], u_b = [\infty,\infty]','southeastoutside','FontSize',30);



%% sigma eps
figure;
sp1 = subplot(2,3,1);plot((t0:delta_t:tf),ecart_A_002,(t0:delta_t:tf),ecart_A_001,(t0:delta_t:tf),ecart_A_003,'LineWidth',2);set(sp1,'YScale','log','FontSize',20);title('||A-A_{ex}||_{L^2(\Omega)}');
sp2 = subplot(2,3,2);plot((t0:delta_t:tf),ecart_B_002,(t0:delta_t:tf),ecart_B_001,(t0:delta_t:tf),ecart_B_003,'LineWidth',2);set(sp2,'YScale','log','FontSize',20);title('||B-B_{ex}||_{L^2(\Omega)}');
sp3 = subplot(2,3,4);plot((t0:delta_t:tf),norm_L1_002,(t0:delta_t:tf),norm_L1_001,(t0:delta_t:tf),norm_L1_003,'LineWidth',2);set(sp3,'YScale','log','FontSize',20);title('||\lambda_1||_{L^2(\Omega)}');
sp4 = subplot(2,3,5);plot((t0:delta_t:tf),norm_L2_002,(t0:delta_t:tf),norm_L2_001,(t0:delta_t:tf),norm_L2_003,'LineWidth',2);set(sp4,'YScale','log','FontSize',20);title('||\lambda_2||_{L^2(\Omega)}');
legend('\sigma_{\epsilon} = 10^{-2}','\sigma_{\epsilon} = 10^{-1}','\sigma_{\epsilon} = 10^{0}','Location','southeastoutside','FontSize',30);
sp5 = subplot(2,3,3);plot((0:N_iter_002-1),max(abs(data_002_eta-Coefs_exact),[],1),'o',(0:N_iter_001-1),max(abs(data_001_eta-Coefs_exact),[],1),'x',(0:N_iter_003-1),max(abs(data_003_eta-Coefs_exact),[],1),'+','MarkerSize',10,'LineWidth',2);set(sp5,'YScale','log','FontSize',20);title('||\vartheta-\vartheta_{ex}||_{\infty}');
st = suptitle('Influence of \sigma_{\epsilon}'); set(st,'FontSize',40);

%% sigma a
figure;
sp1 = subplot(2,3,1);plot((t0:delta_t:tf),ecart_A_005,(t0:delta_t:tf),ecart_A_001,(t0:delta_t:tf),ecart_A_004,'LineWidth',2);set(sp1,'YScale','log','FontSize',20);title('||A-A_{ex}||_{L^2(\Omega)}');
sp2 = subplot(2,3,2);plot((t0:delta_t:tf),ecart_B_005,(t0:delta_t:tf),ecart_B_001,(t0:delta_t:tf),ecart_B_004,'LineWidth',2);set(sp2,'YScale','log','FontSize',20);title('||B-B_{ex}||_{L^2(\Omega)}');
sp3 = subplot(2,3,4);plot((t0:delta_t:tf),norm_L1_005,(t0:delta_t:tf),norm_L1_001,(t0:delta_t:tf),norm_L1_004,'LineWidth',2);set(sp3,'YScale','log','FontSize',20);title('||\lambda_1||_{L^2(\Omega)}');
sp4 = subplot(2,3,5);plot((t0:delta_t:tf),norm_L2_005,(t0:delta_t:tf),norm_L2_001,(t0:delta_t:tf),norm_L2_004,'LineWidth',2);set(sp4,'YScale','log','FontSize',20);title('||\lambda_2||_{L^2(\Omega)}');
legend('\sigma_{a} = 10^{-5}','\sigma_{a} = 10^{-4}','\sigma_{a} = 10^{-3}','southeastoutside','FontSize',30);
sp5 = subplot(2,3,3);plot((0:N_iter_005-1),max(abs(data_005_eta-Coefs_exact),[],1),'o',(0:N_iter_001-1),max(abs(data_001_eta-Coefs_exact),[],1),'x',(0:N_iter_004-1),max(abs(data_004_eta-Coefs_exact),[],1),'+','MarkerSize',10,'LineWidth',2);set(sp5,'YScale','log','FontSize',20);title('||\vartheta-\vartheta_{ex}||_{\infty}');
st = suptitle('Influence of \sigma_{a}'); set(st,'FontSize',40);

%% sigma q
figure;
sp1 = subplot(2,3,1);plot((t0:delta_t:tf),ecart_A_007,(t0:delta_t:tf),ecart_A_001,(t0:delta_t:tf),ecart_A_006,'LineWidth',2);set(sp1,'YScale','log','FontSize',20);title('||A-A_{ex}||_{L^2(\Omega)}');
sp2 = subplot(2,3,2);plot((t0:delta_t:tf),ecart_B_007,(t0:delta_t:tf),ecart_B_001,(t0:delta_t:tf),ecart_B_006,'LineWidth',2);set(sp2,'YScale','log','FontSize',20);title('||B-B_{ex}||_{L^2(\Omega)}');
sp3 = subplot(2,3,4);plot((t0:delta_t:tf),norm_L1_007,(t0:delta_t:tf),norm_L1_001,(t0:delta_t:tf),norm_L1_006,'LineWidth',2);set(sp3,'YScale','log','FontSize',20);title('||\lambda_1||_{L^2(\Omega)}');
sp4 = subplot(2,3,5);plot((t0:delta_t:tf),norm_L2_007,(t0:delta_t:tf),norm_L2_001,(t0:delta_t:tf),norm_L2_006,'LineWidth',2);set(sp4,'YScale','log','FontSize',20);title('||\lambda_2||_{L^2(\Omega)}');
legend('\sigma_{q} = 10^{-5}','\sigma_{q} = 10^{-4}','\sigma_{q} = 10^{-3}','southeastoutside','FontSize',30);
sp5 = subplot(2,3,3);plot((0:N_iter_007-1),max(abs(data_007_eta-Coefs_exact),[],1),'o',(0:N_iter_001-1),max(abs(data_001_eta-Coefs_exact),[],1),'x',(0:N_iter_006-1),max(abs(data_006_eta-Coefs_exact),[],1),'+','MarkerSize',10,'LineWidth',2);set(sp5,'YScale','log','FontSize',20);title('||\vartheta-\vartheta_{ex}||_{\infty}');
st = suptitle('Influence of \sigma_{q}'); set(st,'FontSize',40);

%% \Delta m / \Delta t
figure;
sp1 = subplot(2,3,1);plot((t0:delta_t:tf),ecart_A_009,(t0:delta_t:tf),ecart_A_001,(t0:delta_t:tf),ecart_A_008,'LineWidth',2);set(sp1,'YScale','log','FontSize',20);title('||A-A_{ex}||_{L^2(\Omega)}');
sp2 = subplot(2,3,2);plot((t0:delta_t:tf),ecart_B_009,(t0:delta_t:tf),ecart_B_001,(t0:delta_t:tf),ecart_B_008,'LineWidth',2);set(sp2,'YScale','log','FontSize',20);title('||B-B_{ex}||_{L^2(\Omega)}');
sp3 = subplot(2,3,4);plot((t0:delta_t:tf),norm_L1_009,(t0:delta_t:tf),norm_L1_001,(t0:delta_t:tf),norm_L1_008,'LineWidth',2);set(sp3,'YScale','log','FontSize',20);title('||\lambda_1||_{L^2(\Omega)}');
sp4 = subplot(2,3,5);plot((t0:delta_t:tf),norm_L2_009,(t0:delta_t:tf),norm_L2_001,(t0:delta_t:tf),norm_L2_008,'LineWidth',2);set(sp4,'YScale','log','FontSize',20);title('||\lambda_2||_{L^2(\Omega)}');
legend('\Delta m = 50\Delta t','\Delta m = 100\Delta t','\Delta m = 200\Delta t','southeastoutside','FontSize',30);
sp5 = subplot(2,3,3);plot((0:N_iter_009-1),max(abs(data_009_eta-Coefs_exact),[],1),'o',(0:N_iter_001-1),max(abs(data_001_eta-Coefs_exact),[],1),'x',(0:N_iter_008-1),max(abs(data_008_eta-Coefs_exact),[],1),'+','MarkerSize',10,'LineWidth',2);set(sp5,'YScale','log','FontSize',20);title('||\vartheta-\vartheta_{ex}||_{\infty}');
st = suptitle('Influence of the time between measurements'); set(st,'FontSize',40);

%% m
figure;
sp1 = subplot(2,3,1);plot((t0:delta_t:tf),ecart_A_001,(t0:delta_t:tf),ecart_A_010,(t0:delta_t:tf),ecart_A_011,(t0:delta_t:tf),ecart_A_012,'LineWidth',2);set(sp1,'YScale','log','FontSize',20);title('||A-A_{ex}||_{L^2(\Omega)}');
sp2 = subplot(2,3,2);plot((t0:delta_t:tf),ecart_B_001,(t0:delta_t:tf),ecart_B_010,(t0:delta_t:tf),ecart_B_011,(t0:delta_t:tf),ecart_B_012,'LineWidth',2);set(sp2,'YScale','log','FontSize',20);title('||B-B_{ex}||_{L^2(\Omega)}');
sp3 = subplot(2,3,4);plot((t0:delta_t:tf),norm_L1_001,(t0:delta_t:tf),norm_L1_010,(t0:delta_t:tf),norm_L1_011,(t0:delta_t:tf),norm_L1_012,'LineWidth',2);set(sp3,'YScale','log','FontSize',20);title('||\lambda_1||_{L^2(\Omega)}');
sp4 = subplot(2,3,5);plot((t0:delta_t:tf),norm_L2_001,(t0:delta_t:tf),norm_L2_010,(t0:delta_t:tf),norm_L2_011,(t0:delta_t:tf),norm_L2_012,'LineWidth',2);set(sp4,'YScale','log','FontSize',20);title('||\lambda_2||_{L^2(\Omega)}');
legend('m = 100','m = 50','m = 10','m = 1','southeastoutside','FontSize',30);
sp5 = subplot(2,3,3);plot((0:N_iter_001-1),max(abs(data_001_eta-Coefs_exact),[],1),'o',(0:N_iter_010-1),max(abs(data_010_eta-Coefs_exact),[],1),'x',(0:N_iter_011-1),max(abs(data_011_eta-Coefs_exact),[],1),'+',(0:N_iter_012-1),max(abs(data_012_eta-Coefs_exact),[],1),'d','MarkerSize',10,'LineWidth',2);set(sp5,'YScale','log','FontSize',20);title('||\vartheta-\vartheta_{ex}||_{\infty}');
st = suptitle('Influence of m'); set(st,'FontSize',40);

%% q_tx
figure;
sp1 = subplot(2,3,1);plot((t0:delta_t:tf),ecart_A_001,(t0:delta_t:tf),ecart_A_013,(t0:delta_t:tf),ecart_A_015,(t0:delta_t:tf),ecart_A_014,'LineWidth',2);set(sp1,'YScale','log','FontSize',20);title('||A-A_{ex}||_{L^2(\Omega)}');
sp2 = subplot(2,3,2);plot((t0:delta_t:tf),ecart_B_001,(t0:delta_t:tf),ecart_B_013,(t0:delta_t:tf),ecart_B_015,(t0:delta_t:tf),ecart_B_014,'LineWidth',2);set(sp2,'YScale','log','FontSize',20);title('||B-B_{ex}||_{L^2(\Omega)}');
sp3 = subplot(2,3,4);plot((t0:delta_t:tf),norm_L1_001,(t0:delta_t:tf),norm_L1_013,(t0:delta_t:tf),norm_L1_015,(t0:delta_t:tf),norm_L1_014,'LineWidth',2);set(sp3,'YScale','log','FontSize',20);title('||\lambda_1||_{L^2(\Omega)}');
sp4 = subplot(2,3,5);plot((t0:delta_t:tf),norm_L2_001,(t0:delta_t:tf),norm_L2_013,(t0:delta_t:tf),norm_L2_015,(t0:delta_t:tf),norm_L2_014,'LineWidth',2);set(sp4,'YScale','log','FontSize',20);title('||\lambda_2||_{L^2(\Omega)}');
legend('q_{t,x}=0','q_{t,x}=10^{-4}dW_t','k_4 = 18.55','k_4 = 18.7','southeastoutside','FontSize',30);
sp5 = subplot(2,3,3);plot((0:N_iter_001-1),max(abs(data_001_eta-Coefs_exact),[],1),'o',(0:N_iter_013-1),max(abs(data_013_eta-Coefs_exact),[],1),'x',(0:N_iter_015-1),max(abs(data_015_eta-Coefs_exact),[],1),'+',(0:N_iter_014-1),max(abs(data_014_eta-Coefs_exact),[],1),'d','MarkerSize',10,'LineWidth',2);set(sp5,'YScale','log','FontSize',20);title('||\vartheta-\vartheta_{ex}||_{\infty}');
st = suptitle('Influence of errors in the modelling'); set(st,'FontSize',40);

%% a_x
figure;
sp1 = subplot(2,3,1);plot((t0:delta_t:tf),ecart_A_001,(t0:delta_t:tf),ecart_A_016,(t0:delta_t:tf),ecart_A_017,'LineWidth',2);set(sp1,'YScale','log','FontSize',20);title('||A-A_{ex}||_{L^2(\Omega)}');
sp2 = subplot(2,3,2);plot((t0:delta_t:tf),ecart_B_001,(t0:delta_t:tf),ecart_B_016,(t0:delta_t:tf),ecart_B_017,'LineWidth',2);set(sp2,'YScale','log','FontSize',20);title('||B-B_{ex}||_{L^2(\Omega)}');
sp3 = subplot(2,3,4);plot((t0:delta_t:tf),norm_L1_001,(t0:delta_t:tf),norm_L1_016,(t0:delta_t:tf),norm_L1_017,'LineWidth',2);set(sp3,'YScale','log','FontSize',20);title('||\lambda_1||_{L^2(\Omega)}');
sp4 = subplot(2,3,5);plot((t0:delta_t:tf),norm_L2_001,(t0:delta_t:tf),norm_L2_016,(t0:delta_t:tf),norm_L2_017,'LineWidth',2);set(sp4,'YScale','log','FontSize',20);title('||\lambda_2||_{L^2(\Omega)}');
legend('a_x = 0','a_x = 10^{-4}N(0,1)','a_x = 10^{-3}N(0,1)','southeastoutside','FontSize',30);
sp5 = subplot(2,3,3);plot((0:N_iter_001-1),max(abs(data_001_eta-Coefs_exact),[],1),'o',(0:N_iter_016-1),max(abs(data_016_eta-Coefs_exact),[],1),'x',(0:N_iter_017-1),max(abs(data_017_eta-Coefs_exact),[],1),'+','MarkerSize',10,'LineWidth',2);set(sp5,'YScale','log','FontSize',20);title('||\vartheta-\vartheta_{ex}||_{\infty}');
st = suptitle('Influence of errors in the initial condition'); set(st,'FontSize',40);

%% epsilon
figure;
sp1 = subplot(2,3,1);plot((t0:delta_t:tf),ecart_A_001,(t0:delta_t:tf),ecart_A_018,(t0:delta_t:tf),ecart_A_019,'LineWidth',2);set(sp1,'YScale','log','FontSize',20);title('||A-A_{ex}||_{L^2(\Omega)}');
sp2 = subplot(2,3,2);plot((t0:delta_t:tf),ecart_B_001,(t0:delta_t:tf),ecart_B_018,(t0:delta_t:tf),ecart_B_019,'LineWidth',2);set(sp2,'YScale','log','FontSize',20);title('||B-B_{ex}||_{L^2(\Omega)}');
sp3 = subplot(2,3,4);plot((t0:delta_t:tf),norm_L1_001,(t0:delta_t:tf),norm_L1_018,(t0:delta_t:tf),norm_L1_019,'LineWidth',2);set(sp3,'YScale','log','FontSize',20);title('||\lambda_1||_{L^2(\Omega)}');
sp4 = subplot(2,3,5);plot((t0:delta_t:tf),norm_L2_001,(t0:delta_t:tf),norm_L2_018,(t0:delta_t:tf),norm_L2_019,'LineWidth',2);set(sp4,'YScale','log','FontSize',20);title('||\lambda_2||_{L^2(\Omega)}');
legend('\epsilon = 0','\epsilon = 10^{-1}N(0,1)','\epsilon = 10^{0}N(0,1)','southeastoutside','FontSize',30);
sp5 = subplot(2,3,3);plot((0:N_iter_001-1),max(abs(data_001_eta-Coefs_exact),[],1),'o',(0:N_iter_018-1),max(abs(data_018_eta-Coefs_exact),[],1),'x',(0:N_iter_019-1),max(abs(data_019_eta-Coefs_exact),[],1),'+','MarkerSize',10,'LineWidth',2);set(sp5,'YScale','log','FontSize',20);title('||\vartheta-\vartheta_{ex}||_{\infty}');
st = suptitle('Influence of errors in the measurement process'); set(st,'FontSize',40);

%% lb, ub
figure;
sp1 = subplot(2,3,1);plot((t0:delta_t:tf),ecart_A_001,(t0:delta_t:tf),ecart_A_020,(t0:delta_t:tf),ecart_A_021,'LineWidth',2);set(sp1,'YScale','log','FontSize',20);title('||A-A_{ex}||_{L^2(\Omega)}');
sp2 = subplot(2,3,2);plot((t0:delta_t:tf),ecart_B_001,(t0:delta_t:tf),ecart_B_020,(t0:delta_t:tf),ecart_B_021,'LineWidth',2);set(sp2,'YScale','log','FontSize',20);title('||B-B_{ex}||_{L^2(\Omega)}');
sp3 = subplot(2,3,4);plot((t0:delta_t:tf),norm_L1_001,(t0:delta_t:tf),norm_L1_020,(t0:delta_t:tf),norm_L1_021,'LineWidth',2);set(sp3,'YScale','log','FontSize',20);title('||\lambda_1||_{L^2(\Omega)}');
sp4 = subplot(2,3,5);plot((t0:delta_t:tf),norm_L2_001,(t0:delta_t:tf),norm_L2_020,(t0:delta_t:tf),norm_L2_021,'LineWidth',2);set(sp4,'YScale','log','FontSize',20);title('||\lambda_2||_{L^2(\Omega)}');
legend('l_b = [15,25], u_b = [25,35]','l_b = [10,20], u_b = [30,40]','l_b = [0,0], u_b = [\infty,\infty]','southeastoutside','FontSize',30);
sp5 = subplot(2,3,3);plot((0:N_iter_001-1),max(abs(data_001_eta-Coefs_exact),[],1),'o',(0:N_iter_020-1),max(abs(data_020_eta-Coefs_exact),[],1),'x',(0:N_iter_021-1),max(abs(data_021_eta-Coefs_exact),[],1),'+','MarkerSize',10,'LineWidth',2);set(sp5,'YScale','log','FontSize',20);title('||\vartheta-\vartheta_{ex}||_{\infty}');
st = suptitle('Influence of the range between l_b and u_b'); set(st,'FontSize',40);
