%Matrix Completion
%Authors: Feng-Yi Liao & Yang Zheng
%         SOC Lab @UC San Diego
 
clc;clear;
addpath('.\packages\SBM-Primal');
addpath('.\packages\SBM-Dual');
addpath('.\packages\General');

filename = "G1";
load("examples\MaxCut\"+filename+".mat");

opts.n              = K.s; 
opts.m              = height(At); 
opts.epislon        = 10^-20; 
opts.beta           = 0.25;
opts.mu             = 0.5;
opts.alphamin       = 10^-5;
opts.alphamax       = 100;
opts.alpha          = 1; %does not matter for adaptive case 
opts.adaptive       = true;

%%%%%%%%% Mosek %%%%%%%%%%%%
%     [res,time] = SolveMosek(At,b,c,K);
%     save("results_MaxCut\"+filename+"_Mosek",'res','time');
%     X          = ones(opts.m);
%     Xtril      = tril(X,0);
%     v          = find(Xtril);
%     Y          = zeros(opts.n);
%     Y(v)       = Y(v) = res.sol.itr.barx;
    
   



%%%%%%%%%% [Dual] %%%%%%%%%%

    opts.Maxiter     = 300;
    opts.rho         = height(At)*2+2;
    opts.MaxCols     = 12;
    opts.EvecPast    = 11;
    opts.EvecCurrent = 1;
    opts.solver      = "dual";
    Out_Dual_11_1    = SBM(At,b,c,K,opts);
 
    opts.Maxiter     = 300;
    opts.rho         = height(At)*2+2;
    opts.MaxCols     = 13;
    opts.EvecPast    = 12;
    opts.EvecCurrent = 1;
    opts.solver      = "dual";
    Out_Dual_12_1    = SBM(At,b,c,K,opts);

    opts.Maxiter     = 300;
    opts.rho         = height(At)*2+2;
    opts.MaxCols     = 14;
    opts.EvecPast    = 13;
    opts.EvecCurrent = 1;
    opts.solver      = "dual";
    Out_Dual_13_1    = SBM(At,b,c,K,opts);

    opts.Maxiter     = 300;
    opts.rho         = height(At)*2+2;
    opts.MaxCols     = 12;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 12;
    opts.solver      = "dual";
    Out_Dual_0_12    = SBM(At,b,c,K,opts);
  
    opts.Maxiter     = 300;
    opts.rho         = height(At)*2+2;
    opts.MaxCols     = 13;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 13;
    opts.solver      = "dual";
    Out_Dual_0_13    = SBM(At,b,c,K,opts);

    opts.Maxiter     = 300;
    opts.rho         = height(At)*2+2;
    opts.MaxCols     = 14;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 14;
    opts.solver      = "dual";
    Out_Dual_0_14    = SBM(At,b,c,K,opts);


%%%%%%%%%% [Primal] %%%%%%%%%%
%We do not count the first iteration for SBMP

    rho_p = trace(mat(c))-height(At)*min(eig(mat(c)));

    opts.Maxiter       = 301;
    opts.rho           = 2*rho_p+2;
    opts.MaxCols       = 12;
    opts.EvecPast      = 0;
    opts.EvecCurrent   = 12;
    opts.solver        = "primal";
    Out_Primal_0_12    = SBM(At,b,c,K,opts);
     
    opts.Maxiter       = 301;
    opts.rho           = 2*rho_p+2;
    opts.MaxCols       = 13;
    opts.EvecPast      = 0;
    opts.EvecCurrent   = 13;
    opts.solver        = "primal";
    Out_Primal_0_13    = SBM(At,b,c,K,opts);
    
    opts.Maxiter       = 301;
    opts.rho           = 2*rho_p+2;
    opts.MaxCols       = 14;
    opts.EvecPast      = 0;
    opts.EvecCurrent   = 14;
    opts.solver        = "primal";
    Out_Primal_0_14    = SBM(At,b,c,K,opts);

% save("results_MaxCut\"+filename+"_result.mat",'Out_Dual_11_1','Out_Dual_12_1','Out_Dual_13_1','Out_Dual_0_12','Out_Dual_0_13','Out_Dual_0_14',...
%      'Out_Primal_0_12','Out_Primal_0_13','Out_Primal_0_14');