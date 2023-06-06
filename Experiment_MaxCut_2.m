%Matrix Cut
%Authors: Feng-Yi Liao & Yang Zheng
%         SOC Lab @UC San Diego
 
clc;clear;
addpath('.\packages\SBM-Primal');
addpath('.\packages\SBM-Dual');
addpath('.\packages\General');

filename = "G25";
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
    opts.MaxCols     = 18;
    opts.EvecPast    = 17;
    opts.EvecCurrent = 1;
    opts.solver      = "dual";
    Out_Dual_17_1    = SBM(At,b,c,K,opts);
 
    opts.Maxiter     = 300;
    opts.rho         = height(At)*2+2;
    opts.MaxCols     = 19;
    opts.EvecPast    = 18;
    opts.EvecCurrent = 1;
    opts.solver      = "dual";
    Out_Dual_18_1    = SBM(At,b,c,K,opts);

    opts.Maxiter     = 300;
    opts.rho         = height(At)*2+2;
    opts.MaxCols     = 20;
    opts.EvecPast    = 19;
    opts.EvecCurrent = 1;
    opts.solver      = "dual";
    Out_Dual_19_1    = SBM(At,b,c,K,opts);

    opts.Maxiter     = 300;
    opts.rho         = height(At)*2+2;
    opts.MaxCols     = 18;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 18;
    opts.solver      = "dual";
    Out_Dual_0_18    = SBM(At,b,c,K,opts);
  
    opts.Maxiter     = 300;
    opts.rho         = height(At)*2+2;
    opts.MaxCols     = 19;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 19;
    opts.solver      = "dual";
    Out_Dual_0_19    = SBM(At,b,c,K,opts);

    opts.Maxiter     = 300;
    opts.rho         = height(At)*2+2;
    opts.MaxCols     = 20;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 20;
    opts.solver      = "dual";
    Out_Dual_0_20    = SBM(At,b,c,K,opts);


%%%%%%%%%% [Primal] %%%%%%%%%%
%We do not count the first iteration for SBMP

    rho_p = trace(mat(c))-height(At)*min(eig(mat(c)));

    opts.Maxiter       = 301;
    opts.rho           = 2*rho_p+2;
    opts.MaxCols       = 18;
    opts.EvecPast      = 0;
    opts.EvecCurrent   = 18;
    opts.solver        = "primal";
    Out_Primal_0_18    = SBM(At,b,c,K,opts);
     
    opts.Maxiter       = 301;
    opts.rho           = 2*rho_p+2;
    opts.MaxCols       = 19;
    opts.EvecPast      = 0;
    opts.EvecCurrent   = 19;
    opts.solver        = "primal";
    Out_Primal_0_19    = SBM(At,b,c,K,opts);
    
    opts.Maxiter       = 301;
    opts.rho           = 2*rho_p+2;
    opts.MaxCols       = 20;
    opts.EvecPast      = 0;
    opts.EvecCurrent   = 20;
    opts.solver        = "primal";
    Out_Primal_0_20    = SBM(At,b,c,K,opts);

% save("results_MaxCut\"+filename+"_result.mat",'Out_Dual_17_1','Out_Dual_18_1','Out_Dual_19_1','Out_Dual_0_18','Out_Dual_0_19','Out_Dual_0_20',...
%      'Out_Primal_0_18','Out_Primal_0_19','Out_Primal_0_20');