%Second experiment:
% Goal: Show that the dual formulation is better when X^* is 'low' rank 
% Instance: rank(X^*) = 3

% The penalty parameter is chosen by the fact that we already know the optimal solution 
% Penal term is chosen as tr(X^*)*2+2 and tr(Z^*)*2+2

%Authors: Feng-Yi Liao & Yang Zheng
%         SOC Lab @UC San Diego
 
clc;clear;
addpath('.\packages\SBM-Primal');
addpath('.\packages\SBM-Dual');
addpath('.\packages\General');

filename = "n1000m200dr997";
% Note: the problem data is large. It takes a while to load the data
load("examples\randomSDPs\"+filename+".mat");

opts.epislon        = 10^-20; 
opts.beta           = 0.4; 
opts.mu             = 0.7;
opts.alpha          = 1; %does not matter for adaptive case 
opts.adaptive       = true;


%%%%%%%%%% [Dual] %%%%%%%%%%

    opts.Maxiter     = 300;
    opts.rho         = Optimal.TrX*2+2;
    opts.MaxCols     = 2;
    opts.EvecPast    = 1;
    opts.EvecCurrent = 1;
    opts.solver      = "dual";
    Out_Dual_1_1     = SBM(At_sdp,b_sdp,c_sdp,K_sdp,opts);
 
    opts.Maxiter     = 300;
    opts.rho         = Optimal.TrX*2+2;
    opts.MaxCols     = 3;
    opts.EvecPast    = 2;
    opts.EvecCurrent = 1;
    opts.solver      = "dual";
    Out_Dual_2_1     = SBM(At_sdp,b_sdp,c_sdp,K_sdp,opts);

    opts.Maxiter     = 300;
    opts.rho         = Optimal.TrX*2+2;
    opts.MaxCols     = 4;
    opts.EvecPast    = 3;
    opts.EvecCurrent = 1;
    opts.solver      = "dual";
    Out_Dual_3_1     = SBM(At_sdp,b_sdp,c_sdp,K_sdp,opts);

    opts.Maxiter     = 300;
    opts.rho         = Optimal.TrX*2+2;
    opts.MaxCols     = 2;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 2;
    opts.solver      = "dual";
    Out_Dual_0_2     = SBM(At_sdp,b_sdp,c_sdp,K_sdp,opts);
  
    opts.Maxiter     = 300;
    opts.rho         = Optimal.TrX*2+2;
    opts.MaxCols     = 3;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 3;
    opts.solver      = "dual";
    Out_Dual_0_3     = SBM(At_sdp,b_sdp,c_sdp,K_sdp,opts);

    opts.Maxiter     = 300;
    opts.rho         = Optimal.TrX*2+2;
    opts.MaxCols     = 4;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 4;
    opts.solver      = "dual";
    Out_Dual_0_4     = SBM(At_sdp,b_sdp,c_sdp,K_sdp,opts);


%%%%%%%%%% [Primal] %%%%%%%%%%
%We do not count the first iteration for SBMP

    opts.Maxiter       = 301;
    opts.rho           = Optimal.TrZ*2+2;
    opts.MaxCols       = 2;
    opts.EvecPast      = 0;
    opts.EvecCurrent   = 2;
    opts.solver        = "primal";
    Out_Primal_0_2     = SBM(At_sdp,b_sdp,c_sdp,K_sdp,opts);
     
    opts.Maxiter       = 301;
    opts.rho           = Optimal.TrZ*2+2;
    opts.MaxCols       = 3;
    opts.EvecPast      = 0;
    opts.EvecCurrent   = 3;
    opts.solver        = "primal";
    Out_Primal_0_3     = SBM(At_sdp,b_sdp,c_sdp,K_sdp,opts);
    
    opts.Maxiter       = 301;
    opts.rho           = Optimal.TrZ*2+2;
    opts.MaxCols       = 4;
    opts.EvecPast      = 0;
    opts.EvecCurrent   = 4;
    opts.solver        = "primal";
    Out_Primal_0_4     = SBM(At_sdp,b_sdp,c_sdp,K_sdp,opts);

% save("results_rdSDPs\"+filename+"_result.mat",'Out_Dual_1_1','Out_Dual_2_1','Out_Dual_3_1','Out_Dual_0_2','Out_Dual_0_3','Out_Dual_0_4',...
%      'Out_Primal_0_2','Out_Primal_0_3','Out_Primal_0_4');