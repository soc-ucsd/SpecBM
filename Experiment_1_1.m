%First experiment: small-scale problem.
%Goal: Show that the primal formulation is better when X^* is 'high' rank 
%Instance: rank(X^*) = 97

% The penalty parameter is chosen by the fact that we already know the optimal solution 
% Penal term is chosen as tr(X^*)*2 and tr(Z^*)*2

%Authors: Feng-Yi Liao & Yang Zheng
%         SOC Lab @UC San Diego
 
clc;clear;
addpath('.\packages\SBM-Primal');
addpath('.\packages\SBM-Dual');
addpath('.\packages\General');
filename = "n100m100dr3";
load("examples\randomSDPs\"+filename+".mat");

% At_sdp              = full(At_sdp); 
% b_sdp               = full(b_sdp); 
% c_sdp               = full(c_sdp);
opts.n              = K_sdp.s; 
opts.m              = height(At_sdp); 
opts.epislon        = 10^-20; 

opts.beta           = 0.25; 
opts.mu             = 0.5; 
opts.alpha          = 50; %does not matter for adaptive case 
opts.feasible       = false; 
opts.adaptive       = true;
opts.sparse         = false;
opts.DynamicRho     = false;
opts.DynamicMaxCols = false;

%%%%%%%%%% [Primal] %%%%%%%%%%
%We do not count the first iteration for SBMP

    opts.Maxiter     = 201;
    opts.rho         = Optimal.TrZ*2+2;
    opts.MaxCols     = 2;
    opts.EvecPast    = 1;
    opts.EvecCurrent = 1;
    Out_Primal_1_1   = Copy_of_SBMP(At_sdp,b_sdp,c_sdp,K_sdp,opts);
 
    opts.Maxiter     = 201;
    opts.rho         = Optimal.TrZ*2+2;
    opts.MaxCols     = 3;
    opts.EvecPast    = 2;
    opts.EvecCurrent = 1;
    Out_Primal_2_1   = Copy_of_SBMP(At_sdp,b_sdp,c_sdp,K_sdp,opts);

    opts.Maxiter     = 201;
    opts.rho         = Optimal.TrZ*2+2;
    opts.MaxCols     = 4;
    opts.EvecPast    = 3;
    opts.EvecCurrent = 1;
    Out_Primal_3_1   = Copy_of_SBMP(At_sdp,b_sdp,c_sdp,K_sdp,opts);

    opts.Maxiter     = 201;
    opts.rho         = Optimal.TrZ*2+2;
    opts.MaxCols     = 2;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 2;
    Out_Primal_0_2   = Copy_of_SBMP(At_sdp,b_sdp,c_sdp,K_sdp,opts);
  
    opts.Maxiter     = 201;
    opts.rho         = Optimal.TrZ*2+2;
    opts.MaxCols     = 3;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 3;
    Out_Primal_0_3   = Copy_of_SBMP(At_sdp,b_sdp,c_sdp,K_sdp,opts);

    opts.Maxiter     = 201;
    opts.rho         = Optimal.TrZ*2+2;
    opts.MaxCols     = 4;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 4;
    Out_Primal_0_4   = Copy_of_SBMP(At_sdp,b_sdp,c_sdp,K_sdp,opts);
% 
% 
% %%%%%%%%%% [Dual] %%%%%%%%%%
%     
% 
    opts.Maxiter     = 200;
    opts.rho         = Optimal.TrX*2+2;
    opts.MaxCols     = 2;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 2;
    Out_Dual_0_2     = SBMD(At_sdp,b_sdp,c_sdp,K_sdp,opts);
%      
    opts.Maxiter     = 200;
    opts.rho         = Optimal.TrX*2+2;
    opts.MaxCols     = 3;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 3;
    Out_Dual_0_3     = SBMD(At_sdp,b_sdp,c_sdp,K_sdp,opts);
    
    opts.Maxiter     = 200;
    opts.rho         = Optimal.TrX*2+2;
    opts.MaxCols     = 4;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 4;
    Out_Dual_0_4     = SBMD(At_sdp,b_sdp,c_sdp,K_sdp,opts);
% 
% save("results_rdSDPs_new\"+filename+"_result.mat",'Out_Primal_1_1','Out_Primal_2_1','Out_Primal_3_1','Out_Primal_0_2','Out_Primal_0_3','Out_Primal_0_4',...
%      'Out_Dual_0_2','Out_Dual_0_3','Out_Dual_0_4');