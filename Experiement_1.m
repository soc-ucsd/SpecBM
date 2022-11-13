%First experiment: small-scale problem.
%Goal: Show that the primal formulation is better when X^* is 'high' rank 

% The penalty parameter is chosen by the fact that we already know the optimal solution 
% Penal term is chosen as tr(X^*)*2+1 and tr(Z^*)*2+1

%Authors: Feng-Yi Liao & Yang Zheng
%         SOC Lab @UC San Diego
 
clc;clear;
addpath('.\packages\SBM-Primal');
addpath('.\packages\SBM-Dual');
addpath('.\packages\General');
load('examples\n100m100dr3.mat');

At_sdp        = full(At_sdp); 
b_sdp         = full(b_sdp); 
c_sdp         = full(c_sdp);
opts.n        = K_sdp.s; 
opts.m        = height(At_sdp); 
opts.epislon  = 10^-20; 
opts.beta     = 0.2; 
opts.alpha    = 50; %does not matter for adaptive case 
opts.feasible = false; 
opts.adaptive = true;

%%%%%%%%%% [Primal] %%%%%%%%%%
%We do not count the first iteration for SBMP

    opts.Maxiter     = 201;
    opts.rho         = Optimal.TrZ*2+1;
    opts.MaxCols     = 3;
    opts.EvecPast    = 2;
    opts.EvecCurrent = 1;
    Out_Primal_2_1   = SBMP(At_sdp,b_sdp,c_sdp,K_sdp,opts);
 
    opts.Maxiter     = 201;
    opts.rho         = Optimal.TrZ*2+1;
    opts.MaxCols     = 5;
    opts.EvecPast    = 4;
    opts.EvecCurrent = 1;
    Out_Primal_4_1   = SBMP(At_sdp,b_sdp,c_sdp,K_sdp,opts);

    opts.Maxiter     = 201;
    opts.rho         = Optimal.TrZ*2+1;
    opts.MaxCols     = 10;
    opts.EvecPast    = 9;
    opts.EvecCurrent = 1;
    Out_Primal_9_1   = SBMP(At_sdp,b_sdp,c_sdp,K_sdp,opts);

    opts.Maxiter     = 201;
    opts.rho         = Optimal.TrZ*2+1;
    opts.MaxCols     = 3;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 3;
    Out_Primal_0_3   = SBMP(At_sdp,b_sdp,c_sdp,K_sdp,opts);
  
    opts.Maxiter     = 201;
    opts.rho         = Optimal.TrZ*2+1;
    opts.MaxCols     = 5;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 5;
    Out_Primal_0_5   = SBMP(At_sdp,b_sdp,c_sdp,K_sdp,opts);

    opts.Maxiter     = 201;
    opts.rho         = Optimal.TrZ*2+1;
    opts.MaxCols     = 10;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 10;
    Out_Primal_0_10  = SBMP(At_sdp,b_sdp,c_sdp,K_sdp,opts);


%%%%%%%%%% [Dual] %%%%%%%%%%
    

    opts.Maxiter     = 200;
    opts.rho         = Optimal.TrX*2+1;
    opts.MaxCols     = 3;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 3;
    Out_Dual_0_3     = SBMD(At_sdp,b_sdp,c_sdp,K_sdp,opts);
     
    opts.Maxiter     = 200;
    opts.rho         = Optimal.TrX*2+1;
    opts.MaxCols     = 5;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 5;
    Out_Dual_0_5     = SBMD(At_sdp,b_sdp,c_sdp,K_sdp,opts);
    
    opts.Maxiter     = 200;
    opts.rho         = Optimal.TrX*2+1;
    opts.MaxCols     = 10;
    opts.EvecPast    = 0;
    opts.EvecCurrent = 10;
    Out_Dual_0_10    = SBMD(At_sdp,b_sdp,c_sdp,K_sdp,opts);