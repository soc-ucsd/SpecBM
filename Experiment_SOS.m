%Authors: Feng-Yi Liao & Yang Zheng
%         SOC Lab @UC San Diego

clc;clear;
addpath('.\packages\SBM-Primal');
addpath('.\packages\SBM-Dual');
addpath('.\packages\General');
%addpath('.\packages\InOutApprox');

n = [30,35,40];
filepath = "examples\POP\";
Tran = "";
for i = 1:length(n)
    for  r = [3]
        clearvars -global -except n filepath r i Tran
        close all 
        N = n(i);
        if r == 1
            functionname = "Broyden_sphere";
            filename     = functionname+num2str(N)+"_R1"+Tran;
            file         = filepath + filename;
        elseif r ==2
            functionname = "Rosenbrock_sphere";
            filename     = functionname+num2str(N)+"_R1"+Tran;
            file         = filepath + filename;
        elseif r == 3
            functionname = "Random_sphere";
            filename     = functionname+num2str(N)+"_R1"+Tran;
            file         = filepath + filename;
        end
        load(file);        
        %%%% To recover the original cost value, substract offset1 and offset2 
        
        K_new = K;
        opts.n           = K_new.s;
        
        %opts.m           = height(A_new); 
        opts.m           = height(At); 
        opts.epislon     = 10^-4; 
        opts.beta        = 0.1; 
        opts.mu          = 0.2; 
        opts.alphamin    = 10^-5;
        opts.alphamax    = 100;
        opts.alpha       = 1; %does not matter for adaptive case 
        %opts.feasible    = false; 
        opts.adaptive    = true;
        %opts.sparse      = false;
        %rho = opts.n;
        opts.Maxiter     = 10000;
        
        opts.rho         = 10;
        switch N
            case 30
                opts.MaxCols     = 3;
                opts.EvecPast    = 0;
                opts.EvecCurrent = 3;
            case 35
                opts.MaxCols     = 5;
                opts.EvecPast    = 0;
                opts.EvecCurrent = 5;
            case 40
                opts.MaxCols     = 7;
                opts.EvecPast    = 0;
                opts.EvecCurrent = 7;
        end

        opts.solver      = "primal";
        Out  = SBM(At.',b,c,K_new,opts);

        %save("results_POPs\SBMP\"+filename+"_result.mat",'Out');
    end
end



