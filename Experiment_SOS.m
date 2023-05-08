%Authors: Feng-Yi Liao & Yang Zheng
%         SOC Lab @UC San Diego

clc;clear;
addpath('.\packages\SBM-Primal');
addpath('.\packages\SBM-Dual');
addpath('.\packages\General');
%addpath('.\packages\InOutApprox');

n = [30];
filepath = "examples\POP\";
%Tran     = "_Tran.mat";
Tran = "";
for i = 1:length(n)
    for  r = [1]
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
            functionname = "Quartic_sphere";
            filename     = functionname+num2str(N)+"_R1"+Tran;
            file         = filepath + filename;
        elseif r == 4
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
        opts.feasible    = false; 
        opts.adaptive    = true;
        opts.sparse      = false;
        %rho = opts.n;
        opts.Maxiter     = 10000;
        
        opts.rho         = 10;
        switch N
            case 30
                opts.MaxCols     = 3;%floor(opts.n/10);
                opts.EvecPast    = 0;
                opts.EvecCurrent = 3;%floor(opts.n/10);
            case 35
                opts.MaxCols     = 5;%floor(opts.n/10);
                opts.EvecPast    = 0;
                opts.EvecCurrent = 5;%floor(opts.n/10);
            case 40
                opts.MaxCols     = 7;%floor(opts.n/10);
                opts.EvecPast    = 0;
                opts.EvecCurrent = 7;%floor(opts.n/10);
            case 45
                opts.MaxCols     = 8;%floor(opts.n/10);
                opts.EvecPast    = 0;
                opts.EvecCurrent = 8;%floor(opts.n/10);
            case 50
                opts.MaxCols     = 5;%floor(opts.n/10);
                opts.EvecPast    = 0;
                opts.EvecCurrent = 5;%floor(opts.n/10); 
            otherwise
                opts.MaxCols     = 3;%floor(opts.n/10);
                opts.EvecPast    = 0;
                opts.EvecCurrent = 3;%floor(opts.n/10);
        end



%         if r == 1 
%             opts.rho         = 10;
%             opts.MaxCols     = 5;%floor(opts.n/10);
%             opts.EvecPast    = 0;
%             opts.EvecCurrent = 5;%floor(opts.n/10);
%         elseif r == 2 
%             opts.rho         = 10;
%             opts.MaxCols     = 5;%floor(opts.n/10);
%             opts.EvecPast    = 0;
%             opts.EvecCurrent = 5;%floor(opts.n/10);
%         elseif r == 3
%             opts.rho         = 10;
%             opts.MaxCols     = 5;%floor(opts.n/10);
%             opts.EvecPast    = 0;
%             opts.EvecCurrent = 5;%floor(opts.n/10);
%         elseif r == 4
%             opts.rho         = 4;
%             opts.MaxCols     = 5;%floor(opts.n/10);
%             opts.EvecPast    = 0;
%             opts.EvecCurrent = 5;%floor(opts.n/10);
%         end
        
        if N<40
            opts.epislonphase1  = 10^-3;
            opts.DynamicRho     = false;
            opts.DynamicMaxCols = false; 
            opts.MaxCols2       = 15;%floor(opts.n/10);
            opts.EvecPast2      = 0;
            opts.EvecCurrent2   = 15;%floor(opts.n/10);
        elseif N == 40
            opts.epislonphase1  = 10^-3;
            opts.DynamicRho     = false;
            opts.DynamicMaxCols = false; 
            opts.MaxCols2       = 15;%floor(opts.n/10);
            opts.EvecPast2      = 0;
            opts.EvecCurrent2   = 15;%floor(opts.n/10);
        elseif N == 45
            opts.epislonphase1  = 10^-3;
            opts.DynamicRho     = false;
            opts.DynamicMaxCols = false; 
            opts.MaxCols2       = 15;%floor(opts.n/10);
            opts.EvecPast2      = 0;
            opts.EvecCurrent2   = 15;%floor(opts.n/10);
        elseif N == 50
            opts.epislonphase1  = 10^-3;
            opts.DynamicRho     = false;
            opts.DynamicMaxCols = false; 
            opts.MaxCols2       = 10;%floor(opts.n/10);
            opts.EvecPast2      = 0;
            opts.EvecCurrent2   = 10;%floor(opts.n/10);  
        end
        opts.rho1 = opts.rho; %does not matter 
        opts.rho2 = opts.rho1; %does not matter
%         %Dynamic rho 
%         
%         if N >=30
%             opts.rho1        = 100;
%         else
%             opts.rho1        = 100;
%         end        
%         opts.rho2        = 4;
        
%         fileINVAAT       = filepath + "INV\" + functionname+num2str(N)+"_R1_"+"INVAAT.mat";
%         load(fileINVAAT);
        %opts.AAT_INV  = INVAAT;
        %c_new_nor     = c_new/norm(c_new);

        %Out  = Copy_of_SBMP(A_new,b_new,c_new,K_new,opts);
        Out  = Copy_of_SBMP(At.',b,c,K_new,opts);

        %save("results_POPs\SBMP\"+filename+"_result.mat",'Out');
    end
end



