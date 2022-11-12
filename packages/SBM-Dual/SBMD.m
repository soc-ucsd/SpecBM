function out = SBMD(A_sdp,b_sdp,c_sdp,K_sdp,opts)
    %Spectral Bundle Method - Dual Formulation
    %Author: Feng-Yi Liao
    %Update: 11/10/2022
    %Note: We use x denote y in the SBMD
    
    
    [Paras,OutOption] = Initialize(A_sdp,b_sdp,c_sdp,K_sdp,opts);
   
    %Initial point
    omega = ones(opts.m,1); 
    if opts.adaptive
        %Paras.alpha = norm(omega);
        Paras.alpha = 1;
    end
    
    %Initialize Pt and Wt
    [eig_vec,eig_val] = eig(mat(A_sdp'*omega-c_sdp));
    [~,I]             = sort(diag(eig_val),'descend');
    eig_vec           = eig_vec(:,I);
    Pt                = eig_vec(:,1:opts.MaxCols);  %transformation matrix
    Wt                = reshape(eig_vec(:,1)*eig_vec(:,1)',[],1);

    Obj               = [];
    DualFeasi         = [];
    PrimalFeasi       = [];
    Gap               = [];
    NullCount         = 0;
    
    
    
    AlgorithmTime = tic;
    
    for iter = 1:Paras.Maxiter
        %Master problem
        [Wstar,X_next,Gammastar,Sstar,Pfeasi,gap] = Direction_QP_Dual(omega,Paras,Wt,Pt);
        
        PrimalFeasi = [PrimalFeasi,Pfeasi];
        Gap         = [Gap,gap];
        
        
        %Decompose a small matrix
        [eig_vec,eig_val] = eig(Sstar);
        [~,I]             = sort(diag(eig_val),'descend');
        eig_vec           = eig_vec(:,I);
        eig_val           = eig_val(I,I);
        
        Q1                = eig_vec(:,1:Paras.EvecPast);
        Q2                = eig_vec(:,Paras.EvecPast+1:end);
        Sigma2            = eig_val(Paras.EvecPast+1:end,Paras.EvecPast+1:end);
        
        %Stopping criteria
        f1 = Paras.rho*max([eig(mat(Paras.At_sdp'*omega-Paras.c_sdp));0])-Paras.b_sdp'*omega;
        f2 = (Paras.At_sdp'*X_next-Paras.c_sdp)'*(Gammastar*Wt+vec(Pt*Sstar*Pt'))-Paras.b_sdp'*X_next;
        
        EstimatedDrop = f1 - f2;
        if EstimatedDrop < 0
            warning('something wrong');
            fprintf('f1-f2 = %.6f',EstimatedDrop);
            break;
        elseif EstimatedDrop < Paras.epislon
            Obj       = [Obj,f1];
            DualFeasi = [DualFeasi,eig_val(I(1),I(1))];
            break;
        end
        
        %Serious Step or Null Step
        f3        = Paras.rho*max([eig(mat(Paras.At_sdp'*X_next-Paras.c_sdp));0])-Paras.b_sdp'*X_next;
        CostDrop  = f1 - f3;
        %Threshold = f1-Paras.beta*(f1-f2);
        %if f3 <= Threshold
        if Paras.beta*EstimatedDrop <= CostDrop
            %serious step
            omega = X_next;
            if opts.adaptive   
                if Paras.mu*EstimatedDrop <= CostDrop
                    Paras.alpha = max(Paras.alpha/2,Paras.alphamin);
                end
                NullCount = 0;
            end
        else
            if opts.adaptive %idle
                NullCount = NullCount+1;
                if Paras.ml*EstimatedDrop >= CostDrop && NullCount >=5 && EstimatedDrop > 10^-3
                    Paras.alpha = min(Paras.alpha*2,Paras.alphamax);
                end
            end
        end
        
        %Update Wt and Pt
        Wt = (Gammastar*Wt+ reshape(Pt*Q2*Sigma2*Q2'*Pt',[],1))/(Gammastar+trace(Sigma2));
        
        %improve numerical stability
        Wt([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) = ...
            1/2*(Wt([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) + Wt([Paras.XIndOffDiagCounter,Paras.XIndOffDiag]));
        
        
        [eig_vec,eig_val] = eig(mat(A_sdp'*X_next-c_sdp));
        [~,I]             = sort(diag(eig_val),'descend');
        DualFeasi         = [DualFeasi,eig_val(I(1),I(1))];
        
        eig_vec           = eig_vec(:,I);
        Vt                = eig_vec(:,1:Paras.EvecCurrent);  
        [Pt,~]            = qr([Pt*Q1,Vt],0);
        Obj               = [Obj,f1];        
        
        if  mod(iter,OutOption.step) == 0 
            fprintf('%5d | %7.2e | %7.2e | %9.2e | %9.2e | %8.2e | %8.2e | %8.2e |\n',...
            iter,f1,EstimatedDrop,PrimalFeasi(iter),DualFeasi(iter), Gap(iter),Paras.alpha,toc(AlgorithmTime));
        end
    end
    out.Obj         = Obj;
    out.PrimalFeasi = PrimalFeasi;
    out.DualFeasi   = DualFeasi;
    out.Gap         = Gap;
end

function [Paras,OutOption] = Initialize(At_sdp,b_sdp,c_sdp,K_sdp,opts)
     % start timing
    proctime = tic;

    Paras.epislon    = opts.epislon;
    Paras.beta       = opts.beta;
    Paras.alpha      = opts.alpha;

%     IndicesAll       = BIGPSDpositionAll(opts.n,opts.dx);%The nonzero indices in a n x n matrix (including lower and upper)
%     Indices          = BIGPSDposition(opts.n,opts.dx); %The nonzero indices in a n x n matrix (only symmetric part)
%     Paras.IndicesAll = IndicesAll;
%     Paras.Indices    = Indices;
    Paras.Maxiter    = opts.Maxiter;
%     Paras.dx         = opts.dx;
    Paras.n          = opts.n;
    Paras.m          = opts.m;
    
    Paras.ml         = 0.01;  %for adapative 
    Paras.mu         = 0.7;   %for adapative 
    Paras.alphamin   = 10^-5; %for adapative 
    Paras.alphamax   = 1000;  %for adapative 
    
%     Paras.n_new = 2*Paras.dx;
%     Paras.NumOfVar_new = Paras.n_new^2;
%     Paras.NumOfVar_new_sym = Paras.n_new*(Paras.n_new+1)/2;
%     Paras.NumOfP = nchoosek(Paras.n/Paras.dx,2);%Number of Blocks
    
    [At_sdp,b_sdp,c_sdp,K_sdp,opts] = checkInputs(At_sdp,b_sdp,c_sdp,K_sdp,opts);

    Paras.b_sdp         = b_sdp;
    Paras.At_sdp        = At_sdp;
    Paras.c_sdp         = c_sdp;
    Paras.MaxCols       = opts.MaxCols;
    [xIndSym,~,~,~,~,~] = SymmetricIndices(Paras.MaxCols);
    Paras.IndicesPSD    = xIndSym;
    Paras.rho           = opts.rho;
    Paras.EvecPast      = opts.EvecPast; 
    Paras.EvecCurrent   = opts.EvecCurrent;
    
    [xIndSym,~,xIndOffDiag,~,~,xIndOffDiagCounter] = SymmetricIndices(Paras.MaxCols);
    Paras.IndicesPSD         = xIndSym;
    Paras.IndOffDiagPSD      = xIndOffDiag;
    Paras.IndOffDiagCounter  = xIndOffDiagCounter;
    
    [XIndSym,~,XIndOffDiag,~,~,XIndOffDiagCounter] = SymmetricIndices(Paras.n);
    Paras.XIndSym            = XIndSym;
    Paras.XIndOffDiag        = XIndOffDiag;
    Paras.XIndOffDiagCounter = XIndOffDiagCounter;
    
    %Some Operators in the master problem
    Paras.ATA      = Paras.At_sdp'*Paras.At_sdp;
    Paras.NumOfVar = 1+Paras.MaxCols^2;
    Paras.twoAb    = 2*Paras.At_sdp'*Paras.b_sdp;
    
    
    %Output option
    OutOption.verbose = 1;
    OutOption.K       = K_sdp;
    OutOption.step    = 10;
    OutOption.method  = 'SBMP';
    OutOption.m       = opts.m;
    OutOption.rho     = opts.rho;
    OutOption.past    = opts.EvecPast;
    OutOption.current = opts.EvecCurrent;
    PrintHeader(proctime,OutOption);
end



