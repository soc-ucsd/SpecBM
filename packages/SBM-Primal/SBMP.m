function Out = SBMP(A_sdp,b_sdp,c_sdp,K_sdp,opts)
    % Spectral Bundle Method - Primal Formulation
    % Authors: Feng-Yi Liao & Yang Zheng
    %          SOC Lab @UC San Diego
    % Update : 11/10/2022
    % Note   : We denote Xnext as X_{t+1} in the paper 
    
    
    [Paras,OutOption] = Initialize(A_sdp,b_sdp,c_sdp,K_sdp,opts);
   
    %Initial point
    if opts.feasible %idle
        [~,Omegat] = InnerApproximation(A_sdp,b_sdp,c_sdp,K_sdp,1,1);
        omegat     = reshape(Omegat,[],1);
    else
        Omegat     = eye(Paras.n);
        omegat     = reshape(Omegat,[],1);
    end
    
    if opts.adaptive   
        %Paras.alpha = norm(Omegat,'Fro');
        Paras.alpha = 1;
    end
    
    %Initialize Pt and Wt
    [eig_vec,eig_val] = eig(mat(-omegat));
    [~,I]             = sort(diag(eig_val),'descend');
    eig_vec           = eig_vec(:,I);
    Pt                = eig_vec(:,1:opts.MaxCols);  %transformation matrix
    Wt                = reshape(eig_vec(:,1)*eig_vec(:,1)',[],1);
    
    Obj               = zeros(opts.Maxiter,1); 
    DualFeasi         = zeros(opts.Maxiter,1); 
    PrimalFeasi       = zeros(opts.Maxiter,1); 
    Gap               = zeros(opts.Maxiter,1);
    NullCount         = 0;
    
    AlgorithmTime = tic;
    
    for iter = 1:Paras.Maxiter
        %Master problem
        if iter >1 
           [Wstar,X_next,Gammastar,Sstar,Dfeasi,gap] = Direction_QP_Primal(omegat,Paras,Wt,Pt,true);
        else
           [Wstar,X_next,Gammastar,Sstar,Dfeasi,gap] = Direction_QP_Primal(omegat,Paras,Wt,Pt,false);
        end
        
        %We do not record the first iteration if the initial point is not feasible
        if opts.feasible || iter>1
            DualFeasi(iter) = Dfeasi;
            Gap(iter) = gap;
        end
        
        
        %Decompose a small matrix
        [eig_vec,eig_val] = eig(full(Sstar));
        [~,I]             = sort(diag(eig_val),'descend');
        eig_vec           = eig_vec(:,I);
        eig_val           = eig_val(I,I);
        
        Q1                = eig_vec(:,1:Paras.EvecPast);
        Q2                = eig_vec(:,Paras.EvecPast+1:end);
        Sigma2            = eig_val(Paras.EvecPast+1:end,Paras.EvecPast+1:end);
        
        Omegat = reshape(omegat,Paras.n,Paras.n);
        if issymmetric(Omegat)
            f1 = Paras.c_sdp'*omegat  + Paras.rho*max([eig(-Omegat);0]);
        else
            warning('X Not Symmetric');
        end
        
        f2 = (Paras.c_sdp-Wstar)'*X_next;
        
        EstimatedDrop = f1 - f2;
        
        if EstimatedDrop <0 && iter > 1 
            warning('something wrong');
            fprintf('f1-f2 = %.6f',EstimatedDrop);
            break;
        elseif EstimatedDrop < Paras.epislon && iter >1
            Obj(iter)         = f1;
            PrimalFeasi(iter) = eig_val(I(1),I(1));
            break;
        end
        
        
        %Serious Step or Null Step        
        Y = reshape(-X_next,Paras.n,Paras.n);
        if issymmetric(Y)
            f3 = Paras.c_sdp'*X_next + Paras.rho*max([eig(Y);0]);
        else
            warning('Y Not Symmetric');
        end
        
        CostDrop = f1 - f3;
        %Threshold = f1-Paras.beta*(f1-f2); 
        if opts.feasible || iter>1
            %if f3 <= Threshold
            if Paras.beta*EstimatedDrop <= CostDrop
                %serious step
                omegat = X_next;
                if opts.adaptive %idle
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
        else %only for the first iteration
            omegat = X_next;
        end
        
        %Update Wt and Pt
        Wt = (Gammastar*Wt+ reshape(Pt*Q2*Sigma2*Q2'*Pt',[],1))/(Gammastar+trace(Sigma2));
        
        %improve numerical stability
        Wt([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) = ...
            1/2*(Wt([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) + Wt([Paras.XIndOffDiagCounter,Paras.XIndOffDiag]));
        
        [eig_vec,eig_val] = eig(reshape(-X_next,Paras.n,[]));
        [~,I]             = sort(diag(eig_val),'descend');
        if opts.feasible || iter>1
            PrimalFeasi(iter) = -eig_val(I(1),I(1));
        end
        
        eig_vec = eig_vec(:,I);
        Vt      = eig_vec(:,1:Paras.EvecCurrent);   
        [Pt,~]  = qr([Pt*Q1,Vt],0);
        if opts.feasible || iter>1
            Obj(iter) = f1;
        end
        
        if iter > 1 && mod(iter,OutOption.step) == 0                 
           fprintf('%5d  %7.2e  %7.2e  %9.2e  %9.2e  %8.2e  %8.2e  %8.2e \n',...
           iter,f1,EstimatedDrop,PrimalFeasi(iter-1),DualFeasi(iter-1), Gap(iter-1),Paras.alpha,toc(AlgorithmTime));
        end
    end
    Out.Obj         = Obj;
    Out.PrimalFeasi = PrimalFeasi;
    Out.DualFeasi   = DualFeasi;
    Out.Gap         = Gap;
end

function [Paras,OutOption] = Initialize(At_sdp,b_sdp,c_sdp,K_sdp,opts)
    % start timing
    proctime = tic;
    
    %Parameters intilization
    Paras.epislon     = opts.epislon;
    Paras.beta        = opts.beta;
    Paras.alpha       = opts.alpha;
%     IndicesAll        = BIGPSDpositionAll(opts.n,opts.dx);%The nonzero indices in a n x n matrix (including lower and upper)
%     Indices           = BIGPSDposition(opts.n,opts.dx); %The nonzero indices in a n x n matrix (only symmetric part)
%     Paras.IndicesAll  = IndicesAll;
%     Paras.Indices     = Indices;
    Paras.Maxiter     = opts.Maxiter;
%     Paras.dx               = opts.dx;
    Paras.n           = opts.n;
    Paras.m           = opts.m;
    
    Paras.ml          = 0.01;
    Paras.mu          = 0.7;
    Paras.alphamin    = 10^-5; %for adapative 
    Paras.alphamax    = 1000;  %for adapative 
    
%     Paras.n_new            = 2*Paras.dx;
%     Paras.NumOfVar_new     = Paras.n_new^2;
%     Paras.NumOfVar_new_sym = Paras.n_new*(Paras.n_new+1)/2;
%     Paras.NumOfP           = nchoosek(Paras.n/Paras.dx,2);%Number of Blocks
    
    [At_sdp,b_sdp,c_sdp,K_sdp,opts] = checkInputs(At_sdp,b_sdp,c_sdp,K_sdp,opts);

    Paras.b_sdp            = b_sdp;
    Paras.At_sdp           = At_sdp;
    Paras.c_sdp            = c_sdp;
    Paras.K_sdp            = K_sdp;
    Paras.MaxCols          = opts.MaxCols;
    Paras.NumOfVar         = 1+Paras.MaxCols^2;
    
    [xIndSym,~,xIndOffDiag,~,~,xIndOffDiagCounter] = SymmetricIndices(Paras.MaxCols);
    Paras.IndicesPSD        = xIndSym;
    Paras.IndOffDiagPSD     = xIndOffDiag;
    Paras.IndOffDiagCounter = xIndOffDiagCounter;
    
    
    Paras.rho         = opts.rho;
    Paras.EvecPast    = opts.EvecPast;
    Paras.EvecCurrent = opts.EvecCurrent;
    
    [XIndSym,~,XIndOffDiag,~,~,XIndOffDiagCounter] = SymmetricIndices(Paras.n);
    Paras.XIndSym            = XIndSym;
    Paras.XIndOffDiag        = XIndOffDiag;
    Paras.XIndOffDiagCounter = XIndOffDiagCounter;
    
    [mIndSym,~,mIndOffDiag,~,~,mIndOffDiagCounter] = SymmetricIndices(Paras.m);
    Paras.mIndSym            = mIndSym;
    Paras.mIndOffDiag        = mIndOffDiag;
    Paras.mIndOffDiagCounter = mIndOffDiagCounter;
    
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


