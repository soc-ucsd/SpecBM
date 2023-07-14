function Out = SBMD(At,b,c,K,opts)
    %Spectral Bundle Method - Dual Formulation
    %Authors: Feng-Yi Liao & Yang Zheng
    %         SOC Lab @UC San Diego
    %Note   : We use x denote y in the SBMD
    
    fprintf('Spectral Bundle Method Dual Starts\n');
    [Paras,OutOption] = Initialize(At,b,c,K,opts);
   
    AlgorithmTime = tic;

    %Initial point
    omegat = zeros(opts.m,1); 
    if Paras.adaptive
        Paras.alpha = 0.5;
    end
    
    %Initialize Pt and Wt
    [eig_vec,eig_val] = eig(mat(At'*omegat-c));
    [~,I]             = sort(diag(eig_val),'descend');
    eig_vec           = eig_vec(:,I);
    Pt                = eig_vec(:,1:opts.MaxCols);  %transformation matrix
    Wt                = reshape(eig_vec(:,1)*eig_vec(:,1)',[],1);

    Obj               = zeros(opts.Maxiter,1); 
    DualSemiFeasi     = zeros(opts.Maxiter,1); 
    PrimalFeasi       = zeros(opts.Maxiter,1); 
    Gap               = zeros(opts.Maxiter,1);
    RelativeAccry     = zeros(opts.Maxiter,1);
    NullCount         = 0;
    
    normb             = norm(Paras.b);
    normC             = norm(Paras.c);

    Out.RelativeDFeasi = [];
    Out.RelativePFeasi = [];
    Out.RelativeGap    = [];

    Out.DescentRelativeDFeasi = [];
    Out.DescentRelativePFeasi = [];
    Out.DescentRelativeGap    = [];
    Out.DescentDualSemiFeasi  = [];
    Out.DescentCost           = [];
    
    DescentCount              = 1;
    DescentFlag               = true;

    
    
    for iter = 1:Paras.Maxiter
        if iter == 1
            [Wstar,X_next,Gammastar,Sstar,Pfeasi,gap,Old_G] = Direction_QP_Dual(omegat,Paras,Wt,Pt);
        else
            if DescentFlag
                [Wstar,X_next,Gammastar,Sstar,Pfeasi,gap,Old_G] = Direction_QP_Dual(omegat,Paras,Wt,Pt);
            else
                [Wstar,X_next,Gammastar,Sstar,Pfeasi,gap,Old_G] = Direction_QP_Dual(omegat,Paras,Wt,Pt,Old_G);
            end
        end
        PrimalFeasi(iter) = Pfeasi;
        Gap(iter)         = gap;
        
        
        %Decompose a small matrix
        if Paras.EvecPast ~= 0 
            [eig_vec,eig_val] = eig(Sstar);
            [~,I]             = sort(diag(eig_val),'descend');
            eig_vec           = eig_vec(:,I);
            eig_val           = eig_val(I,I);
            
            Q1                = eig_vec(:,1:Paras.EvecPast);
            Q2                = eig_vec(:,Paras.EvecPast+1:end);
            Sigma2            = eig_val(Paras.EvecPast+1:end,Paras.EvecPast+1:end);
        end
        

        Z_current = Paras.At'*omegat-Paras.c;

        if DescentFlag == true 
            if iter == 1
                f1 = -Paras.b'*omegat+Paras.rho*max([eig(reshape(Z_current,Paras.n,Paras.n));0]);
            else
                f1 = f3;
            end
            f1_Old = f1;
        else
            f1 = f1_Old;
        end

        Z_next    = Paras.At'*X_next-Paras.c;

        b_inner_X_next = -Paras.b'*X_next;% -b^{t}X_next
        f2 = b_inner_X_next+(Z_next)'*(Wstar);

        EstimatedDrop    = f1 - f2;
        RelativeAccuracy = EstimatedDrop/(abs(f1)+1);

        if EstimatedDrop < 0
            warning('something wrong');
            fprintf('f1-f2 = %.6f',EstimatedDrop);
            break;
        end

        [eig_vec,eig_val] = eig(mat(Z_next));
        [eigval,I]        = sort(diag(eig_val),'descend');
        f3                = b_inner_X_next+Paras.rho*max([eigval(1);0]);
        
        DualSemiFeasi(iter) = min(-eig_val(I(1),I(1)),0);
        SemiFeasi           = DualSemiFeasi(iter);
        RelativeAccry(iter) = RelativeAccuracy;
        
        eig_vec           = eig_vec(:,I);
        Vt                = eig_vec(:,1:Paras.EvecCurrent);
        
        CostDrop  = f1 - f3;
        if Paras.beta*EstimatedDrop <= CostDrop
            %serious step
            omegat = X_next;
            if Paras.adaptive   
                if Paras.mu*EstimatedDrop <= CostDrop
                    Paras.alpha = max(Paras.alpha/2,Paras.alphamin);
                end
                NullCount = 0;
            end
            DescentFlag = true;
        else
            if Paras.adaptive %idle
                NullCount = NullCount+1;
                if Paras.ml*EstimatedDrop >= CostDrop && NullCount >=10 %&& EstimatedDrop > 10^-3
                    Paras.alpha = min(Paras.alpha*2,Paras.alphamax);
                    NullCount = 0;
                end
            end
            DescentFlag = false;
        end

        %Update Wt and Pt
        if Paras.EvecPast ~= 0 
            if Gammastar< 0 
                Gammastar = 0;
            end
            Wt = (Gammastar*Wt+ reshape(Pt*Q2*Sigma2*Q2'*Pt',[],1))/(Gammastar+trace(Sigma2));
        else
            Wt = reshape(Wstar,[],1)/trace(mat(Wstar));
        end
        
        if Paras.EvecPast ~= 0
            %[Pt,~]          = qr([Pt*Q1,Vt],0);
            Pt              = orth([Pt*Q1,Vt]);
        else
            Pt              = Vt;
        end
        Obj(iter)           = f1;


        RelativePFeasi     = norm(Paras.At*Wstar - Paras.b)/(1+normb);
        RelativeGap        = gap/(1+abs(Paras.c.'*Wstar)+abs(Paras.b.'*omegat));
     
        Out.RelativeGap(iter)    = RelativeGap;
        Out.RelativePFeasi(iter) = RelativePFeasi;
        if DescentFlag 
            Out.DescentCost(DescentCount)           = Obj(iter);
            Out.DescentRelativePFeasi(DescentCount) = RelativePFeasi;
            Out.DescentRelativeGap(DescentCount)    = RelativeGap;
            Out.DescentDualSemiFeasi(DescentCount)  = DualSemiFeasi(iter);
            Out.X                                   = Wstar;
            Out.y                                   = omegat;
            Out.Z                                   = Paras.c - Paras.At.'*omegat; 
            DescentCount                            = DescentCount+1;
        end

        %stopping criterion
        if iter > 1 && max([RelativeAccuracy,RelativePFeasi,RelativeGap])< Paras.epislon && SemiFeasi > (-Paras.epislon)
            if ~DescentFlag
                %Move to this new point
                Out.DescentDualSemiFeasi(DescentCount)   = DualFeasi(iter);
                Out.DescentCost(DescentCount)            = f3;%Paras.c.'*omegat - Paras.rho*min([eig(Omegat_sdp);0]);
                Out.DescentRelativePFeasi(DescentCount)  = RelativePFeasi;
                Out.DescentRelativeGap(DescentCount)     = RelativeGap;
                Out.X                                    = Wstar;
                Out.y                                    = omegat;
                Out.Z                                    = Paras.c - Paras.At.'*omegat; 
            end
            fprintf('REACH STOPPING CRITERION!!!\n');
            break;
        end

        if  mod(iter,OutOption.step) == 0 
            fprintf('%5d  %7.2e  %7.2e  %9.2e  %9.2e  %9.2e  %10.2e  %7.1e  %7.2e \n',...
            iter,f1,RelativeAccuracy,RelativePFeasi,0, RelativeGap,SemiFeasi,Paras.alpha,toc(AlgorithmTime));
        end
    end

    Out.Obj         = Obj;
    Out.PrimalFeasi = PrimalFeasi;
    Out.DualFeasi   = DualSemiFeasi;
    Out.Gap         = Gap;
    Out.RelativeAccry = RelativeAccry;
    Out.MainAlgTime   = toc(AlgorithmTime);
    Out.PreprosTime   = Paras.PreprosTime;
    Out.Iter          = iter;
    Out.PastEvec      = opts.EvecPast;
    Out.EvecCurrent   = opts.EvecCurrent;
    Out.DescentCount   = DescentCount-1;

    % Print summary
    if Paras.verbose
        [~,myline1,~] = Header();
        fprintf(myline1);
        fprintf(' SOLUTION SUMMARY:\n');
        fprintf('------------------\n');
        fprintf(' Number of iterations : %11.d\n',iter);
        fprintf(' Cost                 : %11.4e\n',Out.DescentCost(end));
        fprintf(' Relative cost gap    : %11.4e\n',Out.DescentRelativeGap(end));
        fprintf(' Primal residual      : %11.4e\n',Out.DescentRelativePFeasi(end));
        fprintf(' Dual residual        : %11.4e\n',0);
        fprintf(' Prepros time   (s)   : %11.4e\n',Out.PreprosTime);
        fprintf(' SBMP  time   (s)     : %11.4e\n',Out.MainAlgTime);
        fprintf(' Avg. master prob (s) : %11.4e\n',Out.MainAlgTime./iter);
        fprintf(' Total time   (s)     : %11.4e\n',Out.PreprosTime+Out.MainAlgTime);
        fprintf(myline1)
    end
end

function [Paras,OutOption] = Initialize(At,b,c,K,opts)
     % start timing
    proctime = tic;

    %Parameters intilization
    if isfield(opts,'epislon')
        Paras.epislon     = opts.epislon;
    else
        Paras.epislon     = 10^(-4);%default
    end

    Paras.epislon    = opts.epislon;
    if isfield(opts,'beta')
        Paras.beta        = opts.beta;
    else
        Paras.beta        = 0.1; %default
    end

    if isfield(opts,'alpha')
        Paras.alpha       = opts.alpha;
    else
        Paras.alpha       = 1;
    end

    Paras.K = K;
    if ~isfield(K,'f')
        Paras.K.f = 0;
    end
    if Paras.K.f > 0 
        if isfield(opts,'alphafree')
            Paras.alphafree = opts.alphafree;
        else
            Paras.alphafree = Paras.alpha;% default
        end
    end

    if ~isfield(Paras,'Maxiter')
        Paras.Maxiter     = opts.Maxiter;
    else
        Paras.Maxiter     = 500;% default
    end

    Paras.n          = opts.n;
    Paras.m          = opts.m;
    
    if  isfield(opts,'ml')
        Paras.ml          = opts.ml;
    else
        Paras.ml          = 0.001;
    end

    if  isfield(opts,'mu')
        Paras.mu    = opts.mu;
    else
        Paras.mu    = min([Paras.beta*1.5,1]); %default
    end 
    Paras.alphamin   = 10^-5; %for adapative 

    if  isfield(opts,'alphamin')
        Paras.alphamin = opts.alphamin;
    else
        Paras.alphamin    = 10^-5; %for adapative 
    end

    if  isfield(opts,'alphamax')
        Paras.alphamax    = opts.alphamax;
    else
       Paras.alphamax    =  1000; %default
    end     
    
    if ~isfield(opts,'verbose')
        Paras.verbose = true;
    else
        Paras.verbose = opts.verbose;
    end

    if isfield(opts,'adaptive')
        Paras.adaptive = opts.adaptive;
    else
        Paras.adaptive = true;%default
    end
    
    [At,b,c,K,opts] = checkInputs(At,b,c,K,opts);

    Paras.b             = b;
    Paras.At            = At;
    Paras.c             = c;
    Paras.K             = K;
    Paras.MaxCols       = opts.MaxCols;
    Paras.rho           = opts.rho;
    Paras.EvecPast      = opts.EvecPast; 
    Paras.EvecCurrent   = opts.EvecCurrent;
    
    [xIndSym,~,xIndOffDiag,~,~,xIndOffDiagCounter] = SymmetricIndices(Paras.MaxCols,false);
    Paras.IndicesPSD         = xIndSym;
    Paras.IndOffDiagPSD      = xIndOffDiag;
    Paras.IndOffDiagCounter  = xIndOffDiagCounter;
    
    [XIndSym,~,XIndOffDiag,~,~,XIndOffDiagCounter] = SymmetricIndices(Paras.n,false);
    Paras.XIndSym            = XIndSym;
    Paras.XIndOffDiag        = XIndOffDiag;
    Paras.XIndOffDiagCounter = XIndOffDiagCounter;
    
    %Some Operators in the master problem
    Paras.NumOfVar = 1+Paras.MaxCols^2;
    Paras.twoATb    = 2*Paras.At'*Paras.b; %2*A.'*b
    

    %Output option
    OutOption.verbose = 1;
    OutOption.K       = K;
    OutOption.step    = 20;
    OutOption.method  = 'SBMD';
    OutOption.m       = opts.m;
    OutOption.rho     = opts.rho;
    OutOption.past    = opts.EvecPast;
    OutOption.current = opts.EvecCurrent;
    
    Paras.PreprosTime = toc(proctime);
    PrintHeader(proctime,OutOption);
end



