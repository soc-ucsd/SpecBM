function Out = Copy_of_SBMP(At,b,c,K,opts)
    % Spectral Bundle Method - Primal Formulation
    % Authors: Feng-Yi Liao & Yang Zheng
    %          SOC Lab @UC San Diego
    % Note   : We denote Xnext as X_{t+1} in the paper 
    
    fprintf('Spectral Bundle Method Primal Starts\n');
    [Paras,OutOption] = Initialize(At,b,c,K,opts);
   
    %Initial point
    if opts.feasible %idle
        [~,Omegat_sdp] = InnerApproximation(At,b,c,K,2,2);
        omegat_sdp     = reshape(Omegat_sdp,[],1);
    else
        Omegat_sdp     = eye(Paras.n);
        %Omegat     = zeros(Paras.n);
        omegat_sdp     = reshape(Omegat_sdp,[],1);
        omegat_free    = zeros(Paras.K.f,1);
        omegat         = [omegat_free;omegat_sdp];
    end
    
    if opts.adaptive   
        %Paras.alpha = norm(Omegat,'Fro');
        Paras.alpha = 0.5;
        if Paras.alpha == 0
            Paras.alpha  = 1;
        end
        %Paras.alpha = 1;
    end
    
    %Initialize Pt and Wt
    [eig_vec,eig_val] = eig(mat(-omegat_sdp));
    [~,I]             = sort(diag(eig_val),'descend');
    eig_vec           = eig_vec(:,I);
    Pt                = eig_vec(:,1:opts.MaxCols);  %transformation matrix
    Wt                = reshape(eig_vec(:,1)*eig_vec(:,1)',[],1);
    
    Obj               = zeros(opts.Maxiter,1); 
    DualFeasi         = zeros(opts.Maxiter,1); 
    PrimalFeasi       = zeros(opts.Maxiter,1); 
    Gap               = zeros(opts.Maxiter,1);
    RelativeAccry     = zeros(opts.Maxiter,1);
    %Complementary     = zeros(opts.Maxiter,1);
    RelativePFeasi    = 0;
    RelativeDFeasi    = 0;
    RelativeGap       = 0;
    NullCount         = 0;
    normb             = norm(Paras.b);
    normC             = norm(Paras.c);
    if norm(Paras.c_free)~= 0
        normc_free            = norm(Paras.c_free);
    end

    Out.RelativeDFeasi = [];
    Out.RelativePFeasi = [];
    Out.RelativeGap    = [];

    Out.DescentRelativeDFeasi  = [];
    Out.DescentRelativePFeasi  = [];
    Out.DescentRelativeGap     = [];
    Out.DescentCost            = [];
    Out.DescentPrimalSemiFeasi = [];


    DescentCount  = 1;

    AlgorithmTime = tic;
    if opts.DynamicRho
        ChangeRho     = false; %A variable that records if Rho has changed or not
    end
    if opts.DynamicMaxCols %If or not changing the number of selected eigenvectors after some threshold
        ChangeMaxCols = false;
    end
    

    for iter = 1:Paras.Maxiter
        %Master problem
        if iter >1 
           if Paras.sparse %idle 
               [Wstar,X_next,Gammastar,Sstar,Dfeasi,gap] = Direction_QP_Primal_New(omegat_sdp,Paras,Wt,Pt,true);
           else
               %[Wstar,X_next,Gammastar,Sstar,Dfeasi,gap,y] = Direction_QP_Primal(omegat,Paras,Wt,Pt,true);
               [Wstar,X_next,Gammastar,Sstar,Dfeasi,Dfeasi_free,gap,y,x_next] = Direction_QP_Primal_Free(omegat_free,omegat_sdp,Paras,Wt,Pt,true);
           end
        else
           if Paras.sparse  %idle 
               [Wstar,X_next,Gammastar,Sstar,Dfeasi,gap] = Direction_QP_Primal_New(omegat_sdp,Paras,Wt,Pt,false);
           else
               %[Wstar,X_next,Gammastar,Sstar,Dfeasi,gap,y] = Direction_QP_Primal(omegat,Paras,Wt,Pt,false);
               [Wstar,X_next,Gammastar,Sstar,Dfeasi,Dfeasi_free,gap,y,x_next] = Direction_QP_Primal_Free(omegat_free,omegat_sdp,Paras,Wt,Pt,false);
           end
        end
        
        %We do not record the first iteration if the initial point is not feasible
        if opts.feasible || iter>1
            DualFeasi(iter) = max(Dfeasi,Dfeasi_free);
            Gap(iter) = gap;
        end
        
        
        %Decompose a small matrix
        [eig_vec,eig_val] = eig(full(Sstar));
        [~,I]             = sort(diag(eig_val),'descend');
        eig_vec           = eig_vec(:,I);
        eig_val           = eig_val(I,I);
        
        
        if opts.EvecPast ~=0
            Q1                = eig_vec(:,1:Paras.EvecPast);
            Q2                = eig_vec(:,Paras.EvecPast+1:end);
            Sigma2            = eig_val(Paras.EvecPast+1:end,Paras.EvecPast+1:end);
        end
        
        Omegat_sdp = reshape(omegat_sdp,Paras.n,Paras.n);
        if issymmetric(Omegat_sdp)
%             f1 =Paras.c0.'*omegat_free + Paras.c_sdp.'*omegat_sdp  + Paras.rho*max([eig(-Omegat_sdp);0]);
            %f1 = Paras.c.'*omegat + Paras.rho*max([eig(-Omegat_sdp);0]);

            f1 = Paras.c.'*omegat - Paras.rho*min([eig(Omegat_sdp);0]);

%             eigs(Omegat_sdp,1)
%             eigs(-Omegat_sdp,1)
            %[T Q] = myLanczos(-Omegat_sdp,5);

%            [v, xi, i] = ApproxMinEvecLanczos(Omegat_sdp, Paras.n, 50);
%            f1 = Paras.c.'*omegat + Paras.rho*max([-xi,0]);
            %f1 = Paras.c.'*omegat + Paras.rho*max([eigs(-(Omegat_sdp),1,'largestabs');0]);
        else
            warning('X Not Symmetric');
        end
        
        f2 = Paras.c_free.'*x_next+(Paras.c_sdp-Wstar).'*X_next;


        EstimatedDrop = f1 - f2;

        %original
        RelativeAccuracy = EstimatedDrop/(abs(f1)+1);
        
        if EstimatedDrop <0 && iter > 1 
%             warning('something wrong');
%             fprintf('f1-f2 = %.6f',EstimatedDrop);
            %break;
        end
        %elseif RelativeAccuracy < Paras.epislon && iter >1 && PrimalFeasi(iter) < 10^-3 && DualFeasi(iter) < 10^-3
%         if RelativeAccuracy < Paras.epislon && iter >1  && PrimalFeasi(iter) < 10^-3 && DualFeasi(iter) < 10^-3 && Complementary(iter) < 10^-3
%             Obj(iter)           = f1;
%             PrimalFeasi(iter)   = eig_val(I(1),I(1));
%             RelativeAccry(iter) = RelativeAccuracy;
%             break;
%         end
        
        
        %Serious Step or Null Step        
        Y = reshape(-X_next,Paras.n,Paras.n);
        if issymmetric(Y)
            [eig_vec,eig_val] = eig(reshape(Y,Paras.n,[])); %This step can be improved to only computing the top r_c number of eigenvectors
            [eigval,I]        = sort(diag(eig_val),'descend');
            f3                = Paras.c_free.'*x_next + Paras.c_sdp.'*X_next + Paras.rho*max([eigval(1);0]);
            eig_vec           = eig_vec(:,I);
            Vt                = eig_vec(:,1:Paras.EvecCurrent);

%             [Vt,eigval] = eigs(reshape(Y,Paras.n,[]),Paras.EvecCurrent,'largestreal');
%             f3           = Paras.c_free.'*x_next + Paras.c_sdp.'*X_next + Paras.rho*max([eigval(1);0]);
            if opts.EvecPast == 0
                Pt = Vt;
            else
                %[Pt,~]  = qr([Pt*Q1,Vt],0);
                OldPt = Pt;
                Pt = orth([Pt*Q1,Vt],1e-12);
            end 
        else
            warning('Y Not Symmetric');
        end
        
        CostDrop = f1 - f3;
        if opts.feasible || iter>1
            if Paras.beta*EstimatedDrop <= CostDrop
                %serious step
                omegat_sdp              = X_next;
                omegat_free             = x_next;
                omegat(1:Paras.K.f)     = omegat_free;
                omegat(Paras.K.f+1:end) = omegat_sdp;
                if opts.adaptive 
                    if Paras.mu*EstimatedDrop <= CostDrop
                        Paras.alpha = max(Paras.alpha/2,Paras.alphamin);
                    end
                    NullCount = 0;
                end
                DescentFlag = true;
            else                
                if opts.adaptive 
                    NullCount = NullCount+1;
                    if Paras.ml*EstimatedDrop >= CostDrop && NullCount >=10 %&& EstimatedDrop > 10^-3
                        Paras.alpha = min(Paras.alpha*2,Paras.alphamax);  
                        NullCount = 0;
                    end
                end
                DescentFlag = false;
            end
        else %only for the first iteration
            omegat_sdp        = X_next;
            omegat_free       = x_next;
            omegat(1:Paras.K.f)     = omegat_free;
            omegat(Paras.K.f+1:end) = omegat_sdp;
            DescentFlag = true;
        end

%         norm(Paras.b_sdp-Paras.At_sdp*omegat-Paras.a*yomegat)
        
        %adjust Maxcols
        if (iter>1 && opts.DynamicMaxCols && ~ChangeMaxCols && (( RelativeAccuracy < opts.epislonphase1 && RelativePFeasi < opts.epislonphase1 && RelativeDFeasi < opts.epislonphase1 && RelativeGap <opts.epislonphase1 || iter >500 ) )) 
             Paras.MaxCols            = opts.MaxCols2;
             Paras.EvecPast           = opts.EvecPast2;
             Paras.EvecCurrent        = opts.EvecCurrent2;
             Paras.IndicesPSD         = Paras.IndicesPSD2;
             Paras.IndOffDiagPSD      = Paras.IndOffDiagPSD2;
             Paras.IndOffDiagCounter  = Paras.IndOffDiagCounter2;
             ChangeMaxCols            = true;
             Paras.Ir2                = eye(Paras.MaxCols^2);
             Paras.NumOfVar           = 1+Paras.MaxCols^2;
             fprintf('MaxCols change \n');
             Q1                = eig_vec(:,1:Paras.EvecPast);
             Q2                = eig_vec(:,Paras.EvecPast+1:end);
             Sigma2            = eig_val(Paras.EvecPast+1:end,Paras.EvecPast+1:end);
             if opts.EvecPast == 0
                Vt      = eig_vec(:,1:Paras.EvecCurrent);
                Pt      = Vt;
             end
        end
        
        %Update Wt and Pt
        if opts.EvecPast> 0 
            if Gammastar< 0 
                Gammastar = 0;
            end
            %Wt = (Gammastar*Wt+reshape(Pt*Q2*Sigma2*Q2'*Pt.',[],1))/(Gammastar+trace(Sigma2));
            Wt = (Gammastar*Wt+reshape(OldPt*Q2*Sigma2*Q2'*OldPt.',[],1))/(Gammastar+trace(Sigma2));
        else
            Wt = (Wstar)/(trace(mat(Wstar)));
        end
        
        %improve numerical stability
%         Wt([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) = ...
%             1/2*(Wt([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) + Wt([Paras.XIndOffDiagCounter,Paras.XIndOffDiag]));
                
        
%         [eig_vec,eig_val] = eig(reshape(-X_next,Paras.n,[]));
%         [~,I]             = sort(diag(eig_val),'descend');
        
        if opts.feasible || iter>1
            PrimalFeasi(iter) = -eig_val(I(1),I(1));
            RelativeAccry(iter) = RelativeAccuracy;
        end
%         eig_vec = eig_vec(:,I);
%         Vt      = eig_vec(:,1:Paras.EvecCurrent);   
%         [Pt,~]  = qr([Pt*Q1,Vt],0);
        if opts.feasible || iter>1
            Obj(iter)            = f1;
            %Complementary(iter)  = abs(omegat_sdp.'*Wt);
        end
        
        RelativePFeasi     = norm(Paras.At*omegat - Paras.b)/(1+normb);
        if norm(Paras.c_free)~= 0
            RelativeDFeasi = max(sqrt(Dfeasi)/(1+normC),sqrt(Dfeasi_free)/(1+normc_free));
        else
            RelativeDFeasi = sqrt(Dfeasi)/(1+normC);
        end
        RelativeGap        = gap/(1+abs(Paras.c.'*omegat)+abs(Paras.b.'*y));
        
        if iter > 1
            temp = iter-1;
            Out.RelativeDFeasi(temp) = RelativeDFeasi;
            Out.RelativeGap(temp)    = RelativeGap;
            Out.RelativePFeasi(temp) = RelativePFeasi; 
            if DescentFlag 
                Out.DescentPrimalSemiFeasi(DescentCount) = PrimalFeasi(iter);
                Out.DescentCost(DescentCount)            = Obj(iter);
                Out.DescentRelativeDFeasi(DescentCount)  = RelativeDFeasi;
                Out.DescentRelativePFeasi(DescentCount)  = RelativePFeasi;
                Out.DescentRelativeGap(DescentCount)     = RelativeGap;
                Out.X                                    = omegat_sdp;
                Out.y                                    = y;
                Out.Z                                    = Wstar; 
                DescentCount                             = DescentCount+1;
            end
        end

        %stopping criterion
        if iter > 1 && RelativeAccuracy < Paras.epislon && RelativePFeasi <Paras.epislon && RelativeDFeasi< Paras.epislon &&  RelativeGap  < Paras.epislon
%        if iter > 1 && RelativeAccuracy < Paras.epislon && RelativePFeasi(iter) < 10^-3 && RelativeDFeasi(iter) < 10^-3 &&  RelativeGap(iter)   < 10^-3
%        if RelativeAccuracy < Paras.epislon && iter >1  && PrimalFeasi(iter) > -Paras.epislon && DualFeasi(iter) < Paras.epislon && Complementary(iter) < Paras.epislon
%             Obj(iter)           = f1;
%             PrimalFeasi(iter)   = eig_val(I(1),I(1));
%             RelativeAccry(iter) = RelativeAccuracy;
            break;
        end


%         %adjusting rho
%         if iter>1 && Paras.DynamicRho && ~ChangeRho && PrimalFeasi(iter) > -10^-2
%             Paras.rho = Paras.rho2; 
%             fprintf('rho change to %f \n', Paras.rho2);
%             ChangeRho = true;
%         end


        if iter > 1 && mod(iter,OutOption.step) == 0                 
%            fprintf('%5d  %7.2e  %7.2e  %9.2e  %9.2e  %8.2e  %8.2e  %8.2e %8.2e \n',...
%            iter,f1,RelativeAccuracy,PrimalFeasi(iter),DualFeasi(iter), Gap(iter),Paras.alpha,toc(AlgorithmTime), Complementary(iter));
            fprintf('%5d  %7.2e  %7.2e  %9.2e  %9.2e  %8.2e  %8.2e  %8.2e \n',...
            iter,f1,RelativeAccuracy,RelativePFeasi,RelativeDFeasi, RelativeGap,Paras.alpha,toc(AlgorithmTime));
        end
    end
    Out.Obj            = Obj;
    Out.PrimalFeasi    = PrimalFeasi;
    Out.DualFeasi      = DualFeasi;
    Out.Gap            = Gap;
    Out.RelativeAccry  = RelativeAccry;
    Out.MainAlgTime    = toc(AlgorithmTime);
    Out.PreprosTime    = Paras.PreprosTime;
    Out.Iter           = iter;
    Out.PastEvec       = opts.EvecPast;
    Out.EvecCurrent    = opts.EvecCurrent;
    %Out.Complementary  = Complementary;
    Out.DynamicRho     = opts.DynamicRho;

%     RelativePFeasi     = norm(Paras.At*omegat - Paras.b)/(1+normb);
%     RelativeDFeasi     = sqrt(Dfeasi)/(1+normC);
%     RelativeGap        = gap/(1+abs(Paras.c.'*omegat)+abs(Paras.b.'*y));
%     Out.RelativeDFeasi = RelativeDFeasi;
%     Out.RelativePFeasi = RelativePFeasi;
%     Out.RelativeGap    = RelativeGap;
    
%     Out.X              = omegat_sdp;
%     Out.y              = y;
%     Out.Z              = Wstar;
    Out.rho            = opts.rho;
    Out.cost           = Out.Obj(iter);
    Out.DescentCount   = DescentCount-1;
    % Print summary
    if Paras.verbose
        [~,myline1,~] = Header();
        fprintf(myline1);
        fprintf(' SOLUTION SUMMARY:\n');
        fprintf('------------------\n');
        %fprintf(' Termination code     : %11.1d\n',info.problem)
        fprintf(' Number of iterations : %11.d\n',iter);
        fprintf(' Cost                 : %11.4e\n',Out.DescentCost(end));
        fprintf(' Relative cost gap    : %11.4e\n',Out.DescentRelativeGap(end));
        fprintf(' Primal residual      : %11.4e\n',Out.DescentRelativePFeasi(end));
        fprintf(' Dual residual        : %11.4e\n',Out.DescentRelativeDFeasi(end));
        fprintf(' Prepros time   (s)   : %11.4e\n',Out.PreprosTime);
        fprintf(' SBMP  time   (s)     : %11.4e\n',Out.MainAlgTime);
        %fprintf(' Avg. conic proj (s)  : %11.4e\n',info.time.subiter(2)./iter)
        fprintf(' Avg. master prob (s) : %11.4e\n',Out.MainAlgTime./iter);
        %fprintf(' Cleanup time (s)     : %11.4e\n',posttime)
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
    
    Paras.n           = opts.n;
    Paras.m           = opts.m;

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
    
    if  isfield(opts,'alphamin')
        Paras.alphamin = opts.alphamin;
    else
        Paras.alphamin    = 10^-5; %for adapative 
    end
    
    if  isfield(opts,'alphamax')
        Paras.alphamax    = opts.alphamax;
    else
        Paras.alphamax    = 1000;  %for adapative 
    end

    if ~isfield(opts,'verbose')
        Paras.verbose = true;
    else
        Paras.verbose = opts.verbose;
    end
    
%     Paras.n_new            = 2*Paras.dx;
%     Paras.NumOfVar_new     = Paras.n_new^2;
%     Paras.NumOfVar_new_sym = Paras.n_new*(Paras.n_new+1)/2;
%     Paras.NumOfP           = nchoosek(Paras.n/Paras.dx,2);%Number of Blocks
    
    [At,b,c,K,opts] = checkInputs(At,b,c,K,opts);

    if  ~isfield(opts,'rescale')
        opts.rescale = true;
    end

    %[At,b,c,K,opts] = rescaleData(At',b,c,K,opts);
    %At = At.';
    
    Paras.b                = b;
    Paras.At               = At;
    Paras.At_free          = At(:,1:K.f);
    Paras.At_sdp           = At(:,K.f+1:end);
    Paras.c                = c;
    Paras.c_free           = c(1:K.f);
    Paras.c_sdp            = c(K.f+1:end);
    Paras.K                = K;
    Paras.MaxCols          = opts.MaxCols;
    Paras.NumOfVar         = 1+Paras.MaxCols^2;

    
    
    %idle now
    Paras.DynamicRho       = opts.DynamicRho;
    if Paras.DynamicRho
        %rho start at rho1 and change rho2 when primal feasibility is low
        Paras.rho1        = opts.rho1;
        Paras.rho2        = opts.rho2;
        Paras.rho         = Paras.rho1;
    else%rho does not change
        Paras.rho         = opts.rho; 
    end
    
    Paras.EvecPast    = opts.EvecPast;
    Paras.EvecCurrent = opts.EvecCurrent;

    %This generate the indices of the lower triangle parts of a symmetric
    % matrix
    [xIndSym,~,xIndOffDiag,~,~,xIndOffDiagCounter] = SymmetricIndices(Paras.MaxCols,false);
    Paras.IndicesPSD        = xIndSym;
    Paras.IndOffDiagPSD     = xIndOffDiag;
    Paras.IndOffDiagCounter = xIndOffDiagCounter;
    
    [XIndSym,~,XIndOffDiag,~,~,XIndOffDiagCounter] = SymmetricIndices(Paras.n,false);
    Paras.XIndSym            = XIndSym;
    Paras.XIndOffDiag        = XIndOffDiag;
    Paras.XIndOffDiagCounter = XIndOffDiagCounter;    
    
    

    if opts.DynamicMaxCols
        [xIndSym2,~,xIndOffDiag2,~,~,xIndOffDiagCounter2] = SymmetricIndices(opts.MaxCols2,false);
        Paras.IndicesPSD2        = xIndSym2;
        Paras.IndOffDiagPSD2     = xIndOffDiag2;
        Paras.IndOffDiagCounter2 = xIndOffDiagCounter2;
    end
    

    %Output option
    OutOption.verbose     = 1;
    OutOption.K           = K;
    OutOption.step        = 10;
    OutOption.method      = 'SBMP';
    OutOption.m           = opts.m;
    OutOption.rho         = opts.rho;
    OutOption.past        = opts.EvecPast;
    OutOption.current     = opts.EvecCurrent;
    
    %Preproces inv(AAT)
    
    if isfield(opts,'AAT_INV')
        Paras.AAT_INV = sparse(opts.AAT_INV); 
    else
        %[Ats,~,~] = svecData(Paras.At_sdp,Paras.c_sdp,Paras.K_sdp);
        [Ats,~,~] = svecData(At,c,K);
        AAT       = Ats*Ats.';
        if isdiag(AAT) 
            Paras.AAT_INV      = sparse(Paras.m,Paras.m);
            idx                = sub2ind([Paras.m,Paras.m],1:Paras.m,1:Paras.m);
            Paras.AAT_INV(idx) = 1./diag(AAT);
        else
            L = chol(AAT);
            invL = inv(L);
            %chol_AAT = chol(AAT);
            %Paras.AAT_INV = inv(AAT);
            Paras.AAT_INV = invL*invL.';
            if issparse(Paras.AAT_INV)
                Paras.AAT_INV = sparse(Paras.AAT_INV);
            end
        end
    end

    Paras.Ir2 = eye(Paras.MaxCols^2);
    Paras.q3  = -2*(Paras.At_sdp*Paras.c_sdp);
    
    %Non-zero elements %not necessary now, because it does not speed up the
    %computation
    Paras.sparse = opts.sparse;
    if Paras.sparse
        [row,idx,val]             = find(Paras.At_sdp);
        Paras.At_sdp_nonzero.row  = row;
        Paras.At_sdp_nonzero.idx  = idx;
        Paras.At_sdp_nonzero.val  = val;

        [r,c]                     = ind2sub([Paras.n,Paras.n],idx);
        Paras.At_sdp_nonzero.Mrow = r;
        Paras.At_sdp_nonzero.Mcol = c;

        [~,idx]                 = find(Paras.c_sdp);
        %Paras.c_sdp_nonzero.row   = row;
        Paras.c_sdp_nonzero.idx   = idx;
    end
    
    Paras.PreprosTime = toc(proctime);
    PrintHeader(proctime,OutOption);
    
end
%% Lanczos method
function [v, xi, i] = ApproxMinEvecLanczos(M, n, q)
    % Approximate minimum eigenvector
    % Vanilla Lanczos method
    
    q = min(q, n-1);                    % Iterations < dimension!
    
    if isnumeric(M), M = @(x) M*x; end
    
    Q = zeros(n, q+1);                  % Lanczos vectors
    
    aleph = zeros(q,1);                 % Diagonal Lanczos coefs
    beth = zeros(q,1);                  % Off-diagonal Lanczos coefs
    
    Q(:,1) = randn(n, 1);               % First Lanczos vector is random
    Q(:,1) = Q(:,1) / norm(Q(:,1));
    
    for i = 1 : q
        Q(:, i+1) = M ( Q(:, i) );				% Apply M to previous Lanczos vector
        aleph(i) = real(Q(:, i)' * Q(:, i+1));		% Compute diagonal coefficients
        
        if (i == 1)                     % Lanczos iteration
            Q(:, i+1) = Q(:, i+1) - aleph(i) * Q(:, i);
        else
            Q(:, i+1) = Q(:, i+1) - aleph(i) * Q(:, i) - beth(i-1) * Q(:, i-1);
        end
        
        beth(i) = norm( Q(:, i+1) );            % Compute off-diagonal coefficients
        
        if ( abs(beth(i)) < sqrt(n)*eps ), break; end
        
        Q(:, i+1) = Q(:, i+1) / beth(i);        % Normalize
        
    end
    
    % i contains number of completed iterations
    
    B = diag(aleph(1:i), 0) + diag(beth(1:(i-1)), +1) + diag(beth(1:(i-1)), -1);
    
    [U, D] = cgal_eig(0.5*(B+B'));
    [xi, ind] = min(D);
    v = Q(:, 1:i) * U(:, ind);
    
    % Next lines are unnecessary in general, but I observed numerical errors in
    % norm(v) at some experiments, so let's normalize it for robustness. 
    nv = norm(v);
    xi = xi*nv;
    v = v/nv;
end
function [V,D] = cgal_eig(X)
% Eig in Lanczos based LMO solver sometimes fall into numerical issues. 
% This function replaces eig with a SVD based solver, in case eig does not
% converge. 
try
    [V,D]       = eig(X,'vector');
catch 
	warning('eig did not work. Using the svd based replacement instead.');
    [V,D,W]     = svd(X);
    D           = diag(D).' .* sign(real(dot(V,W,1)));
    [D,ind]     = sort(D);
    V           = V(:,ind);
end
end


