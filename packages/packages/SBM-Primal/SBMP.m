function Out = SBMP(A_sdp,b_sdp,c_sdp,K_sdp,opts)
    % Only compute the top few eigenvectors and eigenvalues. 
    %%%%%%%This version is not stable!!!!
    % Spectral Bundle Method - Primal Formulation
    % Authors: Feng-Yi Liao & Yang Zheng
    %          SOC Lab @UC San Diego
    % Note   : We denote Xnext as X_{t+1} in the paper 
    
    fprintf('Spectral Bundle Method Primal Starts\n');
    [Paras,OutOption] = Initialize(A_sdp,b_sdp,c_sdp,K_sdp,opts);
   
    %Initial point
    if opts.feasible %idle
        [~,Omegat] = InnerApproximation(A_sdp,b_sdp,c_sdp,K_sdp,2,2);
        omegat     = reshape(Omegat,[],1);
    else
        Omegat     = eye(Paras.n);
        %Omegat     = zeros(Paras.n);
        omegat     = reshape(Omegat,[],1);
    end
    
    if opts.adaptive   
        %Paras.alpha = norm(Omegat,'Fro');
        Paras.alpha = 0.1;
        if Paras.alpha == 0
            Paras.alpha  = 1;
        end
        %Paras.alpha = 1;
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
    RelativeAccry     = zeros(opts.Maxiter,1);
    Complementary     = zeros(opts.Maxiter,1);
    RelativePFeasi    = 0;
    RelativeDFeasi    = 0;
    RelativeGap       = 0;
    NullCount         = 0;
    normb             = norm(b_sdp);
    normC             = norm(c_sdp);
    DescentStep       = false;
    
    AlgorithmTime = tic;
    if opts.DynamicRho
        ChangeRho     = false; %A variable that records if Rho has changed or not
    end
    if opts.DynamicMaxCols
        ChangeMaxCols = false;
    end
    
    flag = false;

    for iter = 1:Paras.Maxiter
        %Master problem
        if iter >1 
           if Paras.sparse
               [Wstar,X_next,Gammastar,Sstar,Dfeasi,gap] = Direction_QP_Primal_New(omegat,Paras,Wt,Pt,true);
           else
               [Wstar,X_next,Gammastar,Sstar,Dfeasi,gap,y] = Direction_QP_Primal(omegat,Paras,Wt,Pt,true);
           end
        else
           if Paras.sparse
               [Wstar,X_next,Gammastar,Sstar,Dfeasi,gap] = Direction_QP_Primal_New(omegat,Paras,Wt,Pt,false);
           else
               [Wstar,X_next,Gammastar,Sstar,Dfeasi,gap,y] = Direction_QP_Primal(omegat,Paras,Wt,Pt,false);
           end
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
        
        
        if opts.EvecPast ~=0
            Q1                = eig_vec(:,1:Paras.EvecPast);
            Q2                = eig_vec(:,Paras.EvecPast+1:end);
            Sigma2            = eig_val(Paras.EvecPast+1:end,Paras.EvecPast+1:end);
        end
        
        Omegat = reshape(omegat,Paras.n,Paras.n);
        if issymmetric(Omegat)
            %f1 = c_sdp.'*omegat  + Paras.rho*max([eig(-Omegat);0]);
            f1 = c_sdp.'*omegat  + Paras.rho*max([eigs(-Omegat,1,'la');0]);
        else
            warning('X Not Symmetric');
            %f1 = c_sdp.'*omegat  + Paras.rho*max([eigs(-Omegat,1,'la');0]);
        end
        
        f2 = (c_sdp-Wstar).'*X_next;
        
        EstimatedDrop = f1 - f2;
%         if f2 > 0
%             RelativeAccuracy = EstimatedDrop/f2;
%         elseif f1<0
%             RelativeAccuracy = EstimatedDrop/(-f1);
%         else
%             RelativeAccuracy = EstimatedDrop/(abs(f1)+1);
%         end

        %original
         RelativeAccuracy = EstimatedDrop/(abs(f1)+1);
        
        if EstimatedDrop <0 && iter > 1 
            warning('something wrong');
            fprintf('f1-f2 = %.6f',EstimatedDrop);
            break;
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
            [eig_vec,eig_val] = eigs(Y,Paras.EvecCurrent,'la');
            %[eig_vec,eig_val] = eigs(reshape(Y,Paras.n,[]),Paras.EvecCurrent,'la');
            %[eig_vec,eig_val] = eig(reshape(Y,Paras.n,[]));
            %[eigval,I]             = sort(diag(eig_val),'descend');
           
            f3 = c_sdp.'*X_next + Paras.rho*max([eig_val(1);0]);
%             eig_vec = eig_vec(:,I);
%             Vt      = eig_vec(:,1:Paras.EvecCurrent);
            Vt = eig_vec;
            if opts.EvecPast == 0
                Pt = Vt;
                %[Pt,~]  = qr()
            else
                [Pt,~]  = qr([Pt*Q1,Vt],0);
            end
            
        else
            %warning('Y Not Symmetric');%%%%%%%%%%%%%%%%%%%%%%

            [eig_vec,eig_val] = eigs(Y,Paras.EvecCurrent,'la');
            %[eig_vec,eig_val] = eigs(reshape(Y,Paras.n,[]),Paras.EvecCurrent,'la');
            %[eig_vec,eig_val] = eig(reshape(Y,Paras.n,[]));
            %[eigval,I]             = sort(diag(eig_val),'descend');
           
            f3 = c_sdp.'*X_next + Paras.rho*max([eig_val(1);0]);
%             eig_vec = eig_vec(:,I);
%             Vt      = eig_vec(:,1:Paras.EvecCurrent);
            Vt = eig_vec;
            if opts.EvecPast == 0
                Pt = Vt;
                %[Pt,~]  = qr()
            else
                [Pt,~]  = qr([Pt*Q1,Vt],0);
            end            
        end
        
        CostDrop = f1 - f3;
        if opts.feasible || iter>1
            if Paras.beta*EstimatedDrop <= CostDrop
                %serious step
                omegat = X_next;
                if opts.adaptive 
                    if Paras.mu*EstimatedDrop <= CostDrop
                        Paras.alpha = max(Paras.alpha/2,Paras.alphamin);
                    end
                    NullCount = 0;
                end
            else                
                if opts.adaptive 
                    NullCount = NullCount+1;
                    if Paras.ml*EstimatedDrop >= CostDrop && NullCount >=10 %&& EstimatedDrop > 10^-3
                        Paras.alpha = min(Paras.alpha*2,Paras.alphamax);  
                        %NullCount = 0;
                    end
                end
            end
        else %only for the first iteration
            omegat = X_next;
        end
        

        
        %adjust Maxcols
        if (iter>1 && opts.DynamicMaxCols && ~ChangeMaxCols && (( RelativeAccuracy < opts.epislonphase1 && RelativePFeasi <opts.epislonphase1 && RelativeDFeasi <opts.epislonphase1 && RelativeGap <opts.epislonphase1  ) || iter > 500)) 
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
               % Vt      = eig_vec(:,1:Paras.EvecCurrent);
               [Vt,eig_val] = eigs(Y,Paras.EvecCurrent,'la');
               Pt      = Vt;
             end
        end

        %Update Wt and Pt
        if opts.EvecPast> 0 
            Wt = (Gammastar*Wt+reshape(Pt*Q2*Sigma2*Q2'*Pt.',[],1))/(Gammastar+trace(Sigma2));
        else
            Wt = (Wstar)/(trace(mat(Wstar)));
        end
        
        %improve numerical stability
%         Wt([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) = ...
%             1/2*(Wt([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) + Wt([Paras.XIndOffDiagCounter,Paras.XIndOffDiag]));
                
        
%         [eig_vec,eig_val] = eig(reshape(-X_next,Paras.n,[]));
%         [~,I]             = sort(diag(eig_val),'descend');
        
        if opts.feasible || iter>1
            %PrimalFeasi(iter) = -eig_val(I(1),I(1));
            PrimalFeasi(iter)   = -eig_val(1);
            RelativeAccry(iter) = RelativeAccuracy;
        end
%         eig_vec = eig_vec(:,I);
%         Vt      = eig_vec(:,1:Paras.EvecCurrent);   
%         [Pt,~]  = qr([Pt*Q1,Vt],0);
        if opts.feasible || iter>1
            Obj(iter)            = f1;
            Complementary(iter)  = abs(omegat.'*Wt);
        end
        
        RelativePFeasi     = norm(A_sdp*omegat - b_sdp)/(1+normb);
        RelativeDFeasi     = sqrt(Dfeasi)/(1+normC);
        RelativeGap        = gap/(1+abs(c_sdp.'*omegat)+abs(b_sdp.'*y));


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

%         if ~flag && RelativeAccuracy < 10^-4 && iter>1
%             flag = true;
%             opts.adaptive = false;
%             Paras.alpha = 1;
%             opts.alpha = 1;
%         end


        if iter > 1 && mod(iter,OutOption.step) == 0                 
%            fprintf('%5d  %7.2e  %7.2e  %9.2e  %9.2e  %8.2e  %8.2e  %8.2e %8.2e \n',...
%            iter,f1,RelativeAccuracy,PrimalFeasi(iter),DualFeasi(iter), Gap(iter),Paras.alpha,toc(AlgorithmTime), Complementary(iter));
            fprintf('%5d  %7.2e  %7.2e  %9.2e  %9.2e  %8.2e  %8.2e  %8.2e %8.2e \n',...
            iter,f1,RelativeAccuracy,RelativePFeasi,RelativeDFeasi, RelativeGap,Paras.alpha,toc(AlgorithmTime), Complementary(iter-1));
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
    Out.Complementary  = Complementary;
    Out.DynamicRho     = opts.DynamicRho;
    RelativePFeasi     = norm(A_sdp*omegat - b_sdp)/(1+normb);
    RelativeDFeasi     = sqrt(Dfeasi)/(1+normC);
    RelativeGap        = gap/(1+abs(c_sdp.'*omegat)+abs(b_sdp.'*y));
    Out.RelativeDFeasi = RelativeDFeasi;
    Out.RelativePFeasi = RelativePFeasi;
    Out.RelativeGap    = RelativeGap;
    Out.X              = omegat;
    Out.y              = y;
    Out.Z              = Wstar;
    Out.rho            = opts.rho;
    Out.cost           = Out.Obj(iter);
end

function [Paras,OutOption] = Initialize(At_sdp,b_sdp,c_sdp,K_sdp,opts)
    % start timing
    proctime = tic;
    
    %Parameters intilization
    Paras.epislon     = opts.epislon;
    
    if isfield(opts,'beta')
        Paras.beta        = opts.beta;
    else
        Paras.beta        = 0.1; %default
    end


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
    
    Paras.ml          = 0.001;

    if  isfield(opts,'mu')
        Paras.mu    = opts.mu;
    else
        Paras.mu    = 0.2; %default
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

    [xIndSym,~,xIndOffDiag,~,~,xIndOffDiagCounter] = SymmetricIndices(Paras.MaxCols,false);
    Paras.IndicesPSD        = xIndSym;
    Paras.IndOffDiagPSD     = xIndOffDiag;
    Paras.IndOffDiagCounter = xIndOffDiagCounter;
    
    [XIndSym,~,XIndOffDiag,~,~,XIndOffDiagCounter] = SymmetricIndices(Paras.n,false);
    Paras.XIndSym            = XIndSym;
    Paras.XIndOffDiag        = XIndOffDiag;
    Paras.XIndOffDiagCounter = XIndOffDiagCounter;
%    Paras.At_sdp_shrink      = At_sdp(:,XIndSym);
    
%     [mIndSym,~,mIndOffDiag,~,~,mIndOffDiagCounter] = SymmetricIndices(Paras.m,false);
%     Paras.mIndSym            = mIndSym;
%     Paras.mIndOffDiag        = mIndOffDiag;
%     Paras.mIndOffDiagCounter = mIndOffDiagCounter;
    
    
    if opts.DynamicMaxCols
        [xIndSym2,~,xIndOffDiag2,~,~,xIndOffDiagCounter2] = SymmetricIndices(opts.MaxCols2,false);
        Paras.IndicesPSD2        = xIndSym2;
        Paras.IndOffDiagPSD2     = xIndOffDiag2;
        Paras.IndOffDiagCounter2 = xIndOffDiagCounter2;
    end
    

    %Output option
    OutOption.verbose     = 1;
    OutOption.K           = K_sdp;
    OutOption.step        = 5;
    OutOption.method      = 'SBMP';
    OutOption.m           = opts.m;
    OutOption.rho         = opts.rho;
    OutOption.past    = opts.EvecPast;
    OutOption.current = opts.EvecCurrent;
    
    
    
    %Preproces inv(AAT)
    
    if isfield(opts,'AAT_INV')
        Paras.AAT_INV = sparse(opts.AAT_INV); 
    else
        [Ats,~,~] = svecData(Paras.At_sdp,Paras.c_sdp,Paras.K_sdp);
        AAT       = Ats*Ats.';
        if isdiag(AAT)
            Paras.AAT_INV      = sparse(Paras.m,Paras.m);
            idx                = sub2ind([Paras.m,Paras.m],1:Paras.m,1:Paras.m);
            Paras.AAT_INV(idx) = 1./diag(AAT);
        else
            %Paras.AAT = zeros(Paras.m);
            Paras.AAT_INV = inv(AAT);
            if issparse(Paras.AAT_INV)
                Paras.AAT_INV = sparse(Paras.AAT_INV);
            end
        end
    end

    Paras.Ir2 = eye(Paras.MaxCols^2);
    Paras.q3  = -2*(Paras.At_sdp*Paras.c_sdp);
    %Paras.rowidx = RowSymmetric(Paras.n);
    %Paras.Ir2 = kron(eye(Paras.MaxCols),eye(Paras.MaxCols));
    
    
    %Non-zero elements
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


