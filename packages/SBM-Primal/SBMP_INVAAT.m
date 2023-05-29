function AAT_INV = SBMP_INVAAT(At_sdp,b_sdp,c_sdp,K_sdp)
    % Spectral Bundle Method - Primal Formulation
    % Authors: Feng-Yi Liao & Yang Zheng
    %          SOC Lab @UC San Diego
    % Note   : We denote Xnext as X_{t+1} in the paper 
    % This file only does the inverse of AAT
    
    fprintf('Inverse AAT in process\n');
    %[Paras,~] = Initialize(A_sdp,b_sdp,c_sdp,K_sdp,opts);
    K_sdp.f = 0;
    K_sdp.q = 0;
    K_sdp.l = 0;
    [Ats,~,~] = svecData(At_sdp,c_sdp,K_sdp);
    AAT       = Ats*Ats.';
    m  = height(At_sdp);
    if isdiag(AAT)
        AAT_INV      = sparse(m,m);
        idx                = sub2ind([m,m],1:m,1:m);
        AAT_INV(idx) = 1./diag(AAT);
    else 
        %We first perform choleskey decomposition and do the inversion.
        %This trick is particularlly efficient for sparse matrix 
        L = chol(AAT);
        invL = inv(L);
        AAT_INV = invL*invL.';
    end
end

% function [Paras,OutOption] = Initialize(At_sdp,b_sdp,c_sdp,K_sdp,opts)
%     % start timing
%     proctime = tic;
%     
%     %Parameters intilization
%     Paras.epislon     = opts.epislon;
%     Paras.beta        = opts.beta;
%     Paras.alpha       = opts.alpha;
% %     IndicesAll        = BIGPSDpositionAll(opts.n,opts.dx);%The nonzero indices in a n x n matrix (including lower and upper)
% %     Indices           = BIGPSDposition(opts.n,opts.dx); %The nonzero indices in a n x n matrix (only symmetric part)
% %     Paras.IndicesAll  = IndicesAll;
% %     Paras.Indices     = Indices;
%     Paras.Maxiter     = opts.Maxiter;
% %     Paras.dx               = opts.dx;
%     Paras.n           = opts.n;
%     Paras.m           = opts.m;
%     
%     Paras.ml          = 0.01;
%     Paras.mu          = 0.7;
%     Paras.alphamin    = 10^-5; %for adapative 
%     Paras.alphamax    = 1000;  %for adapative 
%     
% %     Paras.n_new            = 2*Paras.dx;
% %     Paras.NumOfVar_new     = Paras.n_new^2;
% %     Paras.NumOfVar_new_sym = Paras.n_new*(Paras.n_new+1)/2;
% %     Paras.NumOfP           = nchoosek(Paras.n/Paras.dx,2);%Number of Blocks
%     
%     [At_sdp,b_sdp,c_sdp,K_sdp,opts] = checkInputs(At_sdp,b_sdp,c_sdp,K_sdp,opts);
% 
%     Paras.b_sdp            = b_sdp;
%     Paras.At_sdp           = At_sdp;
%     Paras.c_sdp            = c_sdp;
%     Paras.K_sdp            = K_sdp;
%     Paras.MaxCols          = opts.MaxCols;
%     Paras.NumOfVar         = 1+Paras.MaxCols^2;
%     
%     [xIndSym,~,xIndOffDiag,~,~,xIndOffDiagCounter] = SymmetricIndices(Paras.MaxCols,false);
%     Paras.IndicesPSD        = xIndSym;
%     Paras.IndOffDiagPSD     = xIndOffDiag;
%     Paras.IndOffDiagCounter = xIndOffDiagCounter;
%     
%     
%     Paras.rho         = opts.rho;
%     Paras.EvecPast    = opts.EvecPast;
%     Paras.EvecCurrent = opts.EvecCurrent;
%     
%     [XIndSym,~,XIndOffDiag,~,~,XIndOffDiagCounter] = SymmetricIndices(Paras.n,false);
%     Paras.XIndSym            = XIndSym;
%     Paras.XIndOffDiag        = XIndOffDiag;
%     Paras.XIndOffDiagCounter = XIndOffDiagCounter;
% %    Paras.At_sdp_shrink      = At_sdp(:,XIndSym);
%     
% %     [mIndSym,~,mIndOffDiag,~,~,mIndOffDiagCounter] = SymmetricIndices(Paras.m,false);
% %     Paras.mIndSym            = mIndSym;
% %     Paras.mIndOffDiag        = mIndOffDiag;
% %     Paras.mIndOffDiagCounter = mIndOffDiagCounter;
%     
%     %Output option
%     OutOption.verbose = 1;
%     OutOption.K       = K_sdp;
%     OutOption.step    = 1;
%     OutOption.method  = 'SBMP';
%     OutOption.m       = opts.m;
%     OutOption.rho     = opts.rho;
%     OutOption.past    = opts.EvecPast;
%     OutOption.current = opts.EvecCurrent;
%     
%     
%     
%     %Preprocesing
%      
%     
%     [Ats,~,~]     = svecData(Paras.At_sdp,Paras.c_sdp,Paras.K_sdp);
%     AAT       = Ats*Ats.';
%     if isdiag(AAT)
%         Paras.AAT_INV      = sparse(Paras.m,Paras.m);
%         idx                = sub2ind([Paras.m,Paras.m],1:Paras.m,1:Paras.m);
%         Paras.AAT_INV(idx) = 1./diag(AAT);
%     else
%         %Paras.AAT = zeros(Paras.m);
%         Paras.AAT_INV = inv(AAT);
%         Paras.AAT_INV = sparse(Paras.AAT_INV);
%         clear temp
%     end
%     
% 
%     Paras.Ir2 = speye(Paras.MaxCols^2);
%     Paras.q3  = -2*(Paras.At_sdp*Paras.c_sdp);
%     %Paras.rowidx = RowSymmetric(Paras.n);
%     %Paras.Ir2 = kron(eye(Paras.MaxCols),eye(Paras.MaxCols));
%     
%     
%     %Non-zero elements
%     Paras.sparse = opts.sparse;
%     if Paras.sparse
%         [row,idx,val]             = find(Paras.At_sdp);
%         Paras.At_sdp_nonzero.row  = row;
%         Paras.At_sdp_nonzero.idx  = idx;
%         Paras.At_sdp_nonzero.val  = val;
% 
%         [r,c]                     = ind2sub([Paras.n,Paras.n],idx);
%         Paras.At_sdp_nonzero.Mrow = r;
%         Paras.At_sdp_nonzero.Mcol = c;
% 
% 
%         [~,idx]                 = find(Paras.c_sdp);
%         %Paras.c_sdp_nonzero.row   = row;
%         Paras.c_sdp_nonzero.idx   = idx;
%     end
%     
%     Paras.PreprosTime = toc(proctime);
%     PrintHeader(proctime,OutOption);
%     
% end


