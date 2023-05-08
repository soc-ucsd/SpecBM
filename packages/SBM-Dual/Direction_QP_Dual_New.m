function [Wstar,X_next,Gammastar,Sstar,PrimalFeasibility,gap] = Direction_QP_Dual_New(omegat,Paras,Wt,Pt) 
    %%%%%%%%%%deprecated%%%%%%%%%%%%%%
    %Authors: Feng-Yi Liao & Yang Zheng
    %         SOC Lab @UC San Diego
    %Wt is a fixed atoms
    %Pt is the transformation matrix
        
    %G = -Paras.At_sdp'*omegat;
    
    AW1 = zeros(Paras.m,1);
    G  = zeros(Paras.n^2,1);
    for i = 1:length(Paras.At_sdp_nonzero.row) 
        idx    = Paras.At_sdp_nonzero.idx(i);
        row    = Paras.At_sdp_nonzero.row(i);
        val    = Paras.At_sdp_nonzero.val(i);
        G(idx) = G(idx) - val*omegat(row);
        AW1(row) = AW1(row) + Paras.At_sdp(row,idx)*Wt(idx);
    end
    G = G + Paras.c_sdp;
    
    %AW1        = Paras.At_sdp*Wt;
    M11        = (AW1)'*(AW1);
    %kronPTPT   = kron(Pt',Pt');     
    %kronPTPTAT = kronPTPT*Paras.At_sdp';
    %M21        = kronPTPTAT*AW1;
    
    %M21        = zeros(Paras.MaxCols,1);
    AAW1 = zeros(Paras.n^2,1);
    for i = 1:length(Paras.At_sdp_nonzero.row) 
        idx  = Paras.At_sdp_nonzero.idx(i);
        row  = Paras.At_sdp_nonzero.row(i);
        val    = Paras.At_sdp_nonzero.val(i);
        AAW1(idx) = AAW1(idx) + val*AW1(row); 
    end
    M21        = vec(Pt'*mat(AAW1)*Pt);
    %M21        = vec(Pt'*mat(Paras.At_sdp'*AW1)*Pt);
    
%    M22        = zeros(Paras.MaxCols^2);
%     for i = 1:Paras.m
%        v    = vec(Pt'*mat(Paras.At_sdp(i,:))*Pt);
%        M22  = M22 + v*v';
%     end
    
    Ptt = Pt';
    V = zeros(Paras.MaxCols^2,Paras.m);
    for i = 1:length(Paras.At_sdp_nonzero.row) 
        idx  = Paras.At_sdp_nonzero.idx(i);
        row  = Paras.At_sdp_nonzero.row(i);
        Mrow = Paras.At_sdp_nonzero.Mrow(i);
        Mcol = Paras.At_sdp_nonzero.Mcol(i); 
        V(:,row) = V(:,row) + kron(Ptt(:,Mcol),Ptt(:,Mrow)); 
    end
    M22 = V*V';
    
    
    %M22        = kronPTPTAT*kronPTPTAT';
    %M22        = (M22+M22')/2;
    
    %for numerical stability
    minlambda = min(eig(M22));
    if minlambda <0
        M22   = M22-(minlambda-10^(-9))*eye(Paras.MaxCols^2);
    end
  
    m1 = (-Paras.twoAb+2*Paras.alpha*G)'*Wt;
    m2 = vec(Pt'*mat(-Paras.twoAb+2*Paras.alpha*G)*Pt);
    M  = [M11,M21';
         M21,M22];
    M  = (M+M')/2;
    [eig_vec,eig_val]  = eig(M) ;
    %just to avoid numerical error (complex number)
    eig_val(eig_val<0) = 0;
    M05                = eig_vec*sqrt(eig_val)*eig_vec'; %B = M^{1/2}
    M05                = (M05+M05.')/2;
   
    %%%%%%%%%%%   VERY IMPORTANT!!!!!!!!!!!!   %%%%%%%%%%%%%%%%%%%%
    M05(:,1+[Paras.IndOffDiagPSD;Paras.IndOffDiagCounter]) = ...
      1/2*(M05(:,1+[Paras.IndOffDiagPSD;Paras.IndOffDiagCounter])+M05(:,1+[Paras.IndOffDiagCounter;Paras.IndOffDiagPSD]));
   
   %sedumi constraint matrix
  % toc(test)
   
   %Ir2 = eye(Paras.MaxCols^2);
   At                                             = zeros(1+Paras.MaxCols^2+2,2*Paras.NumOfVar+3);
   At(1:Paras.NumOfVar,2)                         = M05(:,1);
   At(1:Paras.NumOfVar,4+Paras.NumOfVar+1:end)    = M05(:,2:end);
   At(1:Paras.NumOfVar,5:4+Paras.NumOfVar)        = -eye(Paras.NumOfVar);
   At(Paras.NumOfVar+1,1)                         = 1;
   At(Paras.NumOfVar+1,2)                         = 1;
   At(Paras.NumOfVar+1,end-Paras.MaxCols^2+1:end) = reshape(eye(Paras.MaxCols),1,[]);
   At(end,4)                                      = 1;
   
   b                   = zeros(Paras.NumOfVar+2,1);
   b(Paras.NumOfVar+1) = Paras.rho;
   b(end)              = 0.5;
   c                   = [0;m1;1;0;zeros(Paras.NumOfVar,1);m2];
   K.s                 = [Paras.MaxCols];
   K.r                 = 2+Paras.NumOfVar;
   K.l                 = 2;
   
   %Call mosek to solve the master problem
   prob1     = SedumiToMosek_Latest(At,b,c,K);
   [~, res1] = mosekopt('minimize echo(0)', prob1);
   status    = res1.sol.itr.prosta;
   if ~strcmp(status,'PRIMAL_AND_DUAL_FEASIBLE')
      warning('infeasible');
      return;
   end
   
   xlinear                        = res1.sol.itr.xx;
   xPSD                           = res1.sol.itr.barx;

   Gammastar                      = xlinear(2);
   Sstar                          = zeros(Paras.MaxCols);
   Sstar(Paras.IndicesPSD)        = xPSD;
   Sstar(Paras.IndOffDiagCounter) = Sstar(Paras.IndOffDiagPSD);
   Wstar                          = Gammastar*Wt + reshape(Pt*Sstar*Pt',[],1);    
    
  %improve numerical stability
    Wstar([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) = ...
          1/2*(Wstar([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) + Wstar([Paras.XIndOffDiagCounter,Paras.XIndOffDiag]));
    
    PrimalAffine      = Paras.b_sdp;
    for i = 1:length(Paras.At_sdp_nonzero.row) 
        idx  = Paras.At_sdp_nonzero.idx(i);
        row  = Paras.At_sdp_nonzero.row(i);
        PrimalAffine(row) = PrimalAffine(row) - Paras.At_sdp(row,idx)*Wstar(idx);
    end
    %PrimalAffine      = Paras.b_sdp-Paras.At_sdp*vec(Wstar);
    X_next            = omegat + (PrimalAffine)/Paras.alpha;
    PrimalFeasibility = norm(PrimalAffine)^2;
    gap               = Paras.b_sdp'*X_next - Paras.c_sdp'*reshape(Wstar,[],1);
end