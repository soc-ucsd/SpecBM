function [Wstar,X_next,Gammastar,Sstar,PrimalFeasibility,gap,G] = Direction_QP_Dual(omegat,Paras,Wt,Pt,Old_G)   
    %Authors: Feng-Yi Liao & Yang Zheng
    %         SOC Lab @UC San Diego
    %Wt is a fixed atoms
    %Pt is the transformation matrix
    
    switch nargin
        case 4
            G = -Paras.At.'*omegat;
            G = G + Paras.c;
        case 5 
            G = Old_G;
    end

    AW1        = Paras.At*Wt;
    M11        = (AW1).'*(AW1);
    kronPTPT   = kron(Pt',Pt');     
    kronPTPTAT = kronPTPT*Paras.At.';
    M21        = kronPTPTAT*AW1;
    M22        = kronPTPTAT*kronPTPTAT.';
    M22        = (M22+M22')/2;
    
    %for numerical stability
    minlambda = min(eig(M22));
    if minlambda <0
        M22   = M22-(minlambda-10^(-9))*eye(Paras.MaxCols^2);
    end
  
    m1 = (-Paras.twoATb+2*Paras.alpha*G)'*Wt;
    m2 = kronPTPT*(-Paras.twoATb+2*Paras.alpha*G);
    M  = [M11,M21.';
         M21,M22];
    M  = (M+M.')/2;
    [eig_vec,eig_val]  = eig(M) ;
    %just to avoid numerical error (complex number)
    eig_val(eig_val<0) = 0;
    M05                = eig_vec*sqrt(eig_val)*eig_vec.'; %B = M^{1/2}
    M05                = (M05+M05.')/2;
   
    %%%%%%%%%%%   VERY IMPORTANT!!!!!!!!!!!!   %%%%%%%%%%%%%%%%%%%%
    M05(:,1+[Paras.IndOffDiagPSD;Paras.IndOffDiagCounter]) = ...
      1/2*(M05(:,1+[Paras.IndOffDiagPSD;Paras.IndOffDiagCounter])+M05(:,1+[Paras.IndOffDiagCounter;Paras.IndOffDiagPSD]));
   
   %sedumi constraint matrix
  % toc(test)
   
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
   if (~strcmp(status,'UNKNOWN') && ~strcmp(status,'PRIMAL_AND_DUAL_FEASIBLE'))
      warning('infeasible');
      return;
   end
   
   xlinear                        = res1.sol.itr.xx;
   xPSD                           = res1.sol.itr.barx;

   Gammastar                      = xlinear(2);
   Sstar                          = zeros(Paras.MaxCols);
   Sstar(Paras.IndicesPSD)        = xPSD;
   Sstar(Paras.IndOffDiagCounter) = Sstar(Paras.IndOffDiagPSD);
   Wstar                          = Gammastar*Wt + reshape(Pt*Sstar*Pt.',[],1);    
    
  %improve numerical stability
    Wstar([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) = ...
          1/2*(Wstar([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) + Wstar([Paras.XIndOffDiagCounter,Paras.XIndOffDiag]));
    
    PrimalAffine      = Paras.b-Paras.At*vec(Wstar);
    X_next            = omegat + (PrimalAffine)/Paras.alpha;
    PrimalFeasibility = norm(PrimalAffine)^2;
    gap               = abs(-Paras.b.'*X_next + Paras.c.'*reshape(Wstar,[],1));
end