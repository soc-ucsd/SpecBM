function [Wstar,X_next,Gammastar,Sstar,DualFeasibility,gap] = Direction_QP_Primal(omegat,Paras,Wt,Pt,feasible)
    %Author: Feng-Yi Liao
    %Wt is a fixed atoms
    %Pt is the transformation matrix
    %feasible means b-A(\Omega) = 0. 

    
    Q11 = Wt'*Wt;
    Q21 = kron(Pt',Pt')*Wt;
    Q31 = Paras.At_sdp*Wt;
    Q22 = kron(eye(Paras.MaxCols),eye(Paras.MaxCols));
    Q32 = Paras.At_sdp*kron(Pt,Pt);
    Q33 = Paras.At_sdp*Paras.At_sdp';
    %invQ33 = inv(Q33);
    %numerical stability
    %invQ33([Paras.mIndOffDiag,Paras.mIndOffDiagCounter]) = ...
    %      1/2*(invQ33(([Paras.mIndOffDiag,Paras.mIndOffDiagCounter])) + invQ33(([Paras.mIndOffDiagCounter,Paras.mIndOffDiag])));

    Q12 = Q21';
    Q13 = Q31';
    Q23 = Q32';
    
    q1  = 2*Wt'*(-Paras.c_sdp+Paras.alpha*omegat);
    q2  = 2*Paras.alpha*kron(Pt',Pt')*omegat-2*kron(Pt',Pt')*Paras.c_sdp;
    
    if feasible
        q3 = -2*(Paras.At_sdp*Paras.c_sdp);
    else
        q3 = -2*(Paras.alpha*(Paras.b_sdp-Paras.At_sdp*omegat)+Paras.At_sdp*Paras.c_sdp);
    end
    
    M11 = Q11 - Q13*(Q33\(Q13'));
    M22 = Q22 - Q23*(Q33\(Q23'));
    M12 = Q12 - Q13*(Q33\(Q23'));
    m1  = q1-Q13*(Q33\q3);
    m2  = q2-Q23*(Q33\q3);
    M   = [M11,M12;
          M12',M22];
    M   = (M+M.')/2;
    
    [eig_vec,eig_val]  = eig(M) ;
    %just to avoid numerical error (complex number)
    eig_val(eig_val<0) = 0;
    M05                = eig_vec*sqrt(eig_val)*eig_vec'; %B = M^{1/2}
    M05                = (M05+M05.')/2;
   
    %%%%%%%%%%%   VERY IMPORTANT!!!!!!!!!!!!   %%%%%%%%%%%%%%%%%%%%
    M05(:,1+[Paras.IndOffDiagPSD;Paras.IndOffDiagCounter]) = ...
      1/2*(M05(:,1+[Paras.IndOffDiagPSD;Paras.IndOffDiagCounter])+M05(:,1+[Paras.IndOffDiagCounter;Paras.IndOffDiagPSD]));
   
   %sedumi constraint matrix
   
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
    y                              = Q33\(-q3/2-Gammastar*Q13' - Q23'*reshape(Sstar,[],1));
    Wstar                          = Gammastar*Wt + reshape(Pt*Sstar*Pt',[],1);
    
    %improve numerical stability
    Wstar([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) = ...
          1/2*(Wstar([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) + Wstar([Paras.XIndOffDiagCounter,Paras.XIndOffDiag]));
    
    DualAffine      = Wstar-Paras.c_sdp+Paras.At_sdp'*y;
    X_next          = omegat + (DualAffine)/Paras.alpha;
    DualFeasibility = norm(DualAffine,'fro')^2;
    gap             = Paras.b_sdp'*y - Paras.c_sdp'*X_next;
end