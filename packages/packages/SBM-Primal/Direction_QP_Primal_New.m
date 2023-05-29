function [Wstar,X_next,Gammastar,Sstar,DualFeasibility,gap] = Direction_QP_Primal_New(omegat,Paras,Wt,Pt,feasible)
    %%%%%%%%%%deprecated%%%%%%%%%%%%%%
    %Authors: Feng-Yi Liao & Yang Zheng
    %         SOC Lab @UC San Diego
    %Wt is a fixed atoms
    %Pt is the transformation matrix
    %feasible means b-A(\Omega) = 0. 

    
%    kronPtPt = kron(Pt,Pt);
    Q11 = Wt'*Wt;
    %Q21 = kron(Pt',Pt')*Wt;
    %Q21 = kronPtPt'*Wt;
    Q21 = vec(Pt'*mat(Wt)*Pt);
    
    %traditional way
    %Q31 = Paras.At_sdp*Wt;
    
    Q31 = zeros(Paras.m,1);
    for i = 1:length(Paras.At_sdp_nonzero.row)
        idx      = Paras.At_sdp_nonzero.idx(i);
        row      = Paras.At_sdp_nonzero.row(i);
        val      = Paras.At_sdp_nonzero.val(i);
        Q31(row) = Q31(row)+val*Wt(idx);
    end
    
    
    
    
    %Paras.Ir2 = kron(eye(Paras.MaxCols),eye(Paras.MaxCols));
    %Q32 = Paras.At_sdp*kron(Pt,Pt);
    
    %Original strategy
    %Q32_ori = Paras.At_sdp*kronPtPt;
    
    
%     %Shrink version
%     %At_sdp_shrink = Paras.At_sdp(:,Paras.XIndSym);
%     kronPtPt_temp = KronShrink(Pt,Paras.n,Paras.MaxCols,Paras.rowidx);
%     % = KronShrink(Pt',Paras.n,Paras.MaxCols,Paras.rowidx);
%     Q32 = Paras.At_sdp_shrink*kronPtPt_temp';
     
    
    
    Q32 = zeros(Paras.m,Paras.MaxCols^2);
    for i = 1:length(Paras.At_sdp_nonzero.row) 
        idx = Paras.At_sdp_nonzero.idx(i);
        row  = Paras.At_sdp_nonzero.row(i);
        Mrow = Paras.At_sdp_nonzero.Mrow(i);
        Mcol = Paras.At_sdp_nonzero.Mcol(i);
        Q32(row,:)  = Q32(row,:) + Paras.At_sdp(row,idx)*kron(Pt(Mcol,:),Pt(Mrow,:));
        %Q32(row,:) = Q32(row,:) + Paras.At_sdp(row,col)*kronPtPt(col,:);
    end
    
    
    %Q33 = Paras.At_sdp*Paras.At_sdp'; AAT
    %invQ33 = inv(Q33);
    %numerical stability
    %invQ33([Paras.mIndOffDiag,Paras.mIndOffDiagCounter]) = ...
    %      1/2*(invQ33(([Paras.mIndOffDiag,Paras.mIndOffDiagCounter])) + invQ33(([Paras.mIndOffDiagCounter,Paras.mIndOffDiag])));

    Q12 = Q21';
    Q13 = Q31';
    Q23 = Q32';
    
    alpha_omegat = Paras.alpha*omegat;
    q1  = 2*Wt'*(-Paras.c_sdp+alpha_omegat);
    %q2  = 2*Paras.alpha*kron(Pt',Pt')*omegat-2*kron(Pt',Pt')*Paras.c_sdp;
    %q2  = 2*Paras.alpha*kronPtPt'*omegat-2*kronPtPt'*Paras.c_sdp;
    
%     temp = alpha_omegat;
%     for i = 1:length(Paras.c_sdp_nonzero.idx) 
%         posi = Paras.c_sdp_nonzero.idx(i);
%         temp(posi) = temp(posi) - Paras.c_sdp(posi);
%     end
%    q2  = kronPtPt'*(alpha_omegat-Paras.c_sdp)*2;
    q2  = vec(Pt'*mat((alpha_omegat-Paras.c_sdp)*2)*Pt);
    %q2  = kronPtPt'*temp*2;
    
    if feasible
        %q3 = -2*(Paras.At_sdp*Paras.c_sdp);
        q3 = Paras.q3;
    else
        %q3 = -2*(Paras.alpha*(Paras.b_sdp-Paras.At_sdp*omegat)+Paras.At_sdp*Paras.c_sdp);
        q3 = -2*(Paras.alpha*(Paras.b_sdp-Paras.At_sdp*omegat))+Paras.q3;
    end
    
%     M11 = Q11 - Q13*(Paras.AAT\(Q13'));
%     M22 = Paras.Ir2 - Q23*(Paras.AAT\(Q23'));
%     M12 = Q12 - Q13*(Paras.AAT\(Q23'));
%     m1  = q1-Q13*(Paras.AAT\q3);
%     m2  = q2-Q23*(Paras.AAT\q3);
  
    M11 = Q11 - Q13*(Paras.AAT_INV*(Q13.'));
    M22 = Paras.Ir2 - Q23*(Paras.AAT_INV*(Q23'));
    M12 = Q12 - Q13*(Paras.AAT_INV*(Q23'));
    m1  = q1-Q13*(Paras.AAT_INV*q3);
    m2  = q2-Q23*(Paras.AAT_INV*q3);
    
    
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
    %y                              = Paras.AAT\(-q3/2-Gammastar*Q13' - Q23'*reshape(Sstar,[],1));
    y                              = Paras.AAT_INV*(-q3/2-Gammastar*Q13' - Q23'*reshape(Sstar,[],1));
    Wstar                          = Gammastar*Wt + reshape(Pt*Sstar*Pt',[],1);
    
    %improve numerical stability
    Wstar([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) = ...
          1/2*(Wstar([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) + Wstar([Paras.XIndOffDiagCounter,Paras.XIndOffDiag]));
    
    DualAffine      = Wstar-Paras.c_sdp+Paras.At_sdp'*y;
%     DualAffine       = Wstar - Paras.c_sdp;
%     for i = 1:length(Paras.At_sdp_nonzero.row) 
%         idx              = Paras.At_sdp_nonzero.idx(i);
%         val              = Paras.At_sdp_nonzero.val(i);
%         DualAffine(idx)  = DualAffine(idx) + val*y(row);
%     end
    
    
    X_next          = omegat + (DualAffine)/Paras.alpha;
    DualFeasibility = norm(DualAffine,'fro')^2;
    gap             = Paras.b_sdp'*y - Paras.c_sdp'*X_next;
end