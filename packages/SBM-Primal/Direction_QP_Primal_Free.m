function [Wstar,X_next,Gammastar,Sstar,DualFeasibility,DualFeasibility_free,gap,y,x_next] = Direction_QP_Primal_Free(omegat_free,omegat_sdp,Paras,Wt,Pt,feasible)
    %Authors: Feng-Yi Liao & Yang Zheng
    %         SOC Lab @UC San Diego
    %Wt is a fixed atoms
    %Pt is the transformation matrix
    %feasible means b-A(\Omega) = 0. 
    %Consider free variables as well
    %The two changes are Q_{33} and q_{3}
    %But Q_{33} is pre-computed
    
    kronPtPt = kron(Pt,Pt);
    
    Q11 = Wt.'*Wt;
    %Q21 = kron(Pt',Pt')*Wt;
    %Q21 = kronPtPt'*Wt;
    Q12 = Wt.'*kronPtPt;
    Q31 = Paras.At_sdp*Wt;
    %Paras.Ir2 = kron(eye(Paras.MaxCols),eye(Paras.MaxCols));
    %Q32 = Paras.At_sdp*kron(Pt,Pt);
    
    %Original strategy
    Q32 = Paras.At_sdp*kronPtPt;
    
    Q13 = Q31.';
    Q23 = Q32.';
    
    temp = 2*(-Paras.c_sdp+Paras.alpha*omegat_sdp);
    q1   = Wt.'*temp;
    q2T  = temp.'*kronPtPt;
    q2  = q2T.';

%     q1   = 2*Wt.'*(-Paras.c_sdp+Paras.alpha*omegat_sdp);
%     %q2  = 2*Paras.alpha*kron(Pt',Pt')*omegat-2*kron(Pt',Pt')*Paras.c_sdp;
%     %q2  = 2*Paras.alpha*kronPtPt'*omegat-2*kronPtPt'*Paras.c_sdp;
%     %q2  = 2*kronPtPt'*(Paras.alpha*omegat-Paras.c_sdp);
%     q2T = 2*(Paras.alpha*omegat_sdp-Paras.c_sdp).'*kronPtPt;
%     q2  = q2T.';

    if feasible
        %q3 = -2*(Paras.At_sdp*Paras.c_sdp+Paras.a*Paras.c0);
        q3 = -2*(Paras.At*Paras.c);
    else
        %q3 = -2*(Paras.alpha*(Paras.b_sdp-Paras.At_sdp*omegat_sdp-Paras.a*omegat_free)+Paras.At_sdp*Paras.c_sdp+Paras.a*Paras.c0);
        q3 = -2*(Paras.alpha*(Paras.b-Paras.At_sdp*omegat_sdp-Paras.At_free*omegat_free)+Paras.At*Paras.c);
    end
    
%     M11 = Q11 - Q13*(Paras.AAT\(Q13'));
%     M22 = Paras.Ir2 - Q23*(Paras.AAT\(Q23'));
%     M12 = Q12 - Q13*(Paras.AAT\(Q23'));
%     m1  = q1-Q13*(Paras.AAT\q3);
%     m2  = q2-Q23*(Paras.AAT\q3);
    
    M11 = Q11 - Q13*(Paras.AAT_INV*(Q31));
    %time comsuming
    M22 = Paras.Ir2 - Q23*Paras.AAT_INV*Q32;
    M12 = Q12 - Q13*Paras.AAT_INV*Q32;
    m1  = q1-Q13*(Paras.AAT_INV*q3);
    m2  = q2-Q23*(Paras.AAT_INV*q3);   
    
    M   = [M11,M12;
          M12.',M22];
   % M   = (M+M.')/2;
     
    if Paras.EvecPast == 0 && Paras.EvecCurrent == 1 && feasible
         %In this case, we have closed form solution, so we don't need to
         %rely on solvers!!!
         m         = [m1;m2];  
         
         v = -inv(M)*m/2;
        
         if v(1)< 0 || v(2) < 0  || v(1) + v(2) > Paras.rho
         
             candidate = cell(3,1);

             denominator = (2*M11+2*M22-4*M12);
             if denominator ~= 0 
                Guess1 = (2*M22*Paras.rho-2*M12*Paras.rho-m1+m2)/denominator;
             else
                Guess1 =0 ;
             end

             candidate{1} = {};
             if Guess1 >= 0 && Guess1 <= Paras.rho
               candidate{1}.v1 = Guess1;
             elseif Guess1 < 0
               candidate{1}.v1 = 0;
             elseif Guess1 > Paras.rho
               candidate{1}.v1 = Paras.rho;
             end
             candidate{1}.v2 = Paras.rho-candidate{1}.v1;
             candidate{1}.v  = [candidate{1}.v1;candidate{1}.v2];


             Guess2 = -m2/(2*M22);
             candidate{2} = {};
             if Guess2 >= 0 && Guess2 <= Paras.rho
               candidate{2}.v2 = Guess2;
             elseif Guess2 < 0
               candidate{2}.v2 = 0;
             elseif Guess2 > Paras.rho
               candidate{2}.v2 = Paras.rho;
             end        
             candidate{2}.v1 = 0;
             candidate{2}.v  = [candidate{2}.v1;candidate{2}.v2];

             Guess3 = -m1/(2*M11);
             candidate{3} = {};
             if Guess3 >= 0 && Guess3 <= Paras.rho
               candidate{3}.v1 = Guess3;
             elseif Guess3 < 0
               candidate{3}.v1 = 0;
             elseif Guess3 > Paras.rho
               candidate{3}.v1 = Paras.rho;
             end               
             candidate{3}.v2 = 0;
             candidate{3}.v  = [candidate{3}.v1;candidate{3}.v2];

             f1        = candidate{1}.v.'*M*candidate{1}.v + m.'*candidate{1}.v;
             f2        = candidate{2}.v.'*M*candidate{2}.v + m.'*candidate{2}.v;
             f3        = candidate{3}.v.'*M*candidate{3}.v + m.'*candidate{3}.v;
             [obj,idx]   = min([f1,f2,f3]);
             Gammastar = candidate{idx}.v1;
             Sstar     = candidate{idx}.v2;
             %Wstar     = Gammastar*Wt + reshape(Pt*Sstar*Pt.',[],1);
         else
         
            Gammastar = v(1);
            Sstar     = v(2);
            
         end
         
%          Wstar     = Gammastar*Wt + reshape(Pt*Sstar*Pt.',[],1);
%          y         = Paras.AAT_INV*(-q3/2-Gammastar*Q13.' - Q23'*reshape(Sstar,[],1));

    else
        [eig_vec,eig_val]  = eig(M) ;
        %just to avoid numerical error (complex number)
        eig_val(eig_val<0) = 0;
        M05                = eig_vec*sqrt(eig_val)*eig_vec.'; %B = M^{1/2}
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
        %y                              = Paras.AAT\(-q3/2-Gammastar*Q13' - Q23'*reshape(Sstar,[],1));
%         y                              = Paras.AAT_INV*(-q3/2-Gammastar*Q13.' - Q23'*reshape(Sstar,[],1));
% 
%         Wstar                          = Gammastar*Wt + reshape(Pt*Sstar*Pt.',[],1);
    end 
    
    %numerical issue
    if Gammastar< 0
        Gammastar = 0;
    end

%     issymmetric(Sstar)
%     issymmetric(Pt*Sstar*Pt.')
    Wstar     = Gammastar*Wt + reshape(Pt*Sstar*Pt.',[],1);
    y         = Paras.AAT_INV*(-q3/2-Gammastar*Q13.' - Q23'*reshape(Sstar,[],1));
    
     
    DualAffine_sdp      = Wstar-Paras.c_sdp+Paras.At_sdp.'*y;
    DualAffine_free     = Paras.At_free.'*y-Paras.c_free;
    
    %improve numerical stability %very important!!!
    DualAffine_sdp([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) = ...
          1/2*(DualAffine_sdp([Paras.XIndOffDiag,Paras.XIndOffDiagCounter]) + DualAffine_sdp([Paras.XIndOffDiagCounter,Paras.XIndOffDiag]));

    X_next          = omegat_sdp + (DualAffine_sdp)/Paras.alpha;
    x_next          = omegat_free + (DualAffine_free)/Paras.alpha;
    DualFeasibility = norm(DualAffine_sdp,'fro')^2;
    DualFeasibility_free = norm(DualAffine_free)^2;
    gap             = abs(Paras.b'*y - Paras.c_sdp.'*X_next-Paras.c_free.'*x_next);
end