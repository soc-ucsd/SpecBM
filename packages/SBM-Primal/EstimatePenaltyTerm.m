function rho = EstimatePenaltyTerm(At_sdp,b_sdp,c_sdp,K_sdp,opts)
    %%%%%%%%%%deprecated%%%%%%%%%%%%%%
    %Any feasible primal and dual solution can provide an estimation on the trace(Z*) 
%     Atnew = Paras.At_sdp;
%     bnew  = ones(Paras.m,1);
%     cnew  = Paras.c_sdp;
%     Knew  = Paras.K_sdp;
%     for i = 1:Paras.m
%         Atnew(i,:) = Atnew(i,:)/Paras.b_sdp(i);
%     end
    Atnew = At_sdp;
    bnew  = ones(opts.m,1);
    cnew  = c_sdp;
    Knew  = K_sdp;
    for i = 1:opts.m
        if (b_sdp(i)~=0)
            Atnew(i,:) = Atnew(i,:)/b_sdp(i);
        end
    end
    
    %We first inner and out approximate the optimal solution to have an
    %estimate on the penalty term
    dx = 5;
    [LowerBound,~] = ColumnGen_Both(Atnew,bnew,cnew,Knew,dx,1,true);
    [UpperBound,~] = ColumnGen_Both(Atnew,bnew,cnew,Knew,dx,1,false);
    
    lambad_C       = sort(eig(mat(cnew)),'descend');
    lambda_A       = zeros(opts.n,opts.m);

    for i = 1 : opts.m
        [~,eig_val]   = eig(mat(Atnew(i,:)));
        [eig_val,~]   = sort(diag(eig_val),'descend');
        lambda_A(:,i) = eig_val;
    end

    by           = max([abs(LowerBound),abs(UpperBound)]);
    max_lambda_A = max(lambda_A,[],2);
    rho          = lambad_C+max_lambda_A*by;
    rho          = sum(rho);
    
    
end