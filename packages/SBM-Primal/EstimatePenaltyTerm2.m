function rho = EstimatePenaltyTerm(At_sdp,b_sdp,c_sdp,K_sdp,opts)
    %%%%%%%%%%deprecated%%%%%%%%%%%%%%
    %Any feasible primal and dual solution can provide an estimation on the trace(Z*) 

    idx    = find(eye(opts.n));
    Test   = At_sdp(:,idx);
    TraceA =  sum(Test,2);
    
    
    [x_1,y_1,info_1] = sedumi(At_sdp,b_sdp,c_sdp,K_sdp);
    sum(c_sdp(idx)) - TraceA'*y_1
    
    
    Atnew = At_sdp;
    %bnew  = ones(opts.m,1);
    bnew = -TraceA;
    cnew  = c_sdp;
    Knew  = K_sdp;
%     for i = 1:opts.m
%         if (b_sdp(i)~=0)
%             Atnew(i,:) = Atnew(i,:)/b_sdp(i);
%     end
    

    %We first inner and out approximate the optimal solution to have an
    %estimate on the penalty term
    
    [x_1,y_1,info_1] = sedumi(Atnew,bnew,cnew,Knew);
    dx = 1;
    [LowerBound,~] = ColumnGen_Both(Atnew,bnew,cnew,Knew,dx,1,true);
    [UpperBound,~] = ColumnGen_Both(Atnew,bnew,cnew,Knew,dx,1,false);
    
    
    rho = sum(c_sdp(idx)) - LowerBound;
%     lambad_C       = sort(eig(mat(cnew)),'descend');
%     lambda_A       = zeros(opts.n,opts.m);
% 
%     for i = 1 : opts.m
%         [~,eig_val]   = eig(mat(Atnew(i,:)));
%         [eig_val,~]   = sort(diag(eig_val),'descend');
%         lambda_A(:,i) = eig_val;
%     end
% 
%     by           = max([abs(LowerBound),abs(UpperBound)]);
%     max_lambda_A = max(lambda_A,[],2);
%     rho          = lambad_C+max_lambda_A*by;
%     rho          = sum(rho);
    
    
end