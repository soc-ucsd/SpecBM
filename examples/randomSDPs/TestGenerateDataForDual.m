%Generate data
clc;clear;
addpath('..\..\packages\General');

for s = [10]
    fprintf("s = %d in action\n",s);
%     n  = 3000; 
%     m  = 100; 
    n = s*10;%dimension of X
    m = s*2;%number of affine constraints
    for dr = [997,995]%rank of Z
        %dr = 3;
        [At_sdp,b_sdp,c_sdp,X_sdp,y_sdp]= Generate_SDP_Problems_ForDual(n,m,n-dr);
        K_sdp.s = n;
        
%         res = SolveMosek(At_sdp,b_sdp,c_sdp,K_sdp); 
%         X = zeros(n);
%         [symidx]=SymmetricIndices(n,false);
%         X(symidx) = res.sol.itr.barx ;
%         X = X + tril(X,-1).';
%         Z  = c_sdp - At_sdp'*res.sol.itr.y;
%         TrX           = trace(X);
%         TrZ           = trace(mat(Z));
%         Optimal.TrX   = TrX;
%         Optimal.TrZ   = TrZ;
        Optimal.X_sdp = X_sdp;
        Optimal.y_sdp = y_sdp;
        Optimal.TrX   = trace(X_sdp);
        Optimal.TrZ   = trace(mat(c_sdp-At_sdp.'*y_sdp));
        Optimal.Cost  = c_sdp.'*vec(X_sdp);%res.sol.itr.pobjval;
        fprintf('trace(X) = %f\n',Optimal.TrX);
        fprintf('trace(Z) = %f\n',Optimal.TrZ);
        fprintf('Cost= %f\n', Optimal.Cost);
        %save(['n',int2str(n),'m',int2str(m),'dr',int2str(dr),'.mat'],'At_sdp','b_sdp','c_sdp','K_sdp','Optimal');
        %save(['n',int2str(n),'m',int2str(m),'dr',int2str(dr),'.mat'],'At_sdp','b_sdp','c_sdp','K_sdp','Optimal','-v7.3');
    end
end
