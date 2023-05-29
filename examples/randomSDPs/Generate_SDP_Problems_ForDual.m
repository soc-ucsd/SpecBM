%Generate Random SDPs
%We have restriction on the rank of the primal variable
function [At,b,c,X,y]= Generate_SDP_Problems_ForDual(n,m,r_p)
    %n: n by n symmetric matrices
    %m: number of linear constraints
    rng('default');
    Struc= ones(n);
    %fixed X
    scaling = 1;

    X = scaling*sprandsym(Struc);

    eigval = eig(X);
    if (min(eigval)<0)%%make X feasible
        %test = (min(eigval)-1)*speye(size(X));
        X = X - (min(eigval)-1)*speye(size(X));
    end
            
    [eig_vec,eig_val] = eig(full(X));
    support_x = randsample(n,r_p);
%    indices = 1:n;
%     idx_x = indices(support_x);
    
    temp            = zeros(n,1);
    temp(support_x) = 1;
    support_z       = find(~temp);
   
    P1 = eig_vec(:,support_x);
    P2 = eig_vec(:,support_z);
   
    D1 = eig_val(support_x,support_x);
    D2 = eig_val(support_z,support_z);
    
    X = P1*D1*P1';
    Z = P2*D2*P2';
    
    X = (X+X')/2;%numerical stability
    Z = (Z+Z')/2;%numerical stability 
    
    A = cell(m,1);
    %At = [];
    At = zeros(m,n^2);
    idx = find(eye(n));
    for j=1:m
        A{j}      = scaling*sprandsym(Struc);
        k         = trace(A{j}); 
        A{j}(idx) = A{j}(idx)-k/n; %To make sure Trace(A_i) = 0. With this we can make sure there exist a strictly feasible primal solution
        %At        = [At,vec(A{j})];
        At(j,:)   = vec(A{j})';
    end
    %At = full(At);
    b  = At*vec(X);
    y  = rand(m,1);
    
    C = Z;
    for j = 1:m
        C = C - A{j}*y(j);
    end
    %C = C/trace(C);
    c = reshape(C,[],1);
    %At = At';
     
%     fprintf('trace(X) = %f\n',trace(X));
%     fprintf('trace(Z) = %f\n',trace(Z));
%     fprintf('Cost= %f\n', abs(trace(C*X)));
%     fprintf('KKT = %f\n',trace(X*Z));
end