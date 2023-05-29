function Out = KronShrink(Pt,n,r,rowidx)
    %Authors: Feng-Yi Liao & Yang Zheng
    %         SOC Lab @UC San Diego
    %Pt is the transformation matrix for the subproblem
    
    
    %[n,r]    = size(Pt);
    %[rowidx] = RowSymmetric(n);
    reduce_n = n*(n+1)/2;
    Out      = zeros(r^2,reduce_n);
    Ptt      = Pt'; %Transpose of Pt
    start    = 1;
        
    %First extract the symmetric part
    for i = 1:n
        indices = i:n;
        len     = n-i+1;
        %Out(:,start:start+len-1) =  Pt(:,i)*multiPtt(:,indices); 
        Out(:,start:start+len-1) = kron(Ptt(:,i),Ptt(:,indices));
        start   = start + len;
    end
    
    %Add the corresponding part
    start  = 1;
    for i = 1:n-1
        indices = 1:i;
        Out(:,rowidx(start:start+i-1)) = Out(:,rowidx(start:start+i-1)) + kron(Ptt(:,i+1),Ptt(:,indices));
        start = start + i;
    end
   
end