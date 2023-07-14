function [rowidx] =  RowSymmetric(n)
    rowidx = zeros(1,n*(n-1)/2);
    
    M = zeros(n);
    start = 1;
    len = n;
    for i=1:n
       M(i:n,i) = start:start+len-1;
       start = start + len;
       len = len -1 ;
    end
   
    start = 1;
    for i = 2:n
        rowidx(start:start+i-1-1) = M(i,1:i-1);
        start = start+ i-1;
    end
end