function [IndSym,IndDiag,IndOffDiag,ShrinkIndDiag,ShrinkIndOffDiag,IndOffDiagCounter] = SymmetricIndices(n)
    %A is just a matrix variable
    A = ones(n);
    A = tril(A);
    IndSym = find(A);
    
    A = tril(A,-1);
    IndOffDiag = find(A);
    [row,col] = ind2sub([n,n],IndOffDiag);
    IndOffDiagCounter = sub2ind([n,n],col,row);
    
    
    A = diag(ones(1,n));
    IndDiag = find(A);
    
    ShrinkIndOffDiag = [];
    ShrinkIndDiag = [];
    
    k= 1;
    for i = 1:n
       for j=i:n
           if (j==i)
               ShrinkIndDiag = [ShrinkIndDiag,k];
           else
               ShrinkIndOffDiag = [ShrinkIndOffDiag,k];
           end
           k = k+1;
       end
    end

end