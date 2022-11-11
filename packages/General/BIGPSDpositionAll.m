function BigIndices = BIGPSDpositionAll(n,dx)
    %dx means a value in the partition () 
    CardOfP = n/dx;
   % NumOfVars = 2*dx*(2*dx+1)/2;
    NumOfVars = 4*dx^2;
    Blocks = reshape(1:n,dx,[])';
    Comb = nchoosek(1:CardOfP,2);
    [hei,wei] = size(Comb);
    NumOfComb = hei;
    Idices = zeros(NumOfComb,2*dx);
    for num = 1:NumOfComb
        i = Comb(num,1); j = Comb(num,2);
        Idices(num,:) = [Blocks(i,:),Blocks(j,:)];
    end

    BigIndices = zeros(NumOfComb,NumOfVars);
    for num = 1:NumOfComb
        BigIndices(num,:) = FindIndicesAll(Idices(num,:),n);
    end    
    %disp('wait');
end