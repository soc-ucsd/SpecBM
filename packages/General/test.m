%test
n = 4;
rng('default');
%Pt = rand(n);
Pt = reshape(1:n^2,n,n);
KronShrink(Pt);
