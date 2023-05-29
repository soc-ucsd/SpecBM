A = mat(1:16);
At = A + A.';
At = reshape(At,[],1);

C = mat(101:116);
C = C + C.';
C = reshape(C,[],1);

K.s = 4;

K.f = 0;
K.l = 0;
K.q = 0;

[Ats,Cs,Q] = svecData(At,C,K);
Q'*Ats