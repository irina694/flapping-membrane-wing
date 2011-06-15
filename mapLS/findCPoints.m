function  [controlPointsX,controlPointsY] = findCPoints(Px,Py)

nP = length(Px);
nB = nP - 1;

iM(1) = 1;
jM(1) = 1;
sM(1) = -2;

iM(2) = 1;
jM(2) = 2;
sM(2) = 1;

indexM= 2;
for i = 1:nB-1
    j=i-1;
    indexM = indexM + 1;
    iM(indexM) = 2*i;
    jM(indexM) = 2*j+1;
    sM(indexM) = 1;
    indexM = indexM+1;
    iM(indexM) = 2*i;
    jM(indexM) = 2*j+2;
    sM(indexM) = -2;
    indexM = indexM+1;
    iM(indexM) = 2*i;
    jM(indexM) = 2*j+3;
    sM(indexM) = 2;
    indexM = indexM+1;
    iM(indexM) = 2*i;
    jM(indexM) = 2*j+4;
    sM(indexM) = -1;
    indexM = indexM+1;
    iM(indexM) = 2*i+1;
    jM(indexM) = 2*j+2;
    sM(indexM) = 1;
    indexM = indexM+1;
    iM(indexM) = 2*i+1;
    jM(indexM) = 2*j+3;
    sM(indexM) = 1;
end;
indexM = indexM+1;
iM(indexM) = 2*nB;
jM(indexM) = 2*nB-1;
sM(indexM) = 1;
indexM = indexM+1;
iM(indexM) = 2*nB;
jM(indexM) = 2*nB;
sM(indexM) = -2;

AM = sparse(iM,jM,sM);

[LM UM PM QM] = lu(AM);

bM(1) = -Px(1);
for i=1:nB-1
    bM(2*i) = 0;
    bM(2*i+1) = 2*Px(i+1);
end;
bM(2*nB) = -Px(nP);

controlPointsX = AM\bM';

bM(1) = -Py(1);
for i=1:nB-1
    bM(2*i) = 0;
    bM(2*i+1) = 2*Py(i+1);
end;
bM(2*nB) = -Py(nP);

controlPointsY = AM\bM';
    