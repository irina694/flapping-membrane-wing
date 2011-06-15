function bezierC = cellSplinesCurve(xCell,yCell);
nP = length(xCell);
nB = nP -1;
if(nP>1)
    [xCellC,yCellC]=findCPoints(xCell,yCell);
    for indexB = 1:nB
        bezierC{indexB} = curve2([xCell(indexB),xCellC(2*indexB-1),xCellC(2*indexB),xCell(indexB+1)],...
            [yCell(indexB),yCellC(2*indexB-1),yCellC(2*indexB),yCell(indexB+1)]);
    end;
end;