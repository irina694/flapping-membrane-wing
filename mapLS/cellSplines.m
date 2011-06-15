function bezierC = cellSplines(xCell,yCell,tCellIn,fCap);
nP = length(xCell);
sPar = 1:nP;
% apply factor to thickness
tCell=tCellIn*fCap;
% splines
centralCurveX = spline(sPar,xCell);
centralCurveY = spline(sPar,yCell);
% splines derivative
DcentralCurveX = fnder(centralCurveX);
DcentralCurveY = fnder(centralCurveY);
% calculate the tangent vector at each node
tangentX = fnval(DcentralCurveX,sPar);
tangentY = fnval(DcentralCurveY,sPar);
tangentM = sqrt(tangentX.^2+tangentY.^2);
tangentX = tangentX./tangentM;
tangentY = tangentY./tangentM;
% calculate the unit normal vector at each node
normalX = -tangentY;
normalY = tangentX;
% determine upper curve
tCellHalf = tCell/2;
xCellUp = xCell+tCellHalf.*normalX;
yCellUp = yCell+tCellHalf.*normalY;
bezierUp = geomspline([xCellUp;yCellUp]);
% determine lower curve
xCellDown = xCell-tCellHalf.*normalX;
yCellDown = yCell-tCellHalf.*normalY;
bezierDown = geomspline([xCellDown;yCellDown]);
% determine right cap
rCap=tCellHalf(end);
xRightCap11 = xCellUp(end)+rCap*tangentX(end);
xRightCap12 = xCell(end)+rCap*tangentX(end);
yRightCap11 = yCellUp(end)+rCap*tangentY(end);
yRightCap12 = yCell(end)+rCap*tangentY(end);
rightCap1 = curve2([xCellUp(end) xRightCap11 xRightCap12],...
    [yCellUp(end) yRightCap11 yRightCap12],[1 cos(pi/4) 1]);
xRightCap11 = xCell(end)+rCap*tangentX(end);
xRightCap12 = xCellDown(end)+rCap*tangentX(end);
yRightCap11 = yCell(end)+rCap*tangentY(end);
yRightCap12 = yCellDown(end)+rCap*tangentY(end);
rightCap2 = curve2([xRightCap11 xRightCap12 xCellDown(end)],...
    [yRightCap11 yRightCap12 yCellDown(end)],[1 cos(pi/4) 1]);
% determine left cap
rCap=tCellHalf(1);
xLeftCap11 = xCellUp(1)-rCap*tangentX(1);
xLeftCap12 = xCell(1)-rCap*tangentX(1);
yLeftCap11 = yCellUp(1)-rCap*tangentY(1);
yLeftCap12 = yCell(1)-rCap*tangentY(1);
leftCap1 = curve2([xCellUp(1) xLeftCap11 xLeftCap12],...
    [yCellUp(1) yLeftCap11 yLeftCap12],[1 cos(pi/4) 1]);
xLeftCap11 = xCell(1)-rCap*tangentX(1);
xLeftCap12 = xCellDown(1)-rCap*tangentX(1);
yLeftCap11 = yCell(1)-rCap*tangentY(1);
yLeftCap12 = yCellDown(1)-rCap*tangentY(1);
leftCap2 = curve2([xLeftCap11 xLeftCap12 xCellDown(1)],...
    [yLeftCap11 yLeftCap12 yCellDown(1)],[1 cos(pi/4) 1]);

c.objs = {bezierUp,rightCap1,rightCap2,...
    bezierDown,leftCap1,leftCap2};
fem.draw = struct('c',c);
bezierC = geomcsg(fem,'solidify','on');