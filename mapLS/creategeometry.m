function [g3D lengthB]=creategeometry(nTrees,delta0,rF0,dF0,rtF0,tF0,Y0,coef0,coef1,thkRoot,strTree)

global wingvert

lengthB = 0;

% Geometry
nvert=length(wingvert);
carr{1} = geomspline(wingvert');
carr{2}=curve2([wingvert(nvert,1),wingvert(1,1)],[wingvert(nvert,2),wingvert(1,2)]);
gM=geomcoerce('solid',carr);

gRoot = rect2(0,thkRoot,-wingvert(nvert,2)*0.5,wingvert(nvert,2)*1.1);

for k=1:nTrees
    gB{k} = treeTurtleGeoGeneral(strTree{k},0,Y0(k),0,rF0(k),dF0(k),rtF0(k),tF0(k),delta0(k),gRoot,gM);
end

gTrees=gB{1};
for k=2:nTrees
        gTrees = gB{k}+gTrees;
end

gTrees = geomdel(gTrees);
gTrees = fixSmallEdge(gTrees,1e-6);
gTrees = geomcomp({gTrees,gM},'ns',{'gM','gTrees'},'sf','gM*gTrees');

g2D = tesselateSD(gM,gTrees,coef0,coef1);

g3D=embed(g2D,'Wrkpln',[0 1 0;0 0 1;0 0 0]);