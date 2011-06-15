function fixedGeo = delSmallEdge(geom,lMin)

geomE = geomedit(geom);
nEd = flgeomnbs(geom);

for i=1:nEd
    vtx_x = get(geomE{i},'ctrlx');
    vtx_y = get(geomE{i},'ctrly');
    vt_x1=vtx_x(1);
    vt_y1=vtx_y(1);
    vt_x2=vtx_x(end);
    vt_y2=vtx_y(end);
    lengthE(i)=sqrt((vt_x2-vt_x1)^2+(vt_y2-vt_y1)^2);
end
smallEdges = find(lengthE<lMin);
if(length(smallEdges)==0)
    fixedGeo = geom;
else
    fixedGeo = geomcomp(geom,'edge',smallEdges);
end;