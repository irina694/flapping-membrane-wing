function drawConnectivity3(connectivity)

clf;
hold on;

edgesvert = connectivity.edges.edgesvert;
coord = connectivity.verts.coords;
Ned = connectivity.edges.Nedges;

Xmin = min(coord(:,1));
Xmax = max(coord(:,1));
Ymin = min(coord(:,2));
Ymax = max(coord(:,2));


    
for i=1:Ned
    X = coord(edgesvert(i,:),1);
    Y = coord(edgesvert(i,:),2);
    plot(X,Y,'-bo','MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1);
end

axis equal;
axis([Xmin Xmax Ymin Ymax]);

