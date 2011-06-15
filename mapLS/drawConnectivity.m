function drawConnectivity(connectivity)

clf;
hold on;

edges = connectivity.edges;
vertices = connectivity.vertices;

Xmin = min(vertices.coords(:,1));
Xmax = max(vertices.coords(:,1));
Ymin = min(vertices.coords(:,2));
Ymax = max(vertices.coords(:,2));

N=edges.Nedges;
for j=1:N
    X=[vertices.coords(edges.vertices(j,1),1) vertices.coords(edges.vertices(j,2),1)];
    Y=[vertices.coords(edges.vertices(j,1),2) vertices.coords(edges.vertices(j,2),2)];
    plot(X,Y,'-bo','MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',2);
end

axis equal;
axis([Xmin Xmax Ymin Ymax]);

