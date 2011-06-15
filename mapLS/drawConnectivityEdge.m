function drawConnectivityEdge(connectivity)

figure;
hold on;

coord = connectivity.vertices.coords;
edges = connectivity.edges;

Xmin = min(coord(:,1));
Xmax = max(coord(:,1));
Ymin = min(coord(:,2));
Ymax = max(coord(:,2));

wRef = 2;

for j=1:edges.Nedges
    jEdg = edges.vertices(j,:);
    X=[coord(jEdg(1),1) coord(jEdg(2),1)];
    Y=[coord(jEdg(1),2) coord(jEdg(2),2)];
    lineW = edges.thicknessRatio(j)*wRef;
    if edges.type(j) == 'B'
        colorL = 'k'; 
    else
        colorL = 'r';
    end
    plot(X,Y,['-',colorL],'LineWidth',lineW);
end

for j=1:connectivity.vertices.Nvertices
    X=[coord(j,1) coord(j,1)];
    Y=[coord(j,2) coord(j,2)];
    if connectivity.vertices.type(j) == 'B'
        colorV = 'k';
    elseif connectivity.vertices.type(j) == 'I'
        colorV = 'r';
    else
        colorV = 'b';
    end
    plot(X,Y,['-',colorV,'o'],'MarkerFaceColor',colorV);
end

axis equal;
axis([Xmin Xmax Ymin Ymax]);

