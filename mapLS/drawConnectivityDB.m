function drawConnectivityDB(connectivity)

clf;
hold on;

matrix = connectivity.matrix;
coord = connectivity.coord;
Xmin = min(coord(:,1));
Xmax = max(coord(:,1));
Ymin = min(coord(:,2));
Ymax = max(coord(:,2));

N=size(matrix,1);
for i=1:N
    for j=i+1:N
        if(matrix(i,j)==1)
            X=[coord(i,1) coord(j,1)];
            Y=[coord(i,2) coord(j,2)];
            plot(X,Y,'-bo','MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',2);       
            text(X(2),Y(2),num2str(j),'FontSize',16);
        end
    end
end
axis equal;
axis([Xmin Xmax Ymin Ymax]);

