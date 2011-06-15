function drawMapCell(map)

    clf;
    hold on;

    coord = map.cells.edges.coord;
    width = map.cells.edges.thicknessRatio;
    type = map.cells.edges.type;
    Xmin = min(min(coord(:,[1 2])));
    Xmax = max(max(coord(:,[1 2])));
    Ymin = min(min(coord(:,[3 4])));
    Ymax = max(max(coord(:,[3 4])));

    epsx = 0.06;
    epsy = 0.03;
    
    Ncell = map.cells.Ncells;
    for j=1:Ncell
        Ned = length(map.cells.celledges{j});
        for i=1:Ned
            edg = map.cells.celledges{j}(i);            
            X = coord(edg,1:2);
            Y = coord(edg,3:4);
            if(type(edg)=='B')
                EdStyle={'-ob','b','b'};
            else
                EdStyle={'-or','r','r'};
            end
            %4*width(i)
            plot(X,Y,EdStyle{1},'MarkerEdgeColor',EdStyle{2},'MarkerFaceColor',EdStyle{3},'LineWidth',2);
            Xm=(X(1)+X(2))/2;
            Ym=(Y(1)+Y(2))/2;
            
            dir = map.cells.edgedir{j}(i);
            
            if(dir == '+')
                text(Xm-epsx,Ym+epsy,['+',num2str(edg)],'FontSize',12);
            else
                text(Xm+epsx,Ym+epsy,['-',num2str(edg)],'FontSize',12);
            end
        end
    end    
    axis equal;
    axis([Xmin Xmax Ymin Ymax]);

end

