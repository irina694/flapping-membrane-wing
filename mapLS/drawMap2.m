function drawMap2(map,options)

    if(isempty(options))
        lables=false;
    else
        lables=true;        
    end
        

    clf;
    hold on;        

    coord = map.cells.edges.coord;
    type = map.cells.edges.type;
    Xmin = min(min(coord(:,[1 2])));
    Xmax = max(max(coord(:,[1 2])));
    Ymin = min(min(coord(:,[3 4])));
    Ymax = max(max(coord(:,[3 4])));

    Ned=map.cells.edges.Nedges;
    for i=1:Ned
        X=coord(i,1:2);
        Y=coord(i,3:4);
        if(type(i)=='B')
            EdStyle={'-ob','b','b'};
        else
            EdStyle={'-or','r','r'};
        end
        
        plot(X,Y,EdStyle{1},'MarkerEdgeColor',EdStyle{2},'MarkerFaceColor',EdStyle{3},'LineWidth',2);
        
        if(lables)      
            Xm=(X(1)+X(2))/2;        
            Ym=(Y(1)+Y(2))/2;
            text(Xm,Ym,num2str(i),'FontSize',12);
        end
    end
    
    axis equal;
    axis([Xmin Xmax Ymin Ymax]);

end

