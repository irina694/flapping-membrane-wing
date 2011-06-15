function geom = computeAreaAndNormal(connectivity,geom)
    Ncell = connectivity.Ncell;
    
    coord = connectivity.coord;
    overlap = [];
    for i=1:Ncell
        cellInd = connectivity.cell{i};
        aux=0;
        normal=[];
        nn = length(cellInd);
        Xsign=[];
        Ysign=[];
        for j=1:nn-1
            aux=aux+(coord(cellInd(j),1)-coord(cellInd(j+1),1))*...
                        (coord(cellInd(j),2)+coord(cellInd(j+1),2))/2;
            nvect=[(coord(cellInd(j+1),2)-coord(cellInd(j),2)),...
                       -(coord(cellInd(j+1),1)-coord(cellInd(j),1))];
            normal=[normal;nvect/norm(nvect)];
            Xsign=[Xsign,sign(coord(cellInd(j+1),1)-coord(cellInd(j),1))];
            Ysign=[Ysign,sign(coord(cellInd(j+1),2)-coord(cellInd(j),2))];

%plot([coord(cellInd(j),1) coord(cellInd(j+1),1)],[coord(cellInd(j),2) coord(cellInd(j+1),2)]);

        end
        aux=aux+(coord(cellInd(nn),1)-coord(cellInd(1),1))*...
                    (coord(cellInd(nn),2)+coord(cellInd(1),2))/2;
        area(i)=abs(aux);
        nvect=[(coord(cellInd(1),2)-coord(cellInd(nn),2)),...
                   -(coord(cellInd(1),1)-coord(cellInd(nn),1))];
        normal=[normal;nvect/norm(nvect)];
        Xsign=[Xsign,sign(coord(cellInd(1),1)-coord(cellInd(nn),1))];
        Ysign=[Ysign,sign(coord(cellInd(1),2)-coord(cellInd(nn),2))];
        
        geom.normal{i}=normal;
        overlap=[overlap;checkOverlapping(Xsign,Ysign)];
    end
    
    geom.area=area;
    geom.overlap=overlap;
end



function value = checkOverlapping(X,Y)

    nn = length(X);
    X0 = find(X==0);
    Y0 = find(Y==0);
    
    nx = length(X0);
    for i=1:nx
        if(X0(i)==1)
           X(1)=X(2); 
        else
           X(X0(i))=X(X0(i)-1);
        end
    end
    
    ny = length(Y0);
    for i=1:ny
        if(Y0(i)==1)
           Y(1)=Y(2); 
        else
           Y(Y0(i))=Y(Y0(i)-1);
        end
    end

    Xswitch = 0;
    Xsign = X(1);
    for i=2:nn
        if(X(i)~=Xsign)
            Xswitch = Xswitch+1;
            Xsign = X(i);
        end
    end
    

    Yswitch = 0;
    Ysign = Y(1);
    for i=2:nn
        if(Y(i)~=Ysign)
            Yswitch = Yswitch+1;
            Ysign = Y(i);
        end
    end
    
    if(Xswitch>2 || Yswitch>2)
        value=1;
    else
        value=0;
    end

end