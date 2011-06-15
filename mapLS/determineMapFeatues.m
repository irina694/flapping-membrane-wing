function features = determineMapFeatues(map)

    features.Ncells = map.cells.Ncells;
    areas = computeAreas(map.cells);
    features.areas = areas;
    features.MinArea = min(areas);
    features.MaxArea = max(areas);
    features.StdArea = std(areas);
    
    perimeter = computePerimeter(map.cells);
    AR = ((perimeter/4+sqrt(perimeter.^2/16-areas)).^2)./areas;
    cellInternalAngles = computeCellInternalAngles(map.cells);
    

    
    EdgeCoord = map.cells.edges.coord;
    InternalEdges = find(map.cells.edges.type == 'I');
    features.NintEdges = length(InternalEdges);

    coordX = EdgeCoord(InternalEdges,[1,2])';
    coordY = EdgeCoord(InternalEdges,[3,4])';
    features.EdgeLength = sum(sqrt(diff(coordX).^2 + diff(coordY).^2));
    
    edgeAngles = abs(computesEdgesAngles(coordX,coordY))*180/pi;
    feature.edgeAngles = edgeAngles;
    features.MinAngle = min(edgeAngles);
    features.MaxAngle = max(edgeAngles);
    features.StdAngle = std(edgeAngles);    




end




function cellInternalAngles = computeCellInternalAngles(cells)
    Ncell = cells.Ncells;
    
    perimeter = zeros(1,Ncell);
    coord = cells.vertices.coords;    
    for i=1:Ncell
        cellInd = cells.cellvert{i};
        aux=0;
        nn = length(cellInd);
        for j=1:nn-1
            aux=aux+sqrt((coord(cellInd(j),1)-coord(cellInd(j+1),1))^2+...
                        (coord(cellInd(j),2)-coord(cellInd(j+1),2))^2);
        end
        perimeter(i)=aux+sqrt((coord(cellInd(nn),1)-coord(cellInd(1),1))^2+...
                        (coord(cellInd(nn),2)-coord(cellInd(1),2))^2);
    end

end


function edgeAngles = computesEdgesAngles(coordX,coordY)

    N = length(coordX);
    edgeAngles = zeros(1,N);
    
    for i=1:N        
        edgeAngles(i) = atan((coordY(2,i)-coordY(1,i))/(coordX(2,i)-coordX(1,i)));
    end
  
end


function perimeter = computePerimeter(cells)
    Ncell = cells.Ncells;
    
    perimeter = zeros(1,Ncell);
    coord = cells.vertices.coords;    
    for i=1:Ncell
        cellInd = cells.cellvert{i};
        aux=0;
        nn = length(cellInd);
        for j=1:nn-1
            aux=aux+sqrt((coord(cellInd(j),1)-coord(cellInd(j+1),1))^2+...
                        (coord(cellInd(j),2)-coord(cellInd(j+1),2))^2);
        end
        perimeter(i)=aux+sqrt((coord(cellInd(nn),1)-coord(cellInd(1),1))^2+...
                        (coord(cellInd(nn),2)-coord(cellInd(1),2))^2);
    end
    
end


function areas = computeAreas(cells)
    Ncell = cells.Ncells;
    
    areas = zeros(1,Ncell);
    coord = cells.vertices.coords;    
    for i=1:Ncell
        cellInd = cells.cellvert{i};
        aux=0;
        nn = length(cellInd);
        for j=1:nn-1
            aux=aux+(coord(cellInd(j),1)-coord(cellInd(j+1),1))*...
                        (coord(cellInd(j),2)+coord(cellInd(j+1),2))/2;
        end
        aux=aux+(coord(cellInd(nn),1)-coord(cellInd(1),1))*...
                    (coord(cellInd(nn),2)+coord(cellInd(1),2))/2;
        areas(i)=abs(aux);
    end
    
end
