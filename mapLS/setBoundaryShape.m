function map = setBoundaryShape(map)
    Ncell = map.connectivity.Ncell;
    LR = map.cells.vertices.LR;
    boundariesLR = map.cells.boundary.boundariesLR;
    controlvert = map.cells.vertices.LR;
    vertCoord = map.cells.vertices.coords;
    equations = map.cells.boundary.Eq;
    
    for i=1:Ncell
        cellvert = map.connectivity.cell{i};
        
        vrt2ch = find(map.cells.vertices.type(cellvert) == 'B');
        vrt2ch = cellvert(vrt2ch);

        nn = length(vrt2ch);
        for j=1:nn
            vLR = LR(vrt2ch(j),:);
            Bid = find(boundariesLR == (vLR(1)*10+vLR(2)));
            newVert = changeVert(vertCoord(vrt2ch(j),:),Bid);            
    
            map.connectivity.coord(vrt2ch(j),:) = newVert;
        end
        
    end

end



function vert = changeVert(vert,Bid)
    
    if (Bid==1 | Bid == 4)
        vert(2) = -sqrt(1-vert(1)^2);
    elseif (Bid==2 | Bid == 3)
        vert(2) = sqrt(1-vert(1)^2);
    end
end