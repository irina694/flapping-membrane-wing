function  map = updateMapstruct(map)

    C = map.connectivity;
    map.cells.vertices.coords = C.coord;
    Ned = map.connectivity.edges.Nedges;
    coord = zeros(Ned,4);
    
    for i=1:Ned
        v = map.connectivity.edges.edgesvert(i,:);
        coord(i,:) = [C.coord(v(1),1),C.coord(v(2),1),C.coord(v(1),2),C.coord(v(2),2)];
    end
    map.cells.edges.coord = coord;
    map.connectivity.edges.coord = coord;
end



function  map = updateMapstruct_old(map)

    C = map.connectivity;
    map.cells.vertices.coords = C.coord;
    Ned = map.cells.edges.Nedges;
    coord = zeros(Ned,4);
    
    for i=1:Ned
        v = map.cells.edges.edgesvert(i,:);
        coord(i,:) = [C.coord(v(1),1),C.coord(v(2),1),C.coord(v(1),2),C.coord(v(2),2)];
    end
    
    map.cells.edges.coord = coord;


end