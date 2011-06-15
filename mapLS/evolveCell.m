function map = evolveCell(Nlevels,map)
    global globalVal
    
    for k=1:Nlevels    
        map0 = map;
        map = mapLS(map,k);  
    end
    
    %creates the connectivity structure to be used later in the objective
    %function
    C.matrix = map.cells.vertices.connectivity;
    C.coord = map.cells.vertices.coords;
    C.Ncell = map.cells.Ncells;
    C.last = map.cells.vertices.Nvertices;
    C.cell = map.cells.cellvert;
    C.edges.thicknessRatio = map.cells.edges.thicknessRatio;
    C.edges.elastRatio = map.cells.edges.elastRatio;
    C.edges.pressRatio = map.cells.edges.pressRatio;
    C.edges.type = map.cells.edges.type;
    C.edges.coord = map.cells.edges.coord;
    C.edges.edgesvert = map.cells.edges.edgesvert;
    C.edges.Nedges = map.cells.edges.Nedges;

    map.connectivity = C;
    map = joinBoundaryEdges(map);    
    
    %computes the equilibrium
    map = findBoundaryVertSquare(map);
    map = computeCellEquilibrium(map,globalVal.RelaxF);
    map = updateMapstruct(map);
       
    %map.connectivity=removeBoundaryVertices(map.connectivity);

    map.connectivity.coord(:,1)=map.connectivity.coord(:,1)*globalVal.ScaleFx;    
    map.connectivity.coord(:,2)=map.connectivity.coord(:,2)*globalVal.ScaleFy;
    map.connectivity.edges.coord(:,1:2) = map.connectivity.edges.coord(:,1:2)*globalVal.ScaleFx;
    map.connectivity.edges.coord(:,3:4) = map.connectivity.edges.coord(:,3:4)*globalVal.ScaleFy;    
    
end