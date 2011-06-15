function map = findBoundaryVert(map)
    global  globalVal 

    Ned = map.connectivity.edges.Nedges;
    type = map.connectivity.edges.type;
    edvt = map.connectivity.edges.edgesvert;
    vert = zeros(1,map.cells.vertices.Nvertices);
    
    for i=1:Ned
        if(type(i)=='B')
           vert(edvt(i,1)) = 1;
           vert(edvt(i,2)) = 1;                       
        else
           if (vert(edvt(i,1)) ~= 1), vert(edvt(i,1)) = 2; end
           if (vert(edvt(i,2)) ~= 1), vert(edvt(i,2)) = 2; end           
        end                
    end

    map.cells.vertices.boundary = find(vert == 1); %boundary vertices
    map.cells.vertices.interior = find(vert == 2); %interior vertices    
    map.cells.vertices.unusedBV = find(vert == 0); %unused boundary vertices
    
    aux = map.cells.vertices.coords(map.cells.vertices.boundary,:);
    ind = findVertByCoords('x',aux,[globalVal.Xmin,globalVal.Xmax]); 
    map.cells.vertices.boundaryV = map.cells.vertices.boundary(ind);
    ind = findVertByCoords('y',aux,[globalVal.Ymin,globalVal.Ymax]); 
    map.cells.vertices.boundaryH = map.cells.vertices.boundary(ind);
            
end
