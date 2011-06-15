function map = createVertexConnectivity(map)
    
    edges = map.cells.edges;   
    vert = map.cells.vertices;
    Nvertices = vert.Nvertices;
    Nedges = edges.Nedges;
    
    VertConnct = zeros(Nvertices,Nvertices);
    for i=1:Nedges
        ii = edges.edgesvert(i,1);
        jj = edges.edgesvert(i,2);
        
        VertConnct(ii,jj) = 1;
        VertConnct(jj,ii) = 1;                  
    end
    
    map.cells.vertices.connectivity = VertConnct;
    
end