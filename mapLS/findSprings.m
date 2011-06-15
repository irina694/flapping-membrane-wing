function  springs = findSprings(map);
    coord = map.cells.edges.coord;
    Ned = map.cells.edges.edges.Nedges;
    
    springs.index = map.cells.edges.edgesvert;
    springs.length = zeros(Ned,1);
    
    for i=1:Ned
       L = sqrt((coord(i,2)-coord(i,1))^2+(coord(i,4)-coord(i,3))^2);
       springs.index(i) = L;         
    end
      
end