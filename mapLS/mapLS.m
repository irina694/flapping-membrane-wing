function map=mapLS(map,Nlevel)
%this function creates the Nlevel map L-System of an initial map following
%the production rules, given in Prules
%   Input:
%     map : a structure with the necessary fields:
%     Nlevels : the number of time the rules are applied
%   Output:
%     map : the map structure after applying the map L-S. The structure is
%     ready to be viewed with drawMap function.
%   

    Level = map.Level;                      %initial level of the map LS

    %prepares rules to be used later in the code, just for easier coding 
    map.Prules = preparesRules(map.Prules);
    
    for i=Level+1:Nlevel
        cells=map.cells;

        %divides the edges accordingly the prodution rules
        [cells.edges,cells.vertices] = divideEdges(cells.edges,cells.vertices,map.Prules);

        %finds where to divide the cell
        matches = findMatches(cells);

        %creates the children cells
        cells = createNewCells(cells,matches,map.Prules);    
        
        %updates the map structure with the new values
        map.cells = cells;
        map.Level = i;
    end
    
    %creates the vertices matrix of connectivity (redundant information but helpful for other operations)
    map=createVertexConnectivity(map);
end

