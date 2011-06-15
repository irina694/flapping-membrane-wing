function map=initMapStruct(startMap,startDir,input)

    if(nargin < 3)
        input.Nvertices = 4;
        input.vertCoord = [0 0; 1 0; 1 1; 0 1];
    end        
        

    cells.edges.symbol = startMap;
    cells.edges.direction = startDir;
    cells.edgedir{1} = startDir;
    
    Nvert = input.Nvertices;
    Nedges = Nvert;
    cells.vertices.Nvertices = Nvert;
    cells.vertices.coords = input.vertCoord;
    aux1 = [Nvert,1:Nvert-1];
    aux2 = [2:Nvert,1];
    cells.vertices.LR = [aux1',aux2'];
%    cells.vertices.LR = [6 2; 1 3; 2 4; 3 5; 4 6; 5 1];  
    cells.vertices.type = char(double('F')*ones(1,Nedges));
    
    cells.boundary.boundariesLR = [1:Nvert]*10+[2:Nvert,1];
%    cells.boundary.boundariesLR = [12 23 34 45 56 61]; 
    cells.boundary.boundLRFixed = [Nvert,1:Nvert]*10+[2:Nvert,1,2];
%    cells.boundary.boundLRFixed = [62 13 24 35 46 51 62];
   
%    cells.vertices.LR = [4 2; 1 3; 2 4; 3 1];    
%    cells.vertices.type = 'FFFF';
%    cells.vertices.boundariesLR = [12 23 34 41]; 
%    cells.vertices.boundLRFixed = [42 13 24 31 42];
%    cells.boundaryEq = {'y=0','x=1','y=1','x=0'};
    
    cells.edges.Nedges = Nedges;
    
    aux = 1:cells.edges.Nedges+1;
    aux(end) = 1;
    
    for i=1:Nedges
        if(startDir(i)=='+') 
            cells.edges.edgesvert(i,:) = [aux(i),aux(i+1)];
            cells.edges.coord(i,:) = [cells.vertices.coords(aux(i),1) cells.vertices.coords(aux(i+1),1) ...
                                      cells.vertices.coords(aux(i),2) cells.vertices.coords(aux(i+1),2)];
        else
            cells.edges.edgesvert(i,:) = [aux(i+1),aux(i)];
            cells.edges.coord(i,:) = [cells.vertices.coords(aux(i+1),1) cells.vertices.coords(aux(i),1) ...
                                      cells.vertices.coords(aux(i+1),2) cells.vertices.coords(aux(i),2)];            
        end        
    end
    
    aux = double('B')*ones(1,Nedges);
    cells.edges.type = char(aux);
    cells.edges.thicknessRatio = ones(1,Nedges);
    cells.edges.elastRatio = ones(1,Nedges);
    cells.edges.pressRatio = ones(1,Nedges);
    cells.Ncells = 1;
    cells.celledges = {1:Nedges};
    cells.cellvert = {1:Nvert};
    
    map.cells = cells;
    
    map.Level = 0;
end