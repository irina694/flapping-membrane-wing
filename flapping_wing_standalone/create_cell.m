function[connectivity,error] = create_cell(x, Gvars, input)

%global Gvars input

%save('lastInd','x');

Nlevels = Gvars.Nlevels;

error = 0;

try
    lXsmap = Gvars.NinitEdges;
    Xsmap = x(1:lXsmap);
    
    [startMap,startDir] = extractsSmap(Xsmap, Gvars);
    map = initMapStruct(startMap,startDir,input);
    
    [map.Prules.letters,map.Prules.direction,map.Prules.thickRule,...
        map.Prules.elastRule,map.Prules.pressRule] = extractRules(x(lXsmap+1:end-Gvars.nB), Gvars);
    
    map.cells.vertices.type = 'FFFF';
    
    map = mapLS(map,Nlevels, Gvars);
    
    connectivity=removeVertices(map);
    connectivity.tip = x(end)/100;
    
catch
    Info = 'Map L-system error';
    error = 1;
    connectivity = [];
end