for i=1:24
    
    nn = length(savestate{i}.mapPopulation);
    
    savestate3{i} = savestate{i};
    for j=1:nn
        
        aux = savestate{i}.mapPopulation{j};        
        aux = rmfield(aux,'connectivity');
        
        aux.Prules = rmfield(aux.Prules,{'directRule','inverseRule'});
       
        aux.cells.vertices.connectivity = sparse(aux.cells.vertices.connectivity);
        aux.cells = rmfield(aux.cells,{'edgesGenealogy','cellGenealogy'});
        aux.cells.edges = rmfield(aux.cells.edges,{'MkIndex','parentEdS','parentEdD'});
        
        savestate3{i}.mapPopulation{j} = aux;
    end
    
end