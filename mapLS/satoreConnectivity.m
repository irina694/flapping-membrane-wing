for i=1:24
    aux{i}=savestate{i};
    aux{i}=rmfield('mapPopulation',aux{i});
    
    
    aux2 = savestate{1}.mapPopulation{1}.connectivity;
    aux2.matrix = sparse(aux2.matrix);
    aux{i}.mapPopulation