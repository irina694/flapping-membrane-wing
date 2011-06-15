function ratio = determineEdgeAreaProperties(NoEdges,edgeG,vtxG,edges)

    eps = 1e-6;
    thick = edges.thicknessRatio;
    
    ratio = zeros(1,NoEdges);   
    edgeCoord = edges.coord;
        
    for i=1:NoEdges
        vert1 = vtxG(:,edgeG(1,i))';
        vert2 = vtxG(:,edgeG(2,i))';        
        aux1 = [vert1(1),vert2(1),vert1(2),vert2(2)];
        aux2 = [vert2(1),vert1(1),vert2(2),vert1(2)];        
        
        for k=1:NoEdges
            bol1 = abs(aux1-edgeCoord(k,:))<eps;
            bol2 = abs(aux2-edgeCoord(k,:))<eps;

            if (all(bol1) | all(bol2))
                break;
            end                   
        end  
        ratio(i) =  thick(k);
        
    end
end

