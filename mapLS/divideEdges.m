function [edges,vert]=divideEdges(edges,vert,Prules)

    %Applies production rules to the edges
    edges = applyPrules2egdes(edges,Prules);
    
    Nedges = edges.Nedges;
    Nvert = vert.Nvertices;
    for Ned=1:Nedges;
        
        edges.NewCoord{Ned} = [];
        NewEdVt = edges.edgesvert(Ned,1); 

        Ndiv = length(edges.prodsymb{Ned});
       
        dx = (edges.coord(Ned,2)-edges.coord(Ned,1))/Ndiv;
        dy = (edges.coord(Ned,4)-edges.coord(Ned,3))/Ndiv;

        for j=1:Ndiv
            Newedge = [edges.coord(Ned,1)+dx*(j-1), edges.coord(Ned,1)+dx*(j),...
                       edges.coord(Ned,3)+dy*(j-1), edges.coord(Ned,3)+dy*(j)];
            
            edges.NewCoord{Ned} = [edges.NewCoord{Ned}; Newedge];
                         
            
            if(j>1)
                Nvert = Nvert+1;
                NewEdVt = [NewEdVt,Nvert];
                
                vert.Nvertices = Nvert;
                vert.coords(Nvert,:) = [edges.coord(Ned,1)+dx*(j-1) edges.coord(Ned,3)+dy*(j-1)];
                
                if (edges.direction(Ned) == '+') 
                    LRvert = edges.edgesvert(Ned,:);
                else
                    LRvert = edges.edgesvert(Ned,[2,1]);
                end
                                    
                if ((vert.type(LRvert(1)) == 'B' || vert.type(LRvert(1)) == 'F') && ...
                    (vert.type(LRvert(1)) == 'B' || vert.type(LRvert(1)) == 'F') && edges.type(Ned) == 'B')
                    vert.type(Nvert) = 'B';
                    vert.LR(Nvert,:) = [vert.LR(LRvert(2),1),vert.LR(LRvert(1),2)];
                else
                    vert.type(Nvert) = 'I';
                    vert.LR(Nvert,:) = [0,0]; 
                end
            end              
        end
        NewEdVt = [NewEdVt,edges.edgesvert(Ned,2)];
        edges.NewEdgVert{Ned} = reshapeVec(NewEdVt);        
    end   


end




function edges=applyPrules2egdes(edges,Prules)
   
    Nedges = edges.Nedges;
               
    for Ned=1:Nedges;
        edgeRule = edges.symbol(Ned);
        
        edges.prodsymb{Ned}=[];         %edge symblos after applying the production rules
        edges.prodtdir{Ned}=[];         %edge directions after applying the production rules
        edges.MkPosition{Ned}=[];
        edges.MkSymbol{Ned}=[];            
        edges.MkThick{Ned}=[];
        edges.MkIndex{Ned}=[];
        
        if(edgeRule~='x')
            index = double(edgeRule)-64;        %rule index
            if(edges.direction(Ned)=='+')
                aux1 = Prules.directRule{index}.letters;
                aux2 = Prules.directRule{index}.dir;
                aux3 = Prules.directRule{index}.markposition;
                aux4 = Prules.directRule{index}.markletter;
                aux5 = Prules.directRule{index}.markthick;
                aux6 = Prules.directRule{index}.markindex;
            else
                aux1 = Prules.inverseRule{index}.letters;
                aux2 = Prules.inverseRule{index}.dir;
                aux3 = Prules.inverseRule{index}.markposition;
                aux4 = Prules.inverseRule{index}.markletter;
                aux5 = Prules.inverseRule{index}.markthick;
                aux6 = Prules.inverseRule{index}.markindex;
             end
            edges.prodsymb{Ned} = [edges.prodsymb{Ned}, aux1];
            edges.prodtdir{Ned} = [edges.prodtdir{Ned}, aux2];
            edges.MkPosition{Ned} = [edges.MkPosition{Ned},aux3];
            edges.MkSymbol{Ned} = [edges.MkSymbol{Ned},aux4];            
            edges.MkThick{Ned} = [edges.MkThick{Ned},aux5];
            edges.MkIndex{Ned} = [edges.MkIndex{Ned},aux6];
        else
            edges.prodsymb{Ned} = [edges.prodsymb{Ned}, 'x'];
            edges.prodtdir{Ned} = [edges.prodtdir{Ned}, '+'];           
        end
    end       
end




function mat = reshapeVec(vec)
    
    N = length(vec);
    
    for i=1:N-1
        mat(i,:) = vec(i:i+1);
    end

end

