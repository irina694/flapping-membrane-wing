function edges=applyPrules2egdes(edges,Prules)
   
    Nedges = edges.Nedges;
               
    for Ned=1:Nedges;
        edgeRule = edges.symbol(Ned);
        
        edges.prodsymb{Ned}=[];         %edge symblos after applying the production rules
        edges.prodtdir{Ned}=[];         %edge directions after applying the production rules
        edges.MkPosition{Ned}=[];
        edges.MkSymbol{Ned}=[];            
        edges.MkThick{Ned}=[];
        
        if(edgeRule~='x')
            index = double(edgeRule)-64;        %rule index
            if(edges.direction(Ned)=='+')
                aux1 = Prules.directRule{index}.letters;
                aux2 = Prules.directRule{index}.dir;
                aux3 = Prules.directRule{index}.markposition;
                aux4 = Prules.directRule{index}.markletter;
                aux5 = Prules.directRule{index}.markthick;
            else
                aux1 = Prules.inverseRule{index}.letters;
                aux2 = Prules.inverseRule{index}.dir;
                aux3 = Prules.inverseRule{index}.markposition;
                aux4 = Prules.inverseRule{index}.markletter;
                aux5 = Prules.inverseRule{index}.markthick;
             end
            edges.prodsymb{Ned} = [edges.prodsymb{Ned}, aux1];
            edges.prodtdir{Ned} = [edges.prodtdir{Ned}, aux2];
            edges.MkPosition{Ned} = [edges.MkPosition{Ned},aux3];
            edges.MkSymbol{Ned} = [edges.MkSymbol{Ned},aux4];            
            edges.MkThick{Ned} = [edges.MkThick{Ned},aux5];
        else
            edges.prodsymb{Ned} = [edges.prodsymb{Ned}, 'x'];
            edges.prodtdir{Ned} = [edges.prodtdir{Ned}, '+'];           
        end
    end       
end
