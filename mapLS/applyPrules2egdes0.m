function edges=applyPrules2egdes0(edges,Prules)
   
    Nedges = edges.Nedges;
    
       
    %rules
    rules = Prules.letters;
    direction = Prules.direction;
    thickRule = Prules.thickRule;
    
    for Ned=1:Nedges;
        edgeRule = edges.symbol(Ned);
        
        edges.prodsymb{Ned}=[];         %edge symblos after applying the production rules
        edges.prodtdir{Ned}=[];         %edge directions after applying the production rules
        
        if(edgeRule~='x')
            index = double(edgeRule)-64;
            aux1 = rules{index};
            aux2 = direction{index};
            if(edges.direction(Ned)=='-')
               [aux1,aux2]=invertrule(aux1,aux2);
            end
            edges.prodsymb{Ned} = [edges.prodsymb{Ned}, aux1];
            edges.prodtdir{Ned} = [edges.prodtdir{Ned}, aux2];
        else
            edges.prodsymb{Ned} = [edges.prodsymb{Ned}, 'x'];
            edges.prodtdir{Ned} = [edges.prodtdir{Ned}, '+'];
        end
    end
    
    

end


function [invrule,invdir]=invertrule(rule,dir)

    nn=length(rule);
    
    invrule=[];
    invdir=[];
    i=nn;
    while (i>0)
        if(rule(i)==']')
            if(rule(i-2)=='+')    
                invrule=[invrule,'[-',rule(i-1),']'];
            else
                invrule=[invrule,'[+',rule(i-1),']'];
            end
            
            if(dir(i-1)=='+')    
                invdir=[invdir,'00-0'];
            else
                invdir=[invdir,'00+0'];
            end
            %invdir=[invdir,'00',dir(i-1),'0'];
            i=i-4;
        else
            invrule=[invrule,rule(i)];
            if(dir(i)=='+')    
                invdir=[invdir,'-'];
            else
                invdir=[invdir,'+'];
            end
            %invdir=[invdir,dir(i)];
            i=i-1;
        end
    end

end