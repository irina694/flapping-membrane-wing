function cells=createNewCells(cells,matches,Prules)
    
    edges = cells.edges;
    
    %saves the parent edges
    edges.parentEdS = edges.symbol;
    edges.parentEdD = edges.direction;
    
    %updates the edge and vert sructs and crates edge genealogy
    [edgesGenealogy,edges] = updateEdgeStrct(edges);        
    cells.edges = edges;
    cells.edgesGenealogy = edgesGenealogy;
   
    %creates the new cells and the splitting edges
    cells = divideCells(cells,matches);
    
    %orientes the cells edges in a clockwise manner
    %cells = orientscelledges(cells);
    
    %updates cells.edgetype field based on the edge genealogy
    cells = updateEdType(cells);
    
    cells = updateEdProperties(cells,Prules,matches);
end

        

function cells = updateEdProperties(cells,Prules,matches)
    global Gvars

    NparentEd = length(cells.edgesGenealogy);
    Ned = cells.edges.Nedges;
    EdThick = -1*ones(1,Ned);
    EdElast = -1*ones(1,Ned);
    EdPress = -1*ones(1,Ned);
    
    NnewEd = 0;
    for i=1:NparentEd
        ParentThickR = cells.edges.thicknessRatio(i);
        ParentElastR = cells.edges.elastRatio(i);
        ParentPressR = cells.edges.pressRatio(i);
        ParentSymbol = cells.edges.parentEdS(i);
        childrenEd = cells.edgesGenealogy{i};                                
        NnewEd = NnewEd + length(childrenEd);
        
        if(ParentSymbol~='x')
            index = double(ParentSymbol)-64;
            thickRule = Prules.directRule{index}.thick;
            elastRule = Prules.directRule{index}.elast;
            pressRule = Prules.directRule{index}.press;

            EdThick(childrenEd) = limitValues(ParentThickR*thickRule,Gvars.MinThck,Gvars.MaxThck);            
            EdElast(childrenEd) = limitValues(ParentElastR*elastRule,Gvars.MinElst,Gvars.MaxElst);
            EdPress(childrenEd) = limitValues(ParentPressR*pressRule,Gvars.MinPrss,Gvars.MaxPrss);
        else
            EdThick(childrenEd) = ParentThickR;            
            EdElast(childrenEd) = ParentElastR;
            EdPress(childrenEd) = ParentPressR;
        end
    end
   
    %updates the new new edges that divide the cells
    NewEd = length(matches);
    for i=1:NewEd
        if(~isempty(matches{i}))
            NnewEd = NnewEd+1;
%             ParentThickR1 = matches{i}(8);
%             ParentThickR2 = matches{i}(9);            
            ParentThickR1 = cells.edges.thicknessRatio(matches{i}(6))*matches{i}(8);
            ParentThickR2 = cells.edges.thicknessRatio(matches{i}(7))*matches{i}(9);
            averageT = (ParentThickR1+ParentThickR2)/2;            
            if(averageT > Gvars.MaxThck)
                averageT = Gvars.MaxThck;
            elseif(averageT < Gvars.MinThck)
                averageT = Gvars.MinThck;
            end                            
            EdThick(NnewEd) = averageT;                                    

            aux1 = cells.edges.elastRatio(matches{i}(6))*matches{i}(8);
            aux2 = cells.edges.elastRatio(matches{i}(7))*matches{i}(9);
            averageT = (aux1+aux2)/2;            
            if(averageT > Gvars.MaxElst)
                averageT = Gvars.MaxElst;
            elseif(averageT < Gvars.MinElst)
                averageT = Gvars.MinThck;
            end                            
            EdElast(NnewEd) = averageT;                                    

            aux1 = cells.edges.pressRatio(matches{i}(6))*matches{i}(8);
            aux2 = cells.edges.pressRatio(matches{i}(7))*matches{i}(9);
            averageT = (aux1+aux2)/2;            
            if(averageT > Gvars.MaxPrss)
                averageT = Gvars.MaxPrss;
            elseif(averageT < Gvars.MinPrss)
                averageT = Gvars.MinPrss;
            end                            
            EdPress(NnewEd) = averageT;                                    

        end                       
    end
        
    cells.edges.thicknessRatio = EdThick;
    cells.edges.elastRatio = EdElast;
    cells.edges.pressRatio = EdPress;
end



function y = limitValues(y,minV,maxV)

    trim = y > maxV;
    y = ~trim.*y + trim.*maxV;            
    trim = y < minV;
    y = ~trim.*y + trim.*minV;

end



function cells = updateEdType(cells)

    NparentEd = length(cells.edgesGenealogy);
    Ned = cells.edges.Nedges;
    EdType = char(73*ones(1,Ned));  %all initialized as I = interior edges  
    
    for i=1:NparentEd
        ParentType = cells.edges.type(i);
        childrenEd = cells.edgesGenealogy{i};
        EdType(childrenEd) = ParentType;                
    end
   
    cells.edges.type = EdType;
end




function cells = divideCells(cells,matches)

    Ncells = cells.Ncells;
    edges = cells.edges;
    LastEd = edges.Nedges;
    LastCl = Ncells;
    celledges = cells.celledges;
    edgedir = cells.edgedir;
    edgesGenealogy = cells.edgesGenealogy;
    cellGenealogy = cell(1,Ncells);
    oldDir = edges.parentEdD;

    
    for i=1:Ncells
        if (~isempty(matches{i}))            
            %inserts the new edges in the edge structure
            LastEd = LastEd+1;
            LastCl = LastCl+1;
            edges.Nedges = LastEd;
            edges.symbol(LastEd) = char(mod(matches{i}(5),100));
            edges.direction(LastEd) = '+'; %char(floor(matches{i}(5)/100));
            

            %we need the egde index. matches contains the edges index in
            %the celledges list. To get rigth information about the edges
            %and it's vertices we need the edge number in the cell.edges
            %structure
            match = matches{i};
            match(1) = celledges{i}(match(1));
            match(3) = celledges{i}(match(3));
            
            %we retrieve the edges that are connected in the celular
            %division
            p1 = edgesGenealogy{match(1)}(match(2));
            p2 = edgesGenealogy{match(3)}(match(4));
            d1 = edges.direction(p1);
            d2 = edges.direction(p2);            
            
            if(d1=='+')
                x1 = edges.coord(p1,2);
                y1 = edges.coord(p1,4);
                v1 = edges.edgesvert(p1,2);                
            else
                x1 = edges.coord(p1,1);
                y1 = edges.coord(p1,3);
                v1 = edges.edgesvert(p1,1);                                               
            end
                
            if(d2=='+')
                x2 = edges.coord(p2,2);
                y2 = edges.coord(p2,4);
                v2 = edges.edgesvert(p2,2);                
            else
                x2 = edges.coord(p2,1);
                y2 = edges.coord(p2,3);
                v2 = edges.edgesvert(p2,1);                                               
            end
                
                edges.coord(LastEd,:) = [x1,x2,y1,y2];
                edges.edgesvert(LastEd,:) = [v1,v2];                                                              
                              
            %inserts the new cells in the cell structure
            fstcell = [celledges{i}(1:matches{i}(1)),LastEd,celledges{i}(matches{i}(3):end)];
            sndcell = [celledges{i}(matches{i}(1):matches{i}(3)),LastEd];
            fstcdir = [edgedir{i}(1:matches{i}(1)),edges.direction(LastEd),edgedir{i}(matches{i}(3):end)];
            scncdir = [edgedir{i}(matches{i}(1):matches{i}(3)),edges.direction(LastEd)];
            
            [edgedir{i},celledges{i}] = applyGenealogy2fstcell(fstcell,fstcdir,edges.direction,oldDir,LastEd,match,edgesGenealogy);            
            [edgedir{LastCl},celledges{LastCl}] = applyGenealogy2sndcell(sndcell,scncdir,edges.direction,oldDir,LastEd,match,edgesGenealogy);
                                    
            
            cellGenealogy{i}=[i LastCl];
        else                
            [edgedir{i},celledges{i}] = applyGenealogy2barrencell(celledges{i},edgedir{i},edges.direction,oldDir,edgesGenealogy);
            cellGenealogy{i}=i;
        end
    end    
    cells.edges = edges;
    cells.celledges = celledges;
    cells.edgedir = edgedir;
    cells.Ncells = LastCl; 
    cells.cellGenealogy = cellGenealogy;
    
    cellvert=cell(1,cells.Ncells);
    for i = 1:cells.Ncells
       Ned = length(cells.celledges{i});
       cellvert{i} = zeros(1,Ned);
       for j = 1:Ned
           edg = cells.celledges{i}(j);
           verts = cells.edges.edgesvert(edg,:);
           dir = cells.edges.direction(edg);
           dir2 = cells.edgedir{i}(j);
           if ((dir =='+' & dir2=='+') | (dir =='-' & dir2=='+'))
               cellvert{i}(j)=verts(1);
           elseif ((dir =='+' & dir2=='-') | (dir =='-' & dir2=='-'))
               cellvert{i}(j)=verts(2);
           end                      
       end        
    end
    
    cells.cellvert = cellvert;
       
end


function y = invertDirection(x)
   y = double(x);
   y = ~((y-43)/2);
   y = char(2*y+43);
   y = y(end:-1:1);
end



function y = applyParentDirection(x,dir)
    if dir == '+'        
       y=x;
    else
       y = double(x);
       y = ~((y-43)/2);
       y = char(2*y+43);    
    end
end



function [celldir,celledges] = applyGenealogy2sndcell(celledges,parentDir,direction,oldDir,splitEd,match,edgesGenealogy)
    Ned = length(celledges);
    aux=[];
    dir=[];
    
    for j=1:Ned
        ed=celledges(j);
        if(ed==match(1))
             if (parentDir(j) == oldDir(ed))
                childEd = edgesGenealogy{ed}(match(2)+1:end);
                aux = [aux,childEd];                
                dir = [dir,direction(childEd)];
             else                 
                 childEd = edgesGenealogy{ed}(1:match(2));
                 aux = [aux,childEd(end:-1:1)];
                 dir = [dir,invertDirection(direction(childEd))];
             end              
        elseif(ed==match(3))
             if (parentDir(j) == oldDir(ed))
                childEd = edgesGenealogy{ed}(1:match(4));
                aux = [aux,childEd];    
                dir = [dir,direction(childEd)];
             else
                 childEd = edgesGenealogy{ed}(match(4)+1:end);
                 aux = [aux,childEd(end:-1:1)];    
                 dir = [dir,invertDirection(direction(childEd))];
             end
        elseif(ed==splitEd)
            aux=[aux,splitEd];
            if (direction(splitEd)=='+')
                dir = [dir,'-'];
            else
                dir = [dir,'+'];
            end
        else
            childEd = edgesGenealogy{ed};
             if (parentDir(j) == oldDir(ed))                
                aux = [aux,childEd];
                dir = [dir,direction(childEd)];
             else
                 aux = [aux,childEd(end:-1:1)];
                 dir = [dir,invertDirection(direction(childEd))];
             end
        end
    end        
    celledges = aux;
    celldir = dir;
end


function [celldir,celledges] = applyGenealogy2fstcell(celledges,parentDir,direction,oldDir,splitEd,match,edgesGenealogy)
    Ned = length(celledges);
    aux=[];
    dir=[];
    
    for j=1:Ned
        ed=celledges(j);
        if(ed==match(1))
            if (parentDir(j) == oldDir(ed))
                childEd = edgesGenealogy{ed}(1:match(2));
                aux = [aux,childEd];
                dir = [dir,direction(childEd)];
             else
                 childEd = edgesGenealogy{ed}(match(2)+1:end);                
                 aux = [aux,childEd(end:-1:1)];
                 dir = [dir,invertDirection(direction(childEd))];
             end
        elseif(ed==match(3))
             if (parentDir(j) == oldDir(ed))
                childEd = edgesGenealogy{ed}(match(4)+1:end);
                aux = [aux,childEd];
                dir = [dir,direction(childEd)];
             else
                 childEd = edgesGenealogy{ed}(1:match(4));
                 aux = [aux,childEd(end:-1:1)];
                 dir = [dir,invertDirection(direction(childEd))];
             end
        elseif(ed==splitEd)
            aux = [aux,splitEd]; 
            dir = [dir,direction(splitEd)];
        else
            childEd = edgesGenealogy{ed};
            if (parentDir(j) == oldDir(ed))               
               aux = [aux,childEd];
               dir = [dir,direction(childEd)];
            else
               aux = [aux,childEd(end:-1:1)];
               dir = [dir,invertDirection(direction(childEd))];
            end
        end
    end        
    celledges=aux;
    celldir = dir;
end


function [celldir,celledges] = applyGenealogy2barrencell(celledges,parentDir,direction,oldDir,edgesGenealogy)
    Ned = length(celledges);
    aux=[];
    dir=[];
    for j=1:Ned
        ed=celledges(j);
        
        childEd = edgesGenealogy{ed};
        if (parentDir(j) == oldDir(ed))               
           aux = [aux,childEd];
           dir = [dir,direction(childEd)];
        else
           aux = [aux,childEd(end:-1:1)];
           dir = [dir,invertDirection(direction(childEd))];
        end            
    end        
    celledges=aux;
    celldir = dir;
end




function [edgesGenealogy,edges]=updateEdgeStrct(edges)

    LastEdge = edges.Nedges;
    Nedges = edges.Nedges;
    edgesGenealogy = cell(1,Nedges);
    
    for Ned=1:Nedges
        Ndiv = size(edges.NewCoord{Ned},1);
        edgesGenealogy{Ned} = zeros(1,Ndiv);
        parentDir = edges.direction(Ned);
        
%         if (edges.direction(Ned)=='+')
            %updates the frist new edge
            edgesGenealogy{Ned} = Ned;
            edges.symbol(Ned) = edges.prodsymb{Ned}(1);
            edges.direction(Ned) = edges.prodtdir{Ned}(1);
            if(parentDir==edges.prodtdir{Ned}(1))
                edges.coord(Ned,:) = edges.NewCoord{Ned}(1,:);                
                edges.edgesvert(Ned,:) = edges.NewEdgVert{Ned}(1,:);              
            else
                edges.coord(Ned,:) = [edges.NewCoord{Ned}(1,[2,1]),edges.NewCoord{Ned}(1,[4,3])];
                edges.edgesvert(Ned,:) = edges.NewEdgVert{Ned}(1,[2,1]);            
            end
            %if the edge is divided in more than 1, it updates the remainder cells 
            if(Ndiv > 1)            
                for j=2:Ndiv
                    LastEdge = LastEdge+1;
                    edgesGenealogy{Ned}(j)=LastEdge;
                    edges.symbol(LastEdge) = edges.prodsymb{Ned}(j);
                    edges.direction(LastEdge) = edges.prodtdir{Ned}(j);
                    if(parentDir==edges.prodtdir{Ned}(j))
                        edges.coord(LastEdge,:) = edges.NewCoord{Ned}(j,:);
                        edges.edgesvert(LastEdge,:) = edges.NewEdgVert{Ned}(j,:);
                    else
                        edges.coord(LastEdge,:) = [edges.NewCoord{Ned}(j,[2,1]),edges.NewCoord{Ned}(j,[4,3])];
                        edges.edgesvert(LastEdge,:) = edges.NewEdgVert{Ned}(j,[2,1]);            
                    end
                end
            end
            
            if (parentDir=='-')
                ii = length(edgesGenealogy{Ned}):-1:1;
                edgesGenealogy{Ned} = edgesGenealogy{Ned}(ii); 
            end
%         else
%              %updates the frist new edge
%              edgesGenealogy{Ned}(Ndiv) = Ned;
%              edges.symbol(Ned) = edges.prodsymb{Ned}(Ndiv);
%              edges.direction(Ned) = edges.prodtdir{Ned}(Ndiv);
%              if(edges.prodtdir{Ned}(Ndiv)=='+')         %there's a switch in the direction
%                  edges.coord(Ned,:) = [edges.NewCoord{Ned}(Ndiv,[2,1]),edges.NewCoord{Ned}(Ndiv,[4,3])];
%                  edges.edgesvert(Ned,:) = edges.NewEdgVert{Ned}(Ndiv,[2,1]);
%              else
%                  edges.coord(Ned,:) = edges.NewCoord{Ned}(Ndiv,:);
%                  edges.edgesvert(Ned,:) = edges.NewEdgVert{Ned}(Ndiv,:);
%              end
%  %            edges.edgesvert(Ned,:) = edges.NewEdgVert{Ned}(1,:);
%              if(Ndiv > 1)            
%                  for j=Ndiv-1:-1:1
%                      LastEdge = LastEdge+1;
%                      edgesGenealogy{Ned}(j)=LastEdge;
%                      edges.symbol(LastEdge) = edges.prodsymb{Ned}(j);
%                      edges.direction(LastEdge) = edges.prodtdir{Ned}(j);
%                      if(edges.prodtdir{Ned}(j)=='+')
%                          edges.coord(LastEdge,:) = [edges.NewCoord{Ned}(j,[2,1]),edges.NewCoord{Ned}(j,[4,3])];
%                          edges.edgesvert(LastEdge,:) = edges.NewEdgVert{Ned}(j,[2,1]);
%                      else
%                          edges.coord(LastEdge,:) = edges.NewCoord{Ned}(j,:);
%                          edges.edgesvert(LastEdge,:) = edges.NewEdgVert{Ned}(j,:);
%                      end
%                  end
%  
%                   ii = length(edgesGenealogy{Ned}):-1:1;
%                   edgesGenealogy{Ned} = edgesGenealogy{Ned}(ii);                        
%              end
%         end         
        
        
    end    
    edges.Nedges=LastEdge;    
    edges=rmfield(edges,{'NewCoord','MkPosition','prodsymb','prodtdir','MkSymbol','MkThick','NewEdgVert'});
end

