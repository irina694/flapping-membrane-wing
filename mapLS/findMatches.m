function matches=findMatches(cells)
    global Gvars

    Ncells = cells.Ncells;
    edges=cells.edges;
    coords = cells.vertices.coords; 
    
    areas = computeAreas(cells);
    
    for Nc=1:Ncells            
        matches{Nc}=[];
        
        if (areas(Nc) > Gvars.MinArea)
            celledge = cells.celledges{Nc};
            cellvert = cells.cellvert{Nc};
            edgedir = cells.edgedir{Nc};

            Nedges = length(celledge);

            foundmatch=false;
            %for Ned=1:Nedges-1
            Ned = 0;
            while (~foundmatch & Ned < Nedges-1)
                Ned = Ned+1;            
                markers = edges.MkPosition{celledge(Ned)};

                if(~isempty(markers))
                    MkSymbol = edges.MkSymbol{celledge(Ned)};
                    ParentDir = edgedir(Ned); %applyParentDirection(edgedir(Ned),edges.direction(celledge(Ned)));
                    %we need the global index of the edges
                    ed1Index = celledge(Ned);

                    nn = length(markers);
                    for i=1:nn
                        symbolDir = char(floor(MkSymbol(i)/100));                      
                        if (ParentDir == symbolDir)                        
                            symbol = char(mod(MkSymbol(i),100));
                            candidates = lookforcandidates(edges,celledge(Ned+1:end),edgedir(Ned+1:end),symbol,Ned);
                            if(~isempty(candidates))                             
                                match=testCandidates(edges,ed1Index,markers(i),candidates,areas(Nc),coords(cellvert,:));                    
                                if(~isempty(match))
                                    %match return the egde number, but we need the edge index
                                    index = find(celledge==match(1));
                                    %edge thickness rule
                                    aux = [edges.MkThick{ed1Index}(i),edges.MkThick{match(1)}(match(3))];
                                    matches{Nc}=[Ned,markers(i),index,edges.MkPosition{match(1)}(match(3)),...
                                                 double(symbol),ed1Index,match(1),aux];
                                    foundmatch=true;
                                    break;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end



function match=testCandidates(edges,originEd,originMk,candidates,area0,coords)  
    
    global Gvars
    
    match=[];

    epsilon = 2*0.1745;  %the edges cant be colinear and the angle must be bigger than 0.0873 rad = 20 deg
    
    P0 = [edges.NewCoord{originEd}(originMk,2),edges.NewCoord{originEd}(originMk,4)];
    v0 = [edges.NewCoord{originEd}(originMk,2)-edges.NewCoord{originEd}(originMk,1),...
          edges.NewCoord{originEd}(originMk,4)-edges.NewCoord{originEd}(originMk,3)];
    nV0 = norm(v0);
    mp0 = [edges.NewCoord{originEd}(originMk,2)+edges.NewCoord{originEd}(originMk,1),...
           edges.NewCoord{originEd}(originMk,4)+edges.NewCoord{originEd}(originMk,3)]/2;

    Ncandidates=size(candidates,1);
    for i=1:Ncandidates
        ed = candidates(i,1);
        mk = candidates(i,2);        
        
        P1 = [edges.NewCoord{ed}(mk,2),edges.NewCoord{ed}(mk,4)];          
        range = candidates(i,4)+1:candidates(i,5);
        newCellCoords = [P0;coords(range,:);P1];
        Sarea = computeSingleArea(newCellCoords);
        if (Sarea/area0 > Gvars.MinAreaRatio & Sarea/area0 < (1-Gvars.MinAreaRatio))
            v1 = [edges.NewCoord{ed}(mk,2)-edges.NewCoord{ed}(mk,1),...
                  edges.NewCoord{ed}(mk,4)-edges.NewCoord{ed}(mk,3)];          
            nV1 = norm(v1);
            mp1 = [edges.NewCoord{ed}(mk,2)+edges.NewCoord{ed}(mk,1),...
                   edges.NewCoord{ed}(mk,4)+edges.NewCoord{ed}(mk,3)]/2;
            v01 = mp1-mp0;   
            nV01 = norm(v01);

            angle = acos(dot(v0,v1)/(nV1*nV0));
            angle0 = acos(dot(v0,v01)/(nV01*nV0));
            angle1 = acos(dot(v1,v01)/(nV1*nV01));
            if(~(abs(angle0)<epsilon | abs(angle0-pi)<epsilon | abs(angle1)<epsilon | abs(angle1-pi)<epsilon));
                match = [ed,mk,candidates(i,3)];  
                break
            end
        end
    end    
end



function candidates=lookforcandidates(edges,celledges,edgedir,symbol,StartingEd)

    Nedges=length(celledges);
    candidates=[];
    
    for i=1:Nedges
        MkSymbol = char(mod(edges.MkSymbol{celledges(i)},100));
        symbolDir = char(floor(edges.MkSymbol{celledges(i)}/100));
        nn = length(MkSymbol);        
        sndEdDir = edgedir(i); %applyParentDirection(edgedir(i),edges.direction(celledges(i)));
        
        for j=1:nn
            if (symbol==MkSymbol(j) & symbolDir(j)==sndEdDir)
                candidates=[candidates;celledges(i),edges.MkPosition{celledges(i)}(j),edges.MkIndex{celledges(i)}(j),StartingEd,StartingEd+i];
            end           
        end                
    end 
end


function areas = computeAreas(cells)
    Ncell = cells.Ncells;
    
    areas = zeros(1,Ncell);
    coord = cells.vertices.coords;    
    for i=1:Ncell
        cellInd = cells.cellvert{i};
        aux=0;
        nn = length(cellInd);
        for j=1:nn-1
            aux=aux+(coord(cellInd(j),1)-coord(cellInd(j+1),1))*...
                        (coord(cellInd(j),2)+coord(cellInd(j+1),2))/2;
        end
        aux=aux+(coord(cellInd(nn),1)-coord(cellInd(1),1))*...
                    (coord(cellInd(nn),2)+coord(cellInd(1),2))/2;
        areas(i)=abs(aux);
    end
    
end


function y = computeSingleArea(coords)
        aux=0;
        nn = size(coords,1);
        for j=1:nn-1
            aux=aux+(coords(j,1)-coords(j+1,1))*(coords(j,2)+coords(j+1,2))/2;
        end
        aux=aux+(coords(nn,1)-coords(1,1))*(coords(nn,2)+coords(1,2))/2;
        y = abs(aux);

end

