function C=removeBoundaryVertices(C)

    N=C.Ncell;
    Xmax=1.0;
    Xmin=-1.0;
    Ymax=1.0;
    Ymin=-1.0;
    
    for i=1:N        
        vert2del=[];
         
        cellvert = C.cell{i};
        Nvert = length(cellvert);
        
        xcoord = C.coord(cellvert,1)';
        ycoord = C.coord(cellvert,2)';

        indices = cellvert; %C.cells{i}.edges.index(:,1)';
        Bxmin=find(xcoord==Xmin);
        Bxmax=find(xcoord==Xmax);
        Bymin=find(ycoord==Ymin);
        Bymax=find(ycoord==Ymax);
             
        [aux1,fstV,lstV]=findVert2del(Bxmin,Nvert);
        vert2del=[vert2del,indices(aux1)]; 
        C = updateCedges(C,indices([fstV,aux1,lstV]));
        
        [aux1,fstV,lstV]=findVert2del(Bxmax,Nvert);
        vert2del=[vert2del,indices(aux1)];
        C = updateCedges(C,indices([fstV,aux1,lstV]));
        
        [aux1,fstV,lstV]=findVert2del(Bymin,Nvert);
        vert2del=[vert2del,indices(aux1)];
        C = updateCedges(C,indices([fstV,aux1,lstV]));
        
        [aux1,fstV,lstV]=findVert2del(Bymax,Nvert);
        vert2del=[vert2del,indices(aux1)];
        C = updateCedges(C,indices([fstV,aux1,lstV]));
        
        %remove vertices from cells
        C.cell{i}=removefromcell(vert2del,C.cell{i});
        
        %remove connections from the connectivity matrix and add the new ones
        C.matrix=updateCmatrix(vert2del,C.cell{i},C.matrix);        
    end
    C.edges = removeZerothEntries(C.edges);
    
end



function edges = removeZerothEntries(edges)
    Ned = length(edges.thicknessRatio);
    
    for i=Ned:-1:1
       if(edges.thicknessRatio(i) == 0)
          edges.edgesvert(i,:) = [];
          edges.coord(i,:) = [];
          edges.thicknessRatio(i) = [];
          edges.elastRatio(i) = [];
          edges.pressRatio(i) = [];
          edges.type(i) = [];
       end                
    end
    edges.Nedges = length(edges.thicknessRatio);
end


function C = updateCedges(C,vert)
    Nv = length(vert);
    if(Nv == 0) 
        return;
    end
    
    edgevert = C.edges.edgesvert;
    thick = zeros(1,Nv-1);
    elast = zeros(1,Nv-1);
    press = zeros(1,Nv-1);
    edgIn = zeros(1,Nv-1);
    newEdgVert = [vert(1),vert(end)];
    newEdgCoord = [C.coord(vert(1),1),C.coord(vert(end),1),C.coord(vert(1),2),C.coord(vert(end),2)];
    for i=1:Nv-1
        aux1 = find(edgevert(:,1) == vert(i) & edgevert(:,2) == vert(i+1));
        aux2 = find(edgevert(:,1) == vert(i+1) & edgevert(:,2) == vert(i));
        
        edg = [aux1,aux2];
        edgIn(i) = edg;
        thick(i) = C.edges.thicknessRatio(edg); 
        elast(i) = C.edges.elastRatio(edg); 
        press(i) = C.edges.pressRatio(edg);
    end
    
    C.edges.edgesvert(edgIn(1),:) = newEdgVert;
    C.edges.coord(edgIn(1),:) = newEdgCoord;
    C.edges.thicknessRatio(edgIn(1)) = mean(thick);
    C.edges.elastRatio(edgIn(1)) = 1/sum(1./elast);
    C.edges.pressRatio(edgIn(1)) = mean(press);

    

    edgIn = sort(edgIn(2:end));
    nn = length(edgIn);
    for i=nn:-1:1
        C.edges.edgesvert(edgIn(i),:) = [0 0];
        C.edges.coord(edgIn(i),:) = [0 0 0 0]; 
        C.edges.thicknessRatio(edgIn(i)) = 0;
    end
end


function matrix=updateCmatrix(vert2del,cell,matrix)

    if(isempty(vert2del))
        return;
    end

    nv=length(vert2del);
    nc=length(cell);
    nm=size(matrix,1);
    
    zerol=zeros(1,nm);
    zeroc=zeros(nm,1);
    for i=1:nv
        matrix(vert2del(i),:)=zerol;
        matrix(:,vert2del(i))=zeroc;        
    end
    
    for i=1:nc-1
        matrix(cell(i),cell(i+1))=1;
        matrix(cell(i+1),cell(i))=1;
    end
    i=nc;
    matrix(cell(i),cell(1))=1;
    matrix(cell(1),cell(i))=1;
       
end


function cell=removefromcell(x,cell)

    if(isempty(x))
        return;
    end
    
    nx=length(x);
    nc=length(cell);
    
    newcell=[];
    for i=1:nc
        removevert=false;
        for j=1:nx
            if(cell(i)==x(j))
               removevert=true;
               break;
            end
        end
        if(~removevert)
            newcell=[newcell,cell(i)];
        end
    end
    cell=newcell;
end


function [vert2del,fstV,lstV] = findVert2del(x,Nvert)

    vert2del=[];
    fstV = [];
    lstV = [];
    n=length(x);
    if(n <= 2)
        return;
    elseif (n == Nvert)
        vert2del = x(2:n-1);
        fstV = x(1);
        lstV = x(n); 
        return;        
    end

    foundFst = false;
    for i=2:n-1
        if(x(i)-x(i-1)==1 && x(i+1)-x(i)==1)       %delete
            vert2del=[vert2del,x(i)];
        elseif (~foundFst)
            fstV = x(i);
            foundFst = true;
        elseif (foundFst)
            lstV = x(i);
        end
    end
    
    if (isempty(fstV) && isempty(lstV))
        fstV = x(1);
        lstV = x(n);
    end
        
    if(x(1)==1 && x(n)==Nvert)        %cicle
        if(x(n-1)==Nvert-1)
            if (isempty(lstV)), lstV=1;
            else aux = fstV; fstV = lstV; lstV = aux; end
            vert2del=[vert2del,Nvert];            
        end
        if(x(2)==2)
            if (isempty(lstV)), lstV = fstV; fstV=Nvert;
            elseif(x(n-1)~=Nvert-1) aux = fstV; fstV = lstV; lstV = aux; end
            vert2del=[1,vert2del];
        end        
    end          
end

