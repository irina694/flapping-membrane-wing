function map=computeCellEquilibrium(map,relaxF,output)   

    connectivity = map.connectivity;
    epsilon = 1;
    smallN = 1e-5;
    residual = 1e12;
    residual0 = 1e12;
    maxDXYeps = 5e-3;
    maxDXY = 1;
    maxIt = 10; 
    geom=[];
    %finds sprinds nodes ,length and K
    springs = findSprings(map);
    %determines which nodes are located at the boundary
%     geom.externalX = map.cells.vertices.boundaryV;
%     geom.externalY = map.cells.vertices.boundaryH;
     
    connectivity0=connectivity;
        
    geom = computeAreaAndNormal(connectivity,geom);    
    geom = findCellPressureCoef(map.cells,geom);
        
    it=0;
    ntry=0;
    maxtry=5;
    maxXY=max(connectivity.coord);
    minXY=min(connectivity.coord);
    maxXY0 = maxXY;
    minXY0 = minXY;
    
    while (residual>epsilon & maxDXY > maxDXYeps & it <= maxIt)  
        if(nargin == 3)
            if(output.plot==1 & mod(it,output.rate)==0)
                i=it;
                drawConnectivity(connectivity); 
                disp(residual);
                disp(maxDXY);
                disp(it);                          
                pause(0.5);
                if (output.out2file == 1)
                   masterName=[output.fileName,num2str(level)];
                   if i<10
                      file=[masterName,'00',num2str(i),'.tiff']; 
                   elseif i<100
                       file=[masterName,'0',num2str(i),'.tiff'];
                   else
                       file=[masterName,num2str(i),'.tiff'];
                   end
                   X=getframe(gcf);
                   imwrite(X.cdata,file);          
                end
            end 
        end
        
        %limits = (maxXY(1)>(maxXY0(1)+smallN) | maxXY(2)>(maxXY0(2)+smallN) | minXY(1)<(minXY0(1)-smallN) | minXY(2)<(minXY0(2)-smallN));
        if(residual>residual0*10 | ~chechOverlapping(connectivity) | ~checkMinAngle(connectivity))   %it's diverging
            %disp('restart');
            ntry=ntry+1;
            connectivity=connectivity0;
            
            maxCoord = max(max(maxXY),max(abs(minXY)));
            aux = min(max(maxXY0)/maxCoord,0.1);            
            
            relaxF = relaxF*aux; 
            if(ntry>maxtry)  
            %    disp(residual);
            %    disp(it);
            %    disp(sum(geom.area));
                map.connectivity = connectivity;
                return;
            end
        else
            residual0=residual;
            connectivity0=connectivity;
        end
            
        it = it+1;
        %initialize forces
        force.X=zeros(1,connectivity.last);
        force.Y=zeros(1,connectivity.last);
        X0=connectivity.coord(:,1)';
        Y0=connectivity.coord(:,2)';

        %computes the area of each cell and exterior normals
        geom = computeAreaAndNormal(connectivity,geom); 
        force = computePressure(connectivity,geom,force);
        force = computeSpringForce(springs,connectivity.coord,force);
        
%         %zeros in all nodes of the boundary;
%         %force = boundaryReaction(geom,force);
%         force.X(geom.externalX) = 0;
%         force.Y(geom.externalY) = 0;

%         X1 = force.X*relaxF + X0;
%         Y1 = force.Y*relaxF + Y0;    

        [X1,Y1] = updateVert(map,force,X0,Y0,relaxF);
        
        deltaX = abs(X1-X0);
        deltaY = abs(Y1-Y0);
        maxDXY = max(max(deltaX),max(deltaY));
        
        residualVec = [force.X';force.Y'];  %negativeAreaPenalisation*geom.overlap];
        residual = norm(residualVec);        
        connectivity.coord=[X1;Y1]';
        
        maxXY=max(connectivity.coord);
        minXY=min(connectivity.coord);

    end   
    
    map.connectivity = connectivity;
end


%==========================================================================
%======================        FUNCTIONS        ===========================
%==========================================================================
function [X1,Y1] = updateVert0(map,force,X0,Y0,relaxF)
    
    Nvert = map.cells.vertices.Nvertices;
    type = map.cells.vertices.type;
    vLR = map.cells.vertices.LR;
    X1 = zeros(1,Nvert);
    Y1 = zeros(1,Nvert);
    coords = map.cells.vertices.coords;
    for i=1:Nvert
        if (type(i) == 'B')
            eEdge = coords(vLR(i,2),:)-coords(vLR(i,1),:);
            eEdge = eEdge/norm(eEdge);            
            theta2 = atan2(eEdge(2),eEdge(1));
            if(theta2<0), theta2 = 2*pi+theta2;end
            
            if (theta2 == 0 || theta2 == pi)
                deltaX = force.X(i)*relaxF;
                deltaY = 0;
            elseif (theta2 == pi/2 || theta2 == -pi/2)
                deltaX = 0;
                deltaY = force.Y(i)*relaxF;
            else       
                eForce = [force.X(i),force.Y(i)];
                NForce = norm(eForce);
                eForce = eForce/NForce;           %unit force vector
                theta1 = atan2(eForce(2),eForce(1));
                if(theta1<0), theta1 = 2*pi+theta1;end
                theta = theta1-theta2;
%                 if(theta > -pi/2 && theta < pi/2)
                     tanForce = eEdge*cos(theta)*NForce;               
%                 else
%                     tanForce = -eEdge*cos(theta)*NForce;
%                 end

                deltaX = tanForce(1)*relaxF;
                deltaY = tanForce(2)*relaxF;
            end
            X1(i) = deltaX + X0(i);
            Y1(i) = deltaY + Y0(i);
        elseif(type(i) == 'I')
            X1(i) = force.X(i)*relaxF + X0(i);
            Y1(i) = force.Y(i)*relaxF + Y0(i); 
        elseif(type(i) == 'F')
            X1(i) = X0(i);
            Y1(i) = Y0(i);
        end        
    end

end

function [X1,Y1] = updateVert(map,force,X0,Y0,relaxF)
    
    boundariesLR = map.cells.boundary.boundariesLR;
    Nvert = map.cells.vertices.Nvertices;
    type = map.cells.vertices.type;
    vLR = map.cells.vertices.LR;
    X1 = zeros(1,Nvert);
    Y1 = zeros(1,Nvert);
    for i=1:Nvert
        if (type(i) == 'B')
            deltaX = force.X(i)*relaxF;
            deltaY = force.Y(i)*relaxF;
            Bid = find(boundariesLR == (vLR(i,1)*10+vLR(i,2)));
            [X1(i),Y1(i)] = changeVert([X0(i),Y0(i)],[deltaX,deltaY],Bid);
        elseif(type(i) == 'I')
            X1(i) = force.X(i)*relaxF + X0(i);
            Y1(i) = force.Y(i)*relaxF + Y0(i); 
        elseif(type(i) == 'F')
            X1(i) = X0(i);
            Y1(i) = Y0(i);
        end        
    end

end

%for now stays here but this must be external, provided by the user
function [x1,y1] = changeVert(vert,delta,Bid)    
    smallN = 1e-6;
    x1 = vert(1);
    y1 = vert(2);    
    if(norm(vert)>1+smallN)
        return
    end
    
    if (Bid==1)
        auxX = vert(1)+delta(1);
        if (auxX > 0 && auxX < 1)
            vert(1) = auxX;
        end
        auxY = vert(2)+delta(2);
        refY = -sqrt(1-vert(1)^2);
        if (auxY < 0 && auxY >= refY)
            vert(2) = auxY;
        end
    elseif (Bid==2)
        auxX = vert(1)+delta(1);
        if (auxX > 0 && auxX < 1)
            vert(1) = auxX;
        end
        auxY = vert(2)+delta(2);
        refY = sqrt(1-vert(1)^2);
        if (auxY > 0 && auxY <= refY)
            vert(2) = auxY;
        end
    elseif (Bid==3)
        auxX = vert(1)+delta(1);
        if (auxX < 0 && auxX > -1)
            vert(1) = auxX;
        end
        auxY = vert(2)+delta(2);
        refY = sqrt(1-vert(1)^2);
        if (auxY > 0 && auxY <= refY)
            vert(2) = auxY;
        end
    elseif (Bid==4)
        auxX = vert(1)+delta(1);
        if (auxX < 0 && auxX > -1)
            vert(1) = auxX;
        end
        auxY = vert(2)+delta(2);
        refY = -sqrt(1-vert(1)^2);
        if (auxY < 0 && auxY >= refY)
            vert(2) = auxY;
        end
    end
    x1 = vert(1);
    y1 = vert(2);
end


function geom = findCellPressureCoef(cells,geom)
    Ncell = cells.Ncells;
    
    pressCoef = zeros(1,Ncell);
    for i=1:Ncell
        celledg = cells.celledges{i};        
        aux = cells.edges.pressRatio(celledg);        
        pressCoef(i) = mean(aux);        
    end
    
    geom.pressCoef = pressCoef;
end



function force = computeSpringForce(springs,coords,force)
    global globalVal
    
    L0=globalVal.L0;
    K=globalVal.K;
    %springs.angle = [];
    Nsprings = size(springs.index,1);
    index = springs.index; 
    for i=1:Nsprings
        %vector = [coords(index(i,2),1)-coords(index(i,1),1), coords(index(i,2),2)-coords(index(i,1),2)];      
        theta=computeangle([coords(index(i,1),1),coords(index(i,1),2)],[coords(index(i,2),1),coords(index(i,2),2)]);
        L = sqrt((coords(index(i,1),1)-coords(index(i,2),1))^2+(coords(index(i,1),2)-coords(index(i,2),2))^2);        
        springs.angle(i,1) = theta;
        springs.length(i,1) = L;

        jp0 = index(i,1);
        jp1 = index(i,2);

        aux = -K*springs.Kratio(i)*(L-L0)*cos(theta);
        force.X(jp0)=force.X(jp0)-aux;
        force.X(jp1)=force.X(jp1)+aux;
        aux = -K*springs.Kratio(i)*(L-L0)*sin(theta);
        force.Y(jp0)=force.Y(jp0)-aux;
        force.Y(jp1)=force.Y(jp1)+aux;
    
    end
end


function theta=computeangle(P1,P2)
    vector = [P2(1)-P1(1), P2(2)-P1(2)];        
    theta=atan2(vector(2),vector(1));    
end


function force = computePressure(connectivity,geom,force)
    global globalVal
    P0=globalVal.P0;
    
    Ncell = connectivity.Ncell;
    
    for i=1:Ncell
        cellInd = connectivity.cell{i};
        nn = length(cellInd);
        pressure=P0/(geom.area(i))^geom.pressCoef(i);
        for j=1:nn-1
            jp0 = cellInd(j);
            jp1 = cellInd(j+1);
            theta = atan2(geom.normal{i}(j,2),geom.normal{i}(j,1));
            aux = pressure/2*cos(theta);
            force.X(jp0)=force.X(jp0)+aux;
            force.X(jp1)=force.X(jp1)+aux;
            aux = pressure/2*sin(theta);
            force.Y(jp0)=force.Y(jp0)+aux;
            force.Y(jp1)=force.Y(jp1)+aux;
        end
        jp0 = cellInd(nn);
        jp1 = cellInd(1);
        theta = atan2(geom.normal{i}(nn,2),geom.normal{i}(nn,1));
        aux = pressure/2*cos(theta);
        force.X(jp0)=force.X(jp0)+aux;
        force.X(jp1)=force.X(jp1)+aux;
        aux = pressure/2*sin(theta);
        force.Y(jp0)=force.Y(jp0)+aux;
        force.Y(jp1)=force.Y(jp1)+aux;
    end
end



function  springs = findSprings(map)
    coord = map.connectivity.edges.coord;
    Ned = map.connectivity.edges.Nedges;
    
    springs.index = map.connectivity.edges.edgesvert;
    springs.Kratio = map.connectivity.edges.elastRatio;  
    springs.length = zeros(Ned,1);       
    if (max(springs.Kratio)>1000)
        disp('bad');
    end
    
    for i=1:Ned
       L = sqrt((coord(i,2)-coord(i,1))^2+(coord(i,4)-coord(i,3))^2);
       springs.length(i) = L;                
    end
      
end



function geom = computeAreaAndNormal(connectivity,geom)
    Ncell = connectivity.Ncell;
    
    coord = connectivity.coord;    
    for i=1:Ncell
        cellInd = connectivity.cell{i};
        aux=0;
        normal=[];
        nn = length(cellInd);
        for j=1:nn-1
            aux=aux+(coord(cellInd(j),1)-coord(cellInd(j+1),1))*...
                        (coord(cellInd(j),2)+coord(cellInd(j+1),2))/2;
            nvect=[(coord(cellInd(j+1),2)-coord(cellInd(j),2)),...
                       -(coord(cellInd(j+1),1)-coord(cellInd(j),1))];
            normal=[normal;nvect/norm(nvect)];

        end
        aux=aux+(coord(cellInd(nn),1)-coord(cellInd(1),1))*...
                    (coord(cellInd(nn),2)+coord(cellInd(1),2))/2;
        area(i)=abs(aux);
        nvect=[(coord(cellInd(1),2)-coord(cellInd(nn),2)),...
                   -(coord(cellInd(1),1)-coord(cellInd(nn),1))];
        normal=[normal;nvect/norm(nvect)];        
        geom.normal{i}=normal;
    end
    
    geom.area=area;
end



