function fixedGeo = fixSmallEdge(geom,lMin)

smallN = 1e-10;

geomE = geomedit(geom);
adj = flgeomadj(geom,[1;1]);
nEd = flgeomnbs(geom);
iFlag = 0;

for i=1:nEd
    vtx_x = get(geomE{i},'ctrlx');
    vtx_y = get(geomE{i},'ctrly');
    vt_x1=vtx_x(1);
    vt_y1=vtx_y(1);
    vt_x2=vtx_x(end);
    vt_y2=vtx_y(end);
    lengthE=sqrt((vt_x2-vt_x1)^2+(vt_y2-vt_y1)^2);
    if(lengthE<lMin)
        % flag edge
        iFlag = iFlag+1;
        flagE(iFlag) = i;
        % delete the small edge
        geomE([i])={[]};
        % get adjacent edges
        adjE = find(adj{1}(:,i));
        lAdjE=length(adjE);
        if(lAdjE~=0)
            % calculate the reference point of the small edge
            vtxM_x = vt_x1;
            vtxM_y = vt_y1;
            % adjust the adjacent edges
            for j = 1:lAdjE
                jEd = adjE(j);
                vtx_xJ = get(geomE{jEd},'ctrlx');
                vtx_yJ = get(geomE{jEd},'ctrly');
                edgeJ = [vtx_xJ';vtx_yJ'];
                weightJ = get(geomE{jEd},'weight');
                d11 = max(abs(edgeJ(1,1) - vt_x1),abs(edgeJ(2,1) - vt_y1));
                d12 = max(abs(edgeJ(1,1) - vt_x2),abs(edgeJ(2,1) - vt_y2));
                d21 = max(abs(edgeJ(1,end) - vt_x1),abs(edgeJ(2,end) - vt_y1));
                d22 = max(abs(edgeJ(1,end) - vt_x2),abs(edgeJ(2,end) - vt_y2));
                
                dMin = min([d11,d12,d21,d22]);
                
                switch dMin
                    case d11
                        geomE{jEd} = curve2([vtxM_x edgeJ(1,2:end)],...
                            [vtxM_y edgeJ(2,2:end)],...
                            weightJ');
                    case d21
                        geomE{jEd} = beziercurve2([edgeJ(1,1:end-1) vtxM_x],...
                            [edgeJ(2,1:end-1) vtxM_y],...
                            weightJ');
                    case d12
                        geomE{jEd} = curve2([vtxM_x edgeJ(1,2:end)],...
                            [vtxM_y edgeJ(2,2:end)],...
                            weightJ');
                    case d22
                        geomE{jEd} = curve2([edgeJ(1,1:end-1) vtxM_x],...
                            [edgeJ(2,1:end-1) vtxM_y],...
                            weightJ');
                    otherwise
                        save geomBug geom;
                        fprintf('Error in small geometry.\n');
                end
            end
            % correct adjacency matrix for modified edges
            for j=1:lAdjE
                adj{1}(i,adjE(j)) = 0;
                adj{1}(adjE(j),i) = 0;
                for k=1:lAdjE
                    if(k~=j)
                        adj{1}(adjE(k),adjE(j)) = 1;
                    end
                end
            end
        end
    end
end

if(iFlag==0)
    fixedGeo = geom;
else
    fixedGeo = geomedit(geom,geomE);
end
