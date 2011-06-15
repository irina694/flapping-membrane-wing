function connectivity=removeVertices(mapIn)

vertices.type = mapIn.cells.vertices.type;
vertices.coords=mapIn.cells.vertices.coords;
vertices.Nvertices = mapIn.cells.vertices.Nvertices;

edges.type = mapIn.cells.edges.type;
edges.thicknessRatio = mapIn.cells.edges.thicknessRatio;
edges.vertices = mapIn.cells.edges.edgesvert;
edges.Nedges = mapIn.cells.edges.Nedges;

cells.edges = mapIn.cells.celledges;
cells.Ncells = mapIn.cells.Ncells;

Nedges = edges.Nedges;

Nvtx = vertices.Nvertices;
for j = 1:Nvtx
    
    % don't delete fixed vertices (F)
    if vertices.type(j)~='F'
        % find the adjacent edges to vertex j
        adjEdg2Vtx = find(edges.vertices(:,:)==j);
        nAdj=length(adjEdg2Vtx);
        if nAdj < 3
            % flag vertices that have two or less edges adjacent to it
            vertices.type(j)='D';
            % mark edge with higher index to deletion and make the edge
            % collapse to the deleted vertex
            if nAdj == 2
                Iedge1 = mod(adjEdg2Vtx(1),Nedges);
                Jedge1 = ceil(adjEdg2Vtx(1)/Nedges);
                Iedge2 = mod(adjEdg2Vtx(2),Nedges);
                Jedge2 = ceil(adjEdg2Vtx(2)/Nedges);
                if Iedge1 < Iedge2
                    [Iedge1 Iedge2] = swap(Iedge1,Iedge2);
                    [Jedge1 Jedge2] = swap(Jedge1,Jedge2);
                end
                NJedge2=1+~(Jedge2-1);
                edges.vertices(Iedge1,Jedge1) = edges.vertices(Iedge2,NJedge2);
                edges.vertices(Iedge2,NJedge2) = edges.vertices(Iedge2,Jedge2);
                edges.type(Iedge2) = 'D';
            else
                fprintf('Open end vertice %i',j);
            end
        end
    end
end

% update cells, edges and vertices
%
% delete vertices marked as 'D' and create mapping between vertices before
% and after deletion
verticeQ = vertices.type~='D';

vertices.type=vertices.type(verticeQ);
vertices.coords=vertices.coords(verticeQ,:);
vertices.Nvertices = length(vertices.type);

% delete edges marked as 'D'
edgeQ = edges.type~='D';

edges.type = edges.type(edgeQ);
edges.thicknessRatio = mapIn.cells.edges.thicknessRatio(edgeQ);
edges.vertices = edges.vertices(edgeQ,:);
edges.Nedges = length(edges.type);

% renumber the vertices in the updated edges
bef2aftVtx = cumsum(verticeQ);
edges.vertices = bef2aftVtx(edges.vertices(:,:));

% update cells
bef2aftEdg = cumsum(edgeQ);
for k = 1:cells.Ncells
    cells.edges{k} = bef2aftEdg(cells.edges{k}(edgeQ(cells.edges{k})));
end

connectivity.vertices = vertices;
connectivity.edges = edges;
connectivity.cells = cells;

end


