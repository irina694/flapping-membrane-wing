function[x,y,z,LE_pp,TE_pp,weight,points,x_gauss,y_gauss,z_gauss,LE_gauss,b_gauss,...
    connectivity,node,edge,x_FEA,y_FEA,z_FEA,trielements,beam_el,beam_ID] = ...
    wing_geom2(root_chord,wing_length,camber,twist_in,y_stations,LE_knots,TE_knots,N_gauss,N_stations,connectivity)

%% basic geometry:

y = (0:wing_length/50:wing_length);

LE = spline(y_stations,LE_knots,y); LE_pp = spline(y_stations,LE_knots);
TE = spline(y_stations,TE_knots,y); TE_pp = spline(y_stations,TE_knots);

chord = TE-LE;

x = repmat([0:1/(round(length(y)/(wing_length/root_chord))-1):1]',1,length(y));
x = x.*repmat(chord,round(length(y)/(wing_length/root_chord)),1) + repmat(LE,round(length(y)/(wing_length/root_chord)),1);
xb = repmat([-.5:1/(round(length(y)/(wing_length/root_chord))-1):.5]',1,length(y)).*repmat(chord,round(length(y)/(wing_length/root_chord)),1);

y = repmat(y,round(length(y)/(wing_length/root_chord)),1);
z = ((-4*camber)./repmat(chord,round(length(y)/(wing_length/root_chord)),1)).*xb.^2 + repmat(chord,round(length(y)/(wing_length/root_chord)),1).*camber;
twist = twist_in*sin(pi*y/wing_length/2);
z = -xb.*sind(twist) + z.*cosd(twist);

%% gauss information:

t_lob = -1:.001:1;
Leg = legendre(N_gauss,t_lob); Leg = Leg(1,:);
Leg = polyfit(t_lob,Leg,N_gauss);
points = sort(roots(Leg));
weight = pi./((1-points.^2).*polyval(Leg(1:end-1).*[N_gauss:-1:1],points).^2);
points = points*pi/2 + pi/2;

y_gauss = 0:wing_length/N_stations:wing_length; y_gauss = .5*(y_gauss(1:end-1)+y_gauss(2:end));
LE_gauss = y_gauss*0; b_gauss = y_gauss*0;

for i = 1:LE_pp.pieces
    LE_gauss(y_gauss <= y_stations(i+1) & y_gauss >= y_stations(i)) = ppval(LE_pp,y_gauss(y_gauss <= y_stations(i+1) & y_gauss >= y_stations(i)));
    b_gauss(y_gauss <= y_stations(i+1) & y_gauss >= y_stations(i)) = (ppval(TE_pp,y_gauss(y_gauss <= y_stations(i+1) & y_gauss >= y_stations(i))) - ppval(LE_pp,y_gauss(y_gauss <= y_stations(i+1) & y_gauss >= y_stations(i))))/2;
end

x_gauss = cos(points); 
x_gauss = reshape(repmat(x_gauss,1,N_stations).*repmat(b_gauss,N_gauss,1),1,N_gauss*N_stations)';

z_gauss = ((-2*camber)./reshape(repmat(b_gauss,N_gauss,1),1,N_gauss*N_stations)').*x_gauss.^2 + 2*camber*reshape(repmat(b_gauss,N_gauss,1),1,N_gauss*N_stations)';
twist = reshape(repmat(twist_in*sin(pi*y_gauss/wing_length/2),N_gauss,1),1,N_gauss*N_stations)';
z_gauss = -x_gauss.*sind(twist) + z_gauss.*cosd(twist);
z_gauss = -z_gauss;

%% remap the cell vertices

X_nodes = (-root_chord/2)*connectivity.vertices.coords(:,2) + root_chord/2;
Y_nodes = (wing_length/2)*connectivity.vertices.coords(:,1) + wing_length/2;

for i = 1:LE_pp.pieces
    X_nodes(Y_nodes < y_stations(i+1) & Y_nodes  >= y_stations(i)) = ((ppval(TE_pp,Y_nodes(Y_nodes < y_stations(i+1) & Y_nodes  >= y_stations(i))) - ppval(LE_pp,Y_nodes(Y_nodes < y_stations(i+1) & Y_nodes  >= y_stations(i))))/root_chord).*X_nodes(Y_nodes < y_stations(i+1) & Y_nodes  >= y_stations(i)) + ppval(LE_pp,Y_nodes(Y_nodes < y_stations(i+1) & Y_nodes  >= y_stations(i)));
end
X_nodes(Y_nodes == y_stations(end)) = ((ppval(TE_pp,Y_nodes(Y_nodes == y_stations(end))) - ppval(LE_pp,Y_nodes(Y_nodes == y_stations(end))))/root_chord).*X_nodes(Y_nodes == y_stations(end)) + ppval(LE_pp,Y_nodes(Y_nodes == y_stations(end)));
connectivity.vertices.coords_new = [X_nodes,Y_nodes];

%% create control points along TE

line_y = 0:connectivity.gage:wing_length;
if line_y(end) ~= wing_length
    line_y(end+1) = wing_length;
end
line_x = 0*line_y;
for i = 1:TE_pp.pieces
    line_x(line_y <= y_stations(i+1) & line_y  >= y_stations(i)) = ppval(TE_pp,line_y(line_y <= y_stations(i+1) & line_y  >= y_stations(i)));
end
add = connectivity.vertices.coords_new(connectivity.vertices.coords(:,2) == -1,:); no = find(connectivity.vertices.coords(:,2) == -1);
TE = [[line_x',line_y'];add]; ID = [0*line_x';0*add(:,1)+1]; no = [0*line_x';no];
[i,i] = sort(TE(:,2)); TE = TE(i,:); ID = ID(i); no = no(i);
foo = diff(TE(:,2));  TE(foo == 0,:) = []; ID(foo == 0,:) = []; ID(1) = 1; ID(end) = 1; no(foo == 0,:) = [];
dist = min(TE(3:end,2) - TE(2:end-1,2),TE(2:end-1,2) - TE(1:end-2,2)); dist = [connectivity.gage;dist;connectivity.gage];
dist(ID == 1) = connectivity.gage;
TE(dist < connectivity.gage/2,:) = []; no(dist < connectivity.gage/2,:) = [];
no(no == 0) = connectivity.vertices.Nvertices+1:connectivity.vertices.Nvertices+sum(no == 0);
perim = TE; perim_no = no;

TE_pp_linear = interp1(TE(:,2),TE(:,1),'linear','pp');

%% create control points along tip

tip = [TE_knots(end):-connectivity.gage:LE_knots(end)]'; 
if tip(end) ~= LE_knots(end)
    tip(end+1) = LE_knots(end);
end
tip(:,2) = y_stations(end);
add = connectivity.vertices.coords_new(connectivity.vertices.coords(:,1) == 1,:); ID = [0*tip(:,1);0*add(:,1)+1]; no = find(connectivity.vertices.coords(:,1) == 1);
no = [0*tip(:,1);no]; tip = [tip;add]; 
[i,i] = sort(tip(:,1)); tip = flipud(tip(i,:)); ID = flipud(ID(i,:)); no = flipud(no(i));
foo = diff(tip(:,1)); foo = [1;foo]; tip(foo == 0,:) = []; ID(foo == 0,:) = []; ID(1) = 1; ID(end) = 1; no(foo == 0,:) = [];
dist = min(-tip(3:end,1) + tip(2:end-1,1),-tip(2:end-1,1) + tip(1:end-2,1)); dist = [connectivity.gage;dist;connectivity.gage];
dist(ID == 1) = connectivity.gage;
tip(dist < connectivity.gage/2,:) = []; no(dist < connectivity.gage/2,:) = [];
no(no == 0) = max(perim_no)+1:max(perim_no)+sum(no == 0);
perim = [perim;tip(2:end,:)]; perim_no = [perim_no;no(2:end)];

%% create control points along LE

line_y = wing_length:-connectivity.gage:0;
if line_y(end) ~= 0
    line_y(end+1) = 0;
end
line_x = 0*line_y;
for i = 1:LE_pp.pieces
    line_x(line_y <= y_stations(i+1) & line_y  >= y_stations(i)) = ppval(LE_pp,line_y(line_y <= y_stations(i+1) & line_y  >= y_stations(i)));
end
add = connectivity.vertices.coords_new(connectivity.vertices.coords(:,2) == 1,:);  no = find(connectivity.vertices.coords(:,2) == 1);
LE = [[line_x',line_y'];add]; ID = [0*line_x';0*add(:,1)+1];  no = [0*line_x';no];
[i,i] = sort(LE(:,2)); LE = flipud(LE(i,:)); ID = flipud(ID(i,:)); no = flipud(no(i));
foo = diff(LE(:,2)); foo = [1;foo]; LE(foo == 0,:) = []; ID(foo == 0,:) = []; ID(1) = 1; ID(end) = 1; no(foo == 0,:) = [];
dist = min(-LE(3:end,2) + LE(2:end-1,2),-LE(2:end-1,2) + LE(1:end-2,2)); dist = [connectivity.gage;dist;connectivity.gage];
dist(ID == 1) = connectivity.gage;
LE(dist < connectivity.gage/2,:) = []; no(dist < connectivity.gage/2,:) = [];
no(no == 0) = max(perim_no)+1:max(perim_no)+sum(no == 0);
perim = [perim;LE(2:end,:)]; perim_no = [perim_no;no(2:end)];

LE_pp_linear = interp1(LE(:,2),LE(:,1),'linear','pp');

%% create control points along root

root = [LE_knots(1):connectivity.gage:TE_knots(1)]';
if root(end) ~= TE_knots(1)
    root(end+1) = TE_knots(1);
end
root(:,2) = y_stations(1);
add = connectivity.vertices.coords_new(connectivity.vertices.coords(:,1) == -1,:); ID = [0*root(:,1);0*add(:,1)+1]; no = find(connectivity.vertices.coords(:,1) == -1);
no = [0*root(:,1);no]; root = [root;add];
[i,i] = sort(root(:,1)); root = root(i,:); ID = ID(i,:); no = no(i);
foo = diff(root(:,1));  root(foo == 0,:) = []; ID(foo == 0,:) = []; ID(1) = 1; ID(end) = 1; no(foo == 0,:) = [];
dist = min(root(3:end,1) - root(2:end-1,1),root(2:end-1,1) - root(1:end-2,1)); dist = [connectivity.gage;dist;connectivity.gage];
dist(ID == 1) = connectivity.gage;
root(dist < connectivity.gage/2,:) = []; no(dist < connectivity.gage/2,:) = [];
no(no == 0) = max(perim_no)+1:max(perim_no)+sum(no == 0);
perim = [perim;root(2:end-1,:)]; perim_no = [perim_no;no(2:end-1)];

%% define edges and 2 nodes for each

interior = find(connectivity.vertices.type == 'I')'; N_interior = sum(connectivity.vertices.type == 'I');
perim_no = [perim_no;interior];
perim = [perim;connectivity.vertices.coords_new(interior,:)];

[ID,ID] = sort(perim_no);
node = perim(ID,:);
edge = [perim_no(1:end-N_interior),[perim_no(2:end-N_interior);perim_no(1)]];

%% for each edge, identify which edge of the original cell it belongs to

a = 0; b = 0; count = 1;
edge_ID = zeros(length(edge),1);
for i = 1:length(edge)
   if edge(i,1) <= connectivity.vertices.Nvertices
       a = edge(i,1);
   end
   if a > 0 && edge(i,2) <= connectivity.vertices.Nvertices
       b = edge(i,2);
   end
   if a > 0 && b > 0
       cell_edge = union(intersect(find([a] == connectivity.edges.vertices(:,1)),find([b] == connectivity.edges.vertices(:,2))),intersect(find([b] == connectivity.edges.vertices(:,1)),find([a] == connectivity.edges.vertices(:,2))));
       edge_ID(count:i) = cell_edge;
       count = i + 1;
       a = 0; b = 0;
   end
end

%% Create edges along internal boundaries

foo = connectivity.edges.vertices(connectivity.edges.type == 'I',:);
foo2 = 1:connectivity.edges.Nedges; foo2 = foo2(connectivity.edges.type == 'I');
for i = 1:size(foo,1)
    X1 = connectivity.vertices.coords_new(foo(i,1),1);
    Y1 = connectivity.vertices.coords_new(foo(i,1),2);
    X2 = connectivity.vertices.coords_new(foo(i,2),1);
    Y2 = connectivity.vertices.coords_new(foo(i,2),2);
    L = sqrt((X2-X1)^2+(Y2-Y1)^2);
    s = [0:L/max(round(L/connectivity.gage),1):L]'; 
    
    new = max(max(edge))+1:max(max(edge))+length(s)-2;
    node = [node;[X1+(X2-X1)*s(2:end-1)/L,Y1+(Y2-Y1)*s(2:end-1)/L]];
    edge = [edge;[[foo(i,1);new'],[new';foo(i,2)]]];
    edge_ID = [edge_ID;foo2(i)*ones(length(s)-1,1)];    
end

%% Define the faces based on the original cells:

for i = 1:connectivity.cells.Ncells
    foo = sum(repmat(edge_ID,1,length(cell2mat(connectivity.cells.edges(i)))) == repmat(cell2mat(connectivity.cells.edges(i)),length(edge_ID),1),2);
    face{i} = find(foo)';
end

%% Create mesh:

hdata.hmax  = 1;
options.maxit = 20;
options.dhmax = 20;
options.output = false;

if (~isdeployed) 
    p = strcat(pwd, '/mesh2d');
    path(p, path);
    %path('./mesh2d',path); 
end

[trinodes,trielements,fnum] = meshfaces(node,edge,face,hdata,options);
x_FEA = trinodes(:,1); y_FEA = trinodes(:,2);

%% define camber of mesh:

b_FEA = (ppval(TE_pp,y_FEA) - ppval(LE_pp,y_FEA))/2;
LE_FEA = ppval(LE_pp,y_FEA); x_FEA2 = x_FEA-LE_FEA-b_FEA;

z_FEA = (-2*camber./b_FEA).*x_FEA2.^2 + 2*camber*b_FEA;

twist = twist_in*sin(pi*y_FEA/wing_length/2);
z_FEA = -x_FEA2.*sind(twist) + z_FEA.*cosd(twist);

%% outward normal of mesh

% b_y = TE_pp; b_y.coefs = .5*(TE_pp.coefs - LE_pp.coefs)*diag(3:-1:1,1);
% b_y = ppval(b_y,y_FEA);
% tw_y = twist_in*cos(pi*y_FEA/wing_length/2)*pi/wing_length/2;
% zc_x = (-4*camber./b_FEA).*x_FEA2;
% zc_y = (2*camber./b_FEA./b_FEA).*b_y.*x_FEA2.^2 + 2*camber*b_y;
% 
% o1 = sind(twist) - zc_x.*cosd(twist);
% o2 = x_FEA2.*cosd(twist).*tw_y*pi/180 - zc_y.*cosd(twist) + z_FEA.*sind(twist).*tw_y*pi/180;
% o3 = o2*0 + 1;
% 
% outward = [o1,o2,o3]./repmat(sqrt(o1.^2 + o2.^2 + o3.^2),1,3);

%% define the beam elements along the interior:

beam_el = []; beam_ID = [];
tol = 1E-10;

Nmax = 0;
temp = zeros(connectivity.cells.Ncells,100);
for i = 1:connectivity.cells.Ncells
   N = length(cell2mat(connectivity.cells.edges(i)));
   temp(i,1:N) = cell2mat(connectivity.cells.edges(i));
   if N > Nmax
       Nmax = N;
   end
end
temp = temp(:,1:Nmax);

foo = 1:connectivity.edges.Nedges; foo = foo(connectivity.edges.type == 'I');
for i = 1:length(foo)
    cand = sum(temp == foo(i),2);
    face = 1:connectivity.cells.Ncells; face = face(cand == 1); face = face(1);
    
    tri2 = trielements(fnum == face,:);
    cand2 = unique(tri2(:));
    
    X1 = connectivity.vertices.coords_new(connectivity.edges.vertices(foo(i),1),1);
    X2 = connectivity.vertices.coords_new(connectivity.edges.vertices(foo(i),2),1);
    Y1 = connectivity.vertices.coords_new(connectivity.edges.vertices(foo(i),1),2);
    Y2 = connectivity.vertices.coords_new(connectivity.edges.vertices(foo(i),2),2);
    
    theta_desired = atan2(Y2-Y1,X2-X1);
    theta = atan2(y_FEA(cand2)-Y1,x_FEA(cand2)-X1);
    cand3 = cand2(theta-theta_desired < tol & -theta+theta_desired < tol);
    cand3 = [cand3;cand2(abs(x_FEA(cand2)-X1) < tol & abs(y_FEA(cand2)-Y1) < tol)];
    cand3 = unique(cand3);   % if theta_desired is 0, it will add in X1 twice
    
    r = sqrt((x_FEA(cand3) - X1).^2 + (y_FEA(cand3) - Y1).^2); [j,j] = sort(r); cand3 = cand3(j);
    
    beam_el = [beam_el;[cand3(1:end-1),cand3(2:end)]];
    beam_ID = [beam_ID;repmat('I',length(cand3)-1,1)];
    
end

%% define the beam elements along the root:

foo = 1:connectivity.edges.Nedges; 
foo = foo(connectivity.vertices.coords(connectivity.edges.vertices(:,1),1) == -1 & connectivity.vertices.coords(connectivity.edges.vertices(:,2),1) == -1);
for i = 1:length(foo)
    cand = sum(temp == foo(i),2);
    face = 1:connectivity.cells.Ncells; face = face(cand == 1); face = face(1);
    
    tri2 = trielements(fnum == face,:);
    cand2 = unique(tri2(:));
    
    cand3 = cand2(y_FEA(cand2) == 0);[j,j] = sort(x_FEA(cand3)); cand3 = cand3(j);
    
    beam_el = [beam_el;[cand3(1:end-1),cand3(2:end)]];
    beam_ID = [beam_ID;repmat('W',length(cand3)-1,1)];
end

%% define the beam elements along the TE:
foo = 1:connectivity.edges.Nedges; 
foo = foo(connectivity.vertices.coords(connectivity.edges.vertices(:,1),2) == -1 & connectivity.vertices.coords(connectivity.edges.vertices(:,2),2) == -1);
for i = 1:length(foo)
    cand = sum(temp == foo(i),2);
    face = 1:connectivity.cells.Ncells; face = face(cand == 1); face = face(1);
    
    tri2 = trielements(fnum == face,:);
    cand2 = unique(tri2(:));
    
    cand3 = cand2(abs(x_FEA(cand2) - ppval(TE_pp_linear,y_FEA(cand2))) < tol);
    [j,j] = sort(y_FEA(cand3)); cand3 = cand3(j);
    
    beam_el = [beam_el;[cand3(1:end-1),cand3(2:end)]];
    beam_ID = [beam_ID;repmat('S',length(cand3)-1,1)];
end

%% define the beam elements along the tip:

foo = 1:connectivity.edges.Nedges; 
foo = foo(connectivity.vertices.coords(connectivity.edges.vertices(:,1),1) == 1 & connectivity.vertices.coords(connectivity.edges.vertices(:,2),1) == 1);
for i = 1:length(foo)
    cand = sum(temp == foo(i),2);
    face = 1:connectivity.cells.Ncells; face = face(cand == 1); face = face(1);
    
    tri2 = trielements(fnum == face,:);
    cand2 = unique(tri2(:));
    
    cand3 = cand2(y_FEA(cand2) == wing_length);[j,j] = sort(x_FEA(cand3)); cand3 = cand3(j);
    
    beam_el = [beam_el;[cand3(1:end-1),cand3(2:end)]];
    beam_ID = [beam_ID;repmat('E',length(cand3)-1,1)];
end

%% define the beam elements along the LE:
foo = 1:connectivity.edges.Nedges; 
foo = foo(connectivity.vertices.coords(connectivity.edges.vertices(:,1),2) == 1 & connectivity.vertices.coords(connectivity.edges.vertices(:,2),2) == 1);
for i = 1:length(foo)
    cand = sum(temp == foo(i),2);
    face = 1:connectivity.cells.Ncells; face = face(cand == 1); face = face(1);
    
    tri2 = trielements(fnum == face,:);
    cand2 = unique(tri2(:));
    
    cand3 = cand2(abs(x_FEA(cand2) - ppval(LE_pp_linear,y_FEA(cand2))) < tol);
    [j,j] = sort(y_FEA(cand3)); cand3 = cand3(j);
    
    beam_el = [beam_el;[cand3(1:end-1),cand3(2:end)]];
    beam_ID = [beam_ID;repmat('N',length(cand3)-1,1)];
end