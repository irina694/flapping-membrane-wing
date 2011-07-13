function[dP_convert] = pressure_interp(x,y,b,LE,N_gauss,N_stations,LE_knots,TE_knots,x_FEA,y_FEA,trielements)

x = x+reshape(repmat(b,N_gauss,1),1,N_gauss*N_stations)'+reshape(repmat(LE,N_gauss,1),1,N_gauss*N_stations)';
y = reshape(repmat(y,N_gauss,1),1,N_gauss*N_stations)';

ID = [1:length(x)]';
x_temp = reshape(x,N_gauss,N_stations); y_temp = reshape(y,N_gauss,N_stations); ID_temp = reshape(ID,N_gauss,N_stations);

x_root = ((x_temp(:,1)-x_temp(1,1))/(x_temp(end,1)-x_temp(1,1)))*(LE_knots(1)-TE_knots(1))+TE_knots(1);
x_tip = ((x_temp(:,end)-x_temp(1,end))/(x_temp(end,end)-x_temp(1,end)))*(LE_knots(end)-TE_knots(end))+TE_knots(end);

x_temp = x_temp(1:end-2,:); y_temp = y_temp(1:end-2,:); ID_temp = ID_temp(1:end-2,:); 
x_root = x_root(1:end-2); x_tip = x_tip(1:end-2);

x_temp = [x_root,(x_temp(:,2:end) + x_temp(:,1:end-1))/2,x_tip];
y_temp = [x_root*0,(y_temp(:,2:end) + y_temp(:,1:end-1))/2,x_tip*0+max(y_FEA)];

%% define quad corners

quad_ID = ID_temp(:); quad_ID = [quad_ID(1:end-1),quad_ID(2:end)]; quad_ID([N_gauss-2:N_gauss-2:(N_gauss-2)*(N_stations-1)],:) = [];

a = x_temp(1:end-1,1:end-1); b = x_temp(1:end-1,2:end); c = x_temp(2:end,2:end); d = x_temp(2:end,1:end-1);
x_quad = [a(:),b(:),c(:),d(:)];
a = y_temp(1:end-1,1:end-1); b = y_temp(1:end-1,2:end); c = y_temp(2:end,2:end); d = y_temp(2:end,1:end-1);
y_quad = [a(:),b(:),c(:),d(:)];

%% Line equations of the bottom and top quad edges:

bot_A = (x_quad(:,2)-x_quad(:,1))./(y_quad(:,2)-y_quad(:,1));
bot_B = x_quad(:,1)-y_quad(:,1).*bot_A;

top_A = (x_quad(:,3)-x_quad(:,4))./(y_quad(:,3)-y_quad(:,4));
top_B = x_quad(:,4)-y_quad(:,4).*top_A;

%% Center of each triangular FE:

Xm = mean(x_FEA(trielements),2); Ym = mean(y_FEA(trielements),2);

%% Inside?

a = repmat(Xm,1,length(quad_ID))-(repmat(bot_A',length(Xm),1).*repmat(Ym,1,length(quad_ID)) + repmat(bot_B',length(Xm),1));
a = (a < 0);

b = (repmat(top_A',length(Xm),1).*repmat(Ym,1,length(quad_ID)) + repmat(top_B',length(Xm),1)) - repmat(Xm,1,length(quad_ID));
b = (b < 0);

c = repmat(Ym,1,length(quad_ID)) - repmat(y_quad(:,1)',length(Xm),1);
c = (c > 0);

d = repmat(Ym,1,length(quad_ID)) - repmat(y_quad(:,2)',length(Xm),1);
d = (d < 0);

foo = a & b & c & d;
map = foo.*repmat(1:size(foo,2),size(foo,1),1); map = sum(map,2);

%% center of triangle falls within more than one quad?

if max(sum(foo,2)) > 1
    'interpolation problem'
    keyboard
end

%% Some of the triangles don't fall within any center...just use the closest

lost = find(sum(foo,2) == 0);
hoo = (repmat(Xm(lost),1,length(quad_ID))-repmat(mean(x_quad,2)',length(lost),1)).^2 + ...
    (repmat(Ym(lost),1,length(quad_ID))-repmat(mean(y_quad,2)',length(lost),1)).^2;
[i,i] = min(hoo,[],2);
map(lost) = i;

%% Final interpolation matrix

dP_convert = sparse([(1:length(trielements))';(1:length(trielements))'],[quad_ID(map,1);quad_ID(map,2)],0.5*ones(length(trielements)*2,1),length(trielements),N_gauss*N_stations);
