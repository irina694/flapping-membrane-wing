function[u_convert] = disp_interp(x,y,b,LE,N_gauss,N_stations,x_FEA,y_FEA,trielements,bc,fdof)

x = x+reshape(repmat(b,N_gauss,1),1,N_gauss*N_stations)'+reshape(repmat(LE,N_gauss,1),1,N_gauss*N_stations)';
y = reshape(repmat(y,N_gauss,1),1,N_gauss*N_stations)';

T = tsearch(x_FEA,y_FEA,trielements,x,y);

%% Some gauss points lie outside the triangulation:

lost = find(isnan(T));

Xm = mean(x_FEA(trielements),2); Ym = mean(y_FEA(trielements),2);

foo = (repmat(x(lost),1,length(trielements))-repmat(Xm',length(lost),1)).^2 + (repmat(y(lost),1,length(trielements))-repmat(Ym',length(lost),1)).^2;
[i,i] = min(foo,[],2);
T(lost) = i;

%% Use shape functions to interpolate geometry:

x1 = x_FEA(trielements(T,1)); x2 = x_FEA(trielements(T,2)); x3 = x_FEA(trielements(T,3));
y1 = y_FEA(trielements(T,1)); y2 = y_FEA(trielements(T,2)); y3 = y_FEA(trielements(T,3));

A2 = x1.*y2-x1.*y3+x2.*y3-x2.*y1+x3.*y1-x3.*y2;

a1 = (x2.*y3-x3.*y2)./A2; a2 = (x3.*y1-x1.*y3)./A2; a3 = (x1.*y2-x2.*y1)./A2;
b1 = (y2-y3)./A2; b2 = (y3-y1)./A2; b3 = (y1-y2)./A2;
c1 = (x3-x2)./A2; c2 = (x1-x3)./A2; c3 = (x2-x1)./A2;

N1 = a1 + b1.*x + c1.*y; N2 = a2 + b2.*x + c2.*y; N3 = a3 + b3.*x + c3.*y; 

Q1 = sparse([(1:N_stations*N_gauss)';(1:N_stations*N_gauss)';(1:N_stations*N_gauss)'],[trielements(T,1);trielements(T,2);trielements(T,3)],[N1;N2;N3],N_stations*N_gauss,length(x_FEA));
Q2 = sparse(1:length(x_FEA),3:6:length(x_FEA)*6-3,ones(length(x_FEA),1),length(x_FEA),6*length(x_FEA));

Q3 = sparse(6*length(x_FEA),6*length(x_FEA)-length(bc));
Q3(fdof,:) = diag(ones(6*length(x_FEA)-length(bc),1));

u_convert = -Q1*Q2*Q3;
