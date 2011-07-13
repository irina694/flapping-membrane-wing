function[K,F_dP] = stiffness_matrix(X,Y,Z,trielements,beam_elements,membrane_asmb,beam_asmb,beam_ID,GlobalDOF,fdof,mat)

K = zeros(18*18*length(trielements) + 12*12*length(beam_elements),1);
ID = zeros(18*18*length(trielements) + 12*12*length(beam_elements),2);

F_dP = zeros(18*1*length(trielements),1);
ID2 = zeros(18*1*length(trielements),2);

%% Membrane elements:

C_membrane = (mat.membrane.E/(1-mat.membrane.poisson^2))*[1,mat.membrane.poisson,0;mat.membrane.poisson,1,0;0,0,(1-mat.membrane.poisson)/2];

for i = 1:length(trielements)
    
    %% Undeformed element nodes in global system:
    
    Xg1 = [X(trielements(i,1)); Y(trielements(i,1)); Z(trielements(i,1))];
    Xg2 = [X(trielements(i,2)); Y(trielements(i,2)); Z(trielements(i,2))];
    Xg3 = [X(trielements(i,3)); Y(trielements(i,3)); Z(trielements(i,3))];
    
    %% Transformation matrix:
    
    v1 = (Xg2-Xg1)/sqrt((Xg2(1)-Xg1(1))^2 + (Xg2(2)-Xg1(2))^2 + (Xg2(3)-Xg1(3))^2);
    v13 = (Xg3-Xg1)/sqrt((Xg3(1)-Xg1(1))^2 + (Xg3(2)-Xg1(2))^2 + (Xg3(3)-Xg1(3))^2);
    v3 = [v1(2)*v13(3)-v1(3)*v13(2);-v1(1)*v13(3)+v1(3)*v13(1);v1(1)*v13(2)-v1(2)*v13(1)];
    v3 = v3/sqrt(v3(1)^2 + v3(2)^2 + v3(3)^2);
    v2 = [v3(2)*v1(3)-v3(3)*v1(2);-v3(1)*v1(3)+v3(3)*v1(1);v3(1)*v1(2)-v3(2)*v1(1)];
    v2 = v2/sqrt(v2(1)^2 + v2(2)^2 + v2(3)^2);
    Ro = [v1,v2,v3];
    Eo = zeros(18,18); Eo(1:3,1:3) = Ro; Eo(4:6,4:6) = Ro; Eo(7:9,7:9) = Ro; Eo(10:12,10:12) = Ro; Eo(13:15,13:15) = Ro; Eo(16:18,16:18) = Ro;
    
    %% Local system information:
    
    xg2 = Ro'*(Xg2-Xg1); xg3 = Ro'*(Xg3-Xg1);
    x2 = xg2(1); x3 = xg3(1); y3 = xg3(2);
    A = .5*x2*y3;
    
    %% Form and assemble stiffness matrix:
    
    Kel = zeros(18,18);
    
    Buv = [-1/x2,0,1/x2,0,0,0;0,(x3-x2)/(x2*y3),0,-x3/(x2*y3),0,1/y3;(x3-x2)/(x2*y3),-1/x2,-x3/(x2*y3),1/x2,1/y3,0];
    Kel([1,2,7,8,13,14],[1,2,7,8,13,14]) = A*mat.membrane.thickness*Buv'*C_membrane*Buv;
    
    Bw = [-1/x2,1/x2,0;(x3-x2)/(x2*y3),-x3/(x2*y3),1/y3];
    Kel([3,9,15],[3,9,15]) = A*mat.membrane.pre_stress*Bw'*Bw;
    
    Kel = Eo*Kel*Eo';
    
    temp = membrane_asmb(i,:)'; temp = temp(:,[ones(18,1)]); temp2 = temp';
    ID((i-1)*18*18+1:i*18*18,:) = [temp(:),temp2(:)];
    K((i-1)*18*18+1:i*18*18,:) = Kel(:);
    
    %% Pressure-to-force computation:
    
    ID2((i-1)*18*1+1:i*18*1,:) = [membrane_asmb(i,:)',i*ones(18,1)];
    F_dP((i-1)*18*1+1:i*18*1,:) = (A/3)*Eo*[0;0;1;0;0;0;0;0;1;0;0;0;0;0;1;0;0;0];
    
end

%% Beam elements

EA = zeros(length(beam_elements),1);
EA(beam_ID == 'N') = mat.LE.thickness*mat.LE.base*mat.LE.E;
EA(beam_ID == 'W') = mat.ROOT.thickness*mat.ROOT.base*mat.ROOT.E;
EA(beam_ID == 'S') = mat.TE.thickness*mat.TE.base*mat.TE.E;
EA(beam_ID == 'E') = mat.TIP.thickness*mat.TIP.base*mat.TIP.E;
EA(beam_ID == 'I') = mat.INT.thickness*mat.INT.base*mat.INT.E;

EIy = zeros(length(beam_elements),1);
EIy(beam_ID == 'N') = (1/12)*(mat.LE.thickness^3)*mat.LE.base*mat.LE.E;
EIy(beam_ID == 'W') = (1/12)*(mat.ROOT.thickness^3)*mat.ROOT.base*mat.ROOT.E;
EIy(beam_ID == 'S') = (1/12)*(mat.TE.thickness^3)*mat.TE.base*mat.TE.E;
EIy(beam_ID == 'E') = (1/12)*(mat.TIP.thickness^3)*mat.TIP.base*mat.TIP.E;
EIy(beam_ID == 'I') = (1/12)*(mat.INT.thickness^3)*mat.INT.base*mat.INT.E;

EIz = zeros(length(beam_elements),1);
EIz(beam_ID == 'N') = (1/12)*(mat.LE.base^3)*mat.LE.thickness*mat.LE.E;
EIz(beam_ID == 'W') = (1/12)*(mat.ROOT.base^3)*mat.ROOT.thickness*mat.ROOT.E;
EIz(beam_ID == 'S') = (1/12)*(mat.TE.base^3)*mat.TE.thickness*mat.TE.E;
EIz(beam_ID == 'E') = (1/12)*(mat.TIP.base^3)*mat.TIP.thickness*mat.TIP.E;
EIz(beam_ID == 'I') = (1/12)*(mat.INT.base^3)*mat.INT.thickness*mat.INT.E;

GJ = zeros(length(beam_elements),1);
GJ(beam_ID == 'N') = (1/3)*(mat.LE.thickness^3)*mat.LE.base*mat.LE.E/2.6;
GJ(beam_ID == 'W') = (1/3)*(mat.ROOT.thickness^3)*mat.ROOT.base*mat.ROOT.E/2.6;
GJ(beam_ID == 'S') = (1/3)*(mat.TE.thickness^3)*mat.TE.base*mat.TE.E/2.6;
GJ(beam_ID == 'E') = (1/3)*(mat.TIP.thickness^3)*mat.TIP.base*mat.TIP.E/2.6;
GJ(beam_ID == 'I') = (1/3)*(mat.INT.thickness^3)*mat.INT.base*mat.INT.E/2.6;

for i = 1:length(beam_elements)

    %% Undeformed element nodes in global system:
    
    Xg1 = [X(beam_elements(i,1)); Y(beam_elements(i,1)); Z(beam_elements(i,1))];
    Xg2 = [X(beam_elements(i,2)); Y(beam_elements(i,2)); Z(beam_elements(i,2))];
    
    %% Transformation matrix for undeformed geometry:
    
    L = sqrt((Xg2(1)-Xg1(1))^2 + (Xg2(2)-Xg1(2))^2 + (Xg2(3)-Xg1(3))^2);
    e1 = (Xg2-Xg1)/L;
    e2 = [-e1(2)/(e1(1)+1E-10),1,0]'/sqrt(1 + (-e1(2)/(e1(1)+1E-10))^2); % this assumes wing is flat
    e3 = [e1(2)*e2(3)-e1(3)*e2(2);e1(3)*e2(1)-e1(1)*e2(3);e1(1)*e2(2)-e1(2)*e2(1)]/(e1(2)^2*e2(3)^2-2*e1(2)*e2(3)*e1(3)*e2(2)+e1(3)^2*e2(2)^2+e1(3)^2*e2(1)^2-2*e1(3)*e2(1)*e1(1)*e2(3)+e1(1)^2*e2(3)^2+e1(1)^2*e2(2)^2-2*e1(1)*e2(2)*e1(2)*e2(1)+e1(2)^2*e2(1)^2)^(1/2);
    Ro = [e1,e2,e3];
    Eo = zeros(12,12); Eo(1:3,1:3)=Ro; Eo(4:6,4:6)=Ro; Eo(7:9,7:9)=Ro; Eo(10:12,10:12)=Ro;
    
    %% Form and assemble stiffness matrix:
    
    Kel = zeros(12,12);
    Kel([1,7],[1,7]) = (EA(i)/L)*[1,-1;-1,1];
    Kel([4,10],[4,10]) = (GJ(i)/L)*[1,-1;-1,1];
    Kel([2,6,8,12],[2,6,8,12]) = [12*EIz(i)/L^3,6*EIz(i)/L^2,-12*EIz(i)/L^3,6*EIz(i)/L^2;6*EIz(i)/L^2,4*EIz(i)/L,-6*EIz(i)/L^2,2*EIz(i)/L;
        -12*EIz(i)/L^3,-6*EIz(i)/L^2,12*EIz(i)/L^3,-6*EIz(i)/L^2;6*EIz(i)/L^2,2*EIz(i)/L,-6*EIz(i)/L^2,4*EIz(i)/L];
    Kel([3,5,9,11],[3,5,9,11]) = [12*EIy(i)/L^3,-6*EIy(i)/L^2,-12*EIy(i)/L^3,-6*EIy(i)/L^2;-6*EIy(i)/L^2,4*EIy(i)/L,6*EIy(i)/L^2,2*EIy(i)/L;
        -12*EIy(i)/L^3,6*EIy(i)/L^2,12*EIy(i)/L^3,6*EIy(i)/L^2;-6*EIy(i)/L^2,2*EIy(i)/L,6*EIy(i)/L^2,4*EIy(i)/L];
    
    Kel = Eo*Kel*Eo'; 
    
    temp = beam_asmb(i,:)'; temp = temp(:,[ones(12,1)]); temp2 = temp';
    ID([(i-1)*12*12+1:i*12*12]+18*18*length(trielements),:) = [temp(:),temp2(:)];
    K([(i-1)*12*12+1:i*12*12]+18*18*length(trielements),:) = Kel(:);
    
end

K = sparse(ID(:,1),ID(:,2),K,GlobalDOF,GlobalDOF);
K = K(fdof,fdof);

F_dP = sparse(ID2(:,1),ID2(:,2),F_dP,GlobalDOF,length(trielements));
F_dP = F_dP(fdof,:);
