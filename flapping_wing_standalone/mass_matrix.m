function[M] = mass_matrix(X,Y,Z,trielements,beam_elements,membrane_asmb,beam_asmb,beam_ID,GlobalDOF,fdof,mat)

M = zeros(18*18*length(trielements) + 12*12*length(beam_elements),1);
ID = zeros(18*18*length(trielements) + 12*12*length(beam_elements),2);

%% Membrane elements:

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
    
    %% Form and assemble mass matrix:
    
    Mel = zeros(18,18);
    Mel([1,2,3,7,8,9,13,14,15],[1,2,3,7,8,9,13,14,15]) = (mat.membrane.rho*mat.membrane.thickness*A/12)*[2*eye(3),eye(3),eye(3);eye(3),2*eye(3),eye(3);eye(3),eye(3),2*eye(3)];
    Mel = Eo*Mel*Eo';
    
    temp = membrane_asmb(i,:)'; temp = temp(:,[ones(18,1)]); temp2 = temp';
    ID((i-1)*18*18+1:i*18*18,:) = [temp(:),temp2(:)];
    M((i-1)*18*18+1:i*18*18,:) = Mel(:);

end

%% Beam elements

pA = zeros(length(beam_elements),1);
pA(beam_ID == 'N') = mat.LE.rho*mat.LE.thickness*mat.LE.base;
pA(beam_ID == 'W') = mat.ROOT.rho*mat.ROOT.thickness*mat.ROOT.base;
pA(beam_ID == 'S') = mat.TE.rho*mat.TE.thickness*mat.TE.base;
pA(beam_ID == 'E') = mat.TIP.rho*mat.TIP.thickness*mat.TIP.base;
pA(beam_ID == 'I') = mat.INT.rho*mat.INT.thickness*mat.INT.base;

Io = zeros(length(beam_elements),1);
Io(beam_ID == 'N') = mat.LE.rho*mat.LE.thickness*mat.LE.base*(mat.LE.thickness^2+mat.LE.base^2)/12;
Io(beam_ID == 'W') = mat.ROOT.rho*mat.ROOT.thickness*mat.ROOT.base*(mat.ROOT.thickness^2+mat.ROOT.base^2)/12;
Io(beam_ID == 'S') = mat.TE.rho*mat.TE.thickness*mat.TE.base*(mat.TE.thickness^2+mat.TE.base^2)/12;
Io(beam_ID == 'E') = mat.TIP.rho*mat.TIP.thickness*mat.TIP.base*(mat.TIP.thickness^2+mat.TIP.base^2)/12;
Io(beam_ID == 'I') = mat.INT.rho*mat.INT.thickness*mat.INT.base*(mat.INT.thickness^2+mat.INT.base^2)/12;

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
    
    %% Form and assemble mass matrix:
    
    Mel = zeros(12,12);
    Mel([1,7],[1,7]) = (L*pA(i)/6)*[2,1;1,2];
    Mel([4,10],[4,10]) = (L*Io(i)/6)*[2,1;1,2];
    Mel([2,6,8,12],[2,6,8,12]) = (L*pA(i)/420)*[156,22*L,54,-13*L;22*L,4*L^2,13*L,-3*L^2;54,13*L,156,-22*L;-13*L,-3*L^2,-22*L,4*L^2];
    Mel([3,5,9,11],[3,5,9,11]) = (L*pA(i)/420)*[156,-22*L,54,13*L;-22*L,4*L^2,-13*L,-3*L^2;54,-13*L,156,22*L;13*L,-3*L^2,22*L,4*L^2];
    
    Mel = Eo*(Mel)*Eo';
    
    temp = beam_asmb(i,:)'; temp = temp(:,[ones(12,1)]); temp2 = temp';
    ID([(i-1)*12*12+1:i*12*12]+18*18*length(trielements),:) = [temp(:),temp2(:)];
    M([(i-1)*12*12+1:i*12*12]+18*18*length(trielements),:) = Mel(:);
    
end

M = sparse(ID(:,1),ID(:,2),M,GlobalDOF,GlobalDOF);
M = M(fdof,fdof);

