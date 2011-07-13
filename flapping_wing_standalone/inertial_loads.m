function[F,power] = inertial_loads(X,Y,Z,trielements,beam_elements,membrane_asmb,beam_asmb,beam_ID,GlobalDOF,fdof,mat,zpp,f,fp,fpp,th,thp,thpp)

F = zeros(GlobalDOF,length(zpp));
power = zeros(1,length(zpp));

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
    
    %% Form and assemble force vector:
    
    Fel = zeros(18,length(zpp));
    
    e_gauss = [0.5;0.5;0.0]; n_gauss = [0.0;0.5;0.5]; weight = [1/3,1/3,1/3];
    
    for j = 1:length(weight)
        e = e_gauss(j); n = n_gauss(j);
        N = [1-e-n,0,0,0,0,0,e,0,0,0,0,0,n,0,0,0,0,0;0,1-e-n,0,0,0,0,0,e,0,0,0,0,0,n,0,0,0,0;0,0,1-e-n,0,0,0,0,0,e,0,0,0,0,0,n,0,0,0];
        Fel = Fel + mat.membrane.rho*mat.membrane.thickness*A*weight(j)*N'*Ro'*[-zpp.*sin(th);zpp.*cos(th).*sin(f);zpp.*cos(f).*cos(th)] + ...
            mat.membrane.rho*mat.membrane.thickness*A*weight(j)*N'*Ro'*[(thpp.*sin(f) + 2*fp.*thp.*cos(f))*(Xg1(2) + e*Ro(2,1)*x2 + n*Ro(2,1)*x3 + n*Ro(2,2)*y3) + (thpp.*cos(f) - 2*fp.*thp.*sin(f))*(Xg1(3) + e*Ro(3,1)*x2 + n*Ro(3,1)*x3 + n*Ro(3,2)*y3) - (thp.^2.*cos(f).^2 + thp.^2.*sin(f).^2)*(Xg1(1) + e*Ro(1,1)*x2 + n*Ro(1,1)*x3 + n*Ro(1,2)*y3);-(fp.^2 + thp.^2.*sin(f).^2)*(Xg1(2) + e*Ro(2,1)*x2 + n*Ro(2,1)*x3 + n*Ro(2,2)*y3) - (cos(f).*sin(f).*thp.^2 + fpp)*(Xg1(3) + e*Ro(3,1)*x2 + n*Ro(3,1)*x3 + n*Ro(3,2)*y3) - thpp.*sin(f)*(Xg1(1) + e*Ro(1,1)*x2 + n*Ro(1,1)*x3 + n*Ro(1,2)*y3);(fpp - thp.^2.*cos(f).*sin(f))*(Xg1(2) + e*Ro(2,1)*x2 + n*Ro(2,1)*x3 + n*Ro(2,2)*y3) - (fp.^2 + thp.^2.*cos(f).^2)*(Xg1(3) + e*Ro(3,1)*x2 + n*Ro(3,1)*x3 + n*Ro(3,2)*y3) - thpp.*cos(f)*(Xg1(1) + e*Ro(1,1)*x2 + n*Ro(1,1)*x3 + n*Ro(1,2)*y3)];
    end
    
    F(membrane_asmb(i,:),:) = F(membrane_asmb(i,:),:) - Eo*Fel;
    
    %% Inertial power:
    
    vel = [Xg1(2)*thp.*sin(f) - zpp.*sin(th) + Xg1(3)*thp.*cos(f);zpp.*cos(th).*sin(f) - Xg1(1)*thp.*sin(f) - Xg1(3)*fp;Xg1(2)*fp + zpp.*cos(f).*cos(th) - Xg1(1)*thp.*cos(f);zeros(3,length(zpp));Xg2(2)*thp.*sin(f) - zpp.*sin(th) + Xg2(3)*thp.*cos(f);zpp.*cos(th).*sin(f) - Xg2(1)*thp.*sin(f) - Xg2(3)*fp;Xg2(2)*fp + zpp.*cos(f).*cos(th) - Xg2(1)*thp.*cos(f);zeros(3,length(zpp));Xg3(2)*thp.*sin(f) - zpp.*sin(th) + Xg3(3)*thp.*cos(f);zpp.*cos(th).*sin(f) - Xg3(1)*thp.*sin(f) - Xg3(3)*fp;Xg3(2)*fp + zpp.*cos(f).*cos(th) - Xg3(1)*thp.*cos(f);zeros(3,length(zpp))];
    power = power + sum((Eo*Fel).*vel);
    
 end

%  Xg12*thp*sin(fl) - zpp*sin(th) + Xg13*thp*cos(fl)
%  zpp*cos(th)*sin(fl) - Xg11*thp*sin(fl) - Xg13*flp
%  Xg12*flp + zpp*cos(fl)*cos(th) - Xg11*thp*cos(fl)
%  
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
    
    %% Form and assemble force vector:
    
    Fel = zeros(12,length(zpp));
    
    weight = [1,1]*L/2;
    gauss = [-1/sqrt(3)+1,1/sqrt(3)+1]*L/2;
    
    for q = 1:length(weight)
        N = [(L-gauss(q))/L,0,0,0,0,0,gauss(q)/L,0,0,0,0,0;...
            0,1-3*gauss(q)^2/L^2+2*gauss(q)^3/L^3,0,0,0,gauss(q)-2*gauss(q)^2/L+gauss(q)^3/L^2,0,3*gauss(q)^2/L^2-2*gauss(q)^3/L^3,0,0,0,-gauss(q)^2/L+gauss(q)^3/L^2;...
            0,0,1-3*gauss(q)^2/L^2+2*gauss(q)^3/L^3,0,-gauss(q)+2*gauss(q)^2/L-gauss(q)^3/L^2,0,0,0,3*gauss(q)^2/L^2-2*gauss(q)^3/L^3,0,gauss(q)^2/L-gauss(q)^3/L^2,0];
        N_theta = [0,0,0,(L-gauss(q))/L,0,0,0,0,0,gauss(q)/L,0,0];
        
        Fel = Fel + weight(q)*pA(i)*N'*Ro'*[-zpp.*sin(th);zpp.*cos(th).*sin(f);zpp.*cos(f).*cos(th)] + ...
            weight(q)*pA(i)*N'*Ro'* [(Xg1(2) + Ro(2,1)*gauss(q))*(thpp.*sin(f) + 2*fp.*thp.*cos(f)) + (Xg1(3) + Ro(3,1)*gauss(q))*(thpp.*cos(f) - 2*fp.*thp.*sin(f)) - (Xg1(1) + Ro(1,1)*gauss(q))*(thp.^2.*cos(f).^2 + thp.^2.*sin(f).^2);-(Xg1(2) + Ro(2,1)*gauss(q))*(fp.^2 + thp.^2.*sin(f).^2) - (Xg1(3) + Ro(3,1)*gauss(q))*(cos(f).*sin(f).*thp.^2 + fpp) - thpp.*sin(f)*(Xg1(1) + Ro(1,1)*gauss(q));(Xg1(2) + Ro(2,1)*gauss(q))*(fpp - thp.^2.*cos(f).*sin(f)) - (Xg1(3) + Ro(3,1)*gauss(q))*(fp.^2 + thp.^2.*cos(f).^2) - thpp.*cos(f)*(Xg1(1) + Ro(1,1)*gauss(q))] + ...
            weight(q)*Io(i)*N_theta'*[1,0,0]*Ro'*[fpp;thpp.*cos(f)-thp.*sin(f).*fp;-thpp.*sin(f)-thp.*cos(f).*fp];
    end
    
    F(beam_asmb(i,:),:) = F(beam_asmb(i,:),:) - Eo*Fel;
    
    %% Inertial power:
    
    vel = [Xg1(2)*thp.*sin(f) - zpp.*sin(th) + Xg1(3)*thp.*cos(f);zpp.*cos(th).*sin(f) - Xg1(1)*thp.*sin(f) - Xg1(3)*fp;Xg1(2)*fp + zpp.*cos(f).*cos(th) - Xg1(1)*thp.*cos(f);zeros(3,length(zpp));Xg2(2)*thp.*sin(f) - zpp.*sin(th) + Xg2(3)*thp.*cos(f);zpp.*cos(th).*sin(f) - Xg2(1)*thp.*sin(f) - Xg2(3)*fp;Xg2(2)*fp + zpp.*cos(f).*cos(th) - Xg2(1)*thp.*cos(f);zeros(3,length(zpp))];
    power = power + sum((Eo*Fel).*vel);
    
end

F = F(fdof,:);