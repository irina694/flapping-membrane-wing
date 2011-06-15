function edgStruct = setEdgeAreaProperties(NoEdges,diam,edgStruct,ratio)

    A0 = diam^2*pi/4;
    H0 = diam;
    I0 = pi*(diam/2)^4/4;
    J0 = pi*(diam/2)^4/2;
    
    edgStruct.A = num2cell(A0*(ratio.^2));
    edgStruct.heightz = num2cell(H0*ratio);
    edgStruct.heighty = num2cell(H0*ratio);
    edgStruct.Izz = num2cell(I0*(ratio.^4));
    edgStruct.Iyy = num2cell(I0*(ratio.^4));
    edgStruct.J = num2cell(J0*(ratio.^4));
    
    Fz = cell(1,NoEdges);
    for i=1:NoEdges
        Fz{i} = ['pressure*diam*',num2str(ratio(i),6)];        
    end
    edgStruct.Fz = Fz;
    
end