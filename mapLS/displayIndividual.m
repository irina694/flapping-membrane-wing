function fitnessvalue=displayIndividual(x)
    global Gvars 
    
    %extracts Smap
    lXsmap = Gvars.NinitEdges;
    Xsmap = x(1:lXsmap);
    C=initStruct;
    [C.cells{1}.smap,C.cells{1}.direction]=extractsSmap(Xsmap);
    
%    [Gvars.rules,Gvars.direction] = extractRules(x(lXsmap+1:end));
    [Gvars.rules,Gvars.direction] = extractRules0(x(lXsmap+1:end));
    try
        C=evolveCell(Gvars.Nlevels,C);
        drawConnectivity2(C);
        pause(0.5);
        fitnessvalue = 1;
%        fitnessvalue = ObjectiveFunctionLoad1(C);   
    catch
        fitnessvalue=1e6;
    end
end
