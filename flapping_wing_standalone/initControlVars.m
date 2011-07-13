function [Gvars, input] = initControlVars(Gvars, input)

    

    if(nargin < 1)
        input.Nvertices = 4;
        input.vertCoord = [0 0; 1 0; 1 1; 0 1];
        input.L0 = 0;
        input.P0 = 1;
        input.K = 100;
        input.Nletters = 6;
        input.NLetterperRule = 6;
        input.Nlevels = 5;
        input.RelaxF = 0.001;
        input.ScaleFx = 0.01;
        input.Ncodes = 9;
    end        

    
    Gvars.NinitEdges = input.Nvertices;
    Gvars.Nletters = input.Nletters;
    Gvars.NLetterperRule = input.NLetterperRule;
    Gvars.Nlevels = input.Nlevels;
   
    Gvars.Ncodes = input.Ncodes;
    Gvars.pBrackt = 0.9;                   %the probability of being bracketed
    Gvars.pDirect = 0.5;                   %the probability of having a '+' direction
    Gvars.pOrient = 0.5;                   %the probability of having a '->' orientation
    Gvars.pThickn = 0.33;                  %controls the thickness rules. 0-pThickn  -> edge thins
    Gvars.pElastc = 0.33;
    Gvars.pPressr = 0.33;
    Gvars.pNoSymb = 0.2;                   %the probability of being a blank space
    Gvars.pTermin = 0.03;                  %the probability of being a terminal symbol 'x'
    Gvars.pAltern = 0.1;                   %changes the value Gvars.pDirect 
    Gvars.ThickRt = 1.0;                  
    Gvars.MaxThck = 3;                  
    Gvars.MinThck = 1/3;                  
    Gvars.ElastRt = 2.0;  
    Gvars.MaxElst = 4;                  
    Gvars.MinElst = 1/4;                  
    Gvars.PressRt = 2.0;                  
    Gvars.MaxPrss = 2;                  
    Gvars.MinPrss = 0.5;                  
    Gvars.MinAngle = 10;                   %min angle in degrees
    Gvars.MinAreaRatio = 0.3;              %the resultanting children areas cannot be to small
    Gvars.MinArea = 0.01;                  %is the cell area is less than this value doesnt split

end