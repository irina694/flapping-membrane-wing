function [smap,direction]=extractsSmap(X)
    global Gvars    
    
    nX=length(X);
    smap = []; 
    direction = []; 
    
    for i=1:nX
        X0=mod(abs(X(i)),Gvars.Nletters);
        smap=[smap,char(floor(X0)+65)];
        X0=X0-floor(X0);
        if(X0<0.5)
            direction=[direction,'-'];
        else
            direction=[direction,'+'];
        end
    end
    %smap(4)=smap(2);
    %direction(4)=direction(2);
    %smap(3)=smap(1);
    %direction(3)=direction(1);

end
