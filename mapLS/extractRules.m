function [rules,direct,thick,elast,press]=extractRules(x)
    global Gvars
    
    Nrules=Gvars.Nletters;
    NLetterperRule=Gvars.NLetterperRule;
        
    
    for i=1:Nrules
        rules{i}=[];
        direct{i}=[];
        thick{i}=[];
        elast{i}=[];
        press{i}=[];
        j0=(i-1)*NLetterperRule+1;
        j1=i*NLetterperRule;
        for j=j0:j1
            jj0 = (j-1)*Gvars.Ncodes+1;
            jj1 = j*Gvars.Ncodes;
            [letter,dir,thck,elst,prss]=getletter(abs(x(jj0:jj1)),mod(j,2));
            if(letter~=' ')            
                rules{i}=[rules{i},letter];
                direct{i}=[direct{i},dir];
                thick{i}=[thick{i},thck];
                elast{i}=[elast{i},elst];
                press{i}=[press{i},prss];                
            end           
        end
        
        
        %removes the brackes from the first and the last symbols, if they
        %are bracketed
        n=length(rules{i});
        if(n>0)
            if(rules{i}(n)==']')
                rules{i}(n-3)=rules{i}(n-1);
                rules{i}(n-2:n)=[];
                direct{i}(n-3)=direct{i}(n-1);
                direct{i}(n-2:n)=[];
                thick{i}(n-3)=thick{i}(n-1);
                thick{i}(n-2:n)=[];                
                elast{i}(n-3)=elast{i}(n-1);
                elast{i}(n-2:n)=[];                
                press{i}(n-3)=press{i}(n-1);
                press{i}(n-2:n)=[];                
            end
            if(rules{i}(1)=='[')
                rules{i}(1)=rules{i}(3);
                rules{i}(2:4)=[];
                direct{i}(1)=direct{i}(3);
                direct{i}(2:4)=[];
                thick{i}(1)=thick{i}(3);
                thick{i}(2:4)=[];
                elast{i}(1)=elast{i}(3);
                elast{i}(2:4)=[];
                press{i}(1)=press{i}(3);
                press{i}(2:4)=[];
            end
        end        
    end
       
end



function [y,d,t,e,p]=getletter(x,f)
    global Gvars
     
    Nsyb = Gvars.Nletters;
    
    d=[];
    y=[];
    t=[];
    e=[];
    p=[];
    
    if (x(1)-floor(x(1)) < Gvars.pNoSymb)
        y = ' ';
        return;
    end
        
    if (x(2)-floor(x(2)) < Gvars.pTermin)
        y = 'x';
    else
        x0=mod(x(3),Nsyb);
        x1=floor(x0);
        y = char(x1+65);        
    end

    if (x(4)-floor(x(4)) < Gvars.pBrackt)
        y = ['[0',y,']'];
        d = '0000';
        t = [0 0 0 0];
        e = [0 0 0 0];
        p = [0 0 0 0];
                        
        if(x(5)-floor(x(5)) < (Gvars.pDirect + (1-2*f)*Gvars.pAltern))            
            y(2) = '+';
        else
            y(2) = '-';
        end
    
        if(x(6)-floor(x(6)) < Gvars.pOrient)
            d(3) = '+';
        else
            d(3) = '-';
        end
        
        x1 = x(7)-floor(x(7));        
        if(x1 < Gvars.pThickn)
            t(3) = Gvars.ThickRt;            
        elseif(x1 > (1-Gvars.pThickn))
            t(3) = 1/Gvars.ThickRt;                                    
        else
            t(3) = 1;
        end 

        x1 = x(8)-floor(x(8));
        if(x1 < Gvars.pElastc)
            e(3) = Gvars.ElastRt;            
        elseif(x1 > (1-Gvars.pElastc))
            e(3) = 1/Gvars.ElastRt;                                    
        else
            e(3) = 1;
        end

        x1 = x(9)-floor(x(9));
        if(x1 < Gvars.pPressr)
            p(3) = Gvars.PressRt;            
        elseif(x1 > (1-Gvars.pPressr))
            p(3) = 1/Gvars.PressRt;                                    
        else
            p(3) = 1;
        end 
        
    else
        if(x(6)-floor(x(6)) < Gvars.pOrient)
            d = '+';
        else
            d = '-';
        end
        
        x1 = x(7)-floor(x(7));
        if(x1 < Gvars.pThickn)
            t = Gvars.ThickRt;
        elseif(x1 > (1-Gvars.pThickn))
            t = 1/Gvars.ThickRt;
        else
            t = 1;
        end

        x1 = x(8)-floor(x(8));
        if(x1 < Gvars.pElastc)
            e = Gvars.ElastRt;            
        elseif(x1 > (1-Gvars.pElastc))
            e = 1/Gvars.ElastRt;                                    
        else
            e = 1;
        end

        x1 = x(9)-floor(x(9));
        if(x1 < Gvars.pPressr)
            p = Gvars.PressRt;            
        elseif(x1 > (1-Gvars.pPressr))
            p = 1/Gvars.PressRt;                                    
        else
            p = 1;
        end 
        
    end
    
end