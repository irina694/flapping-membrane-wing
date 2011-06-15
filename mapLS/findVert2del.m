function [vert2del,fstV,lstV] = findVert2del(x,Nvert)

    vert2del=[];
    fstV = [];
    lstV = [];
    n=length(x);
    if(length(x)<=2)
        return;
    elseif (n == Nvert)
        vert2del = x(2:n-1);
        fstV = x(1);
        lstV = x(n); 
        return;        
    end

    foundFst = false;
    for i=2:n-1
        if(x(i)-x(i-1)==1 & x(i+1)-x(i)==1)       %delete
            vert2del=[vert2del,x(i)];
        elseif (~foundFst)
            fstV = x(i);
            foundFst = true;
        elseif (foundFst)
            lstV = x(i);
        end
    end
    
    if (isempty(fstV) & isempty(lstV))
        fstV = x(1);
        lstV = x(n);
    end
        
    if(x(1)==1 & x(n)==Nvert)        %cicle
        if(x(n-1)==Nvert-1)
            if (isempty(lstV)), lstV=1;
            else aux = fstV; fstV = lstV; lstV = aux; end
            vert2del=[vert2del,Nvert];            
        end
        if(x(2)==2)
            if (isempty(lstV)), lstV = fstV; fstV=Nvert;
            elseif(x(n-1)~=Nvert-1) aux = fstV; fstV = lstV; lstV = aux; end
            vert2del=[1,vert2del];
        end        
    end          
end