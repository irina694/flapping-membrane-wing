function noOverlapping=chechOverlapping(C)
    
    noOverlapping=true;
    
    ncell=C.Ncell;
    
    for i=1:ncell
        cell=C.cell{i};
        nvert=length(cell);
        
        x=C.coord(cell,1)';
        y=C.coord(cell,2)';         
        noOverlapping=checkwindingNumber(x,y);
        if(~noOverlapping)
            return;
        end
    end


end