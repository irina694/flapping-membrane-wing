function good=checkMinAngle(C)
    
    good=true;
    
    ncell=C.Ncell;
    
    for i=1:ncell
        cell=C.cell{i};
        nvert=length(cell);
        
        x=C.coord(cell,1)';
        y=C.coord(cell,2)';         
        good=checkAngles(x,y);
        if(~good)
            return;
        end
    end


end


function good=checkAngles(x,y);
    
    global Gvars
    
    MinAngle=Gvars.MinAngle;        %angle in degrees
    good=true;
    N=length(x);
    
    for i=1:N
        if(i==1)
            e1=[x(N)-x(i),y(N)-y(i)]/sqrt((x(N)-x(i))^2+(y(N)-y(i))^2);
            e2=[x(i+1)-x(i),y(i+1)-y(i)]/sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
        elseif(i==N)
            e1=[x(i-1)-x(i),y(i-1)-y(i)]/sqrt((x(i-1)-x(i))^2+(y(i-1)-y(i))^2);
            e2=[x(1)-x(i),y(1)-y(i)]/sqrt((x(1)-x(i))^2+(y(1)-y(i))^2);
        else
            e1=[x(i-1)-x(i),y(i-1)-y(i)]/sqrt((x(i-1)-x(i))^2+(y(i-1)-y(i))^2);
            e2=[x(i+1)-x(i),y(i+1)-y(i)]/sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
        end
        theta1=atan2(e1(2),e1(1));
        theta2=atan2(e2(2),e2(1));
        theta=abs(theta1-theta2)*180/pi;

        if(theta<MinAngle | (360-theta)<MinAngle)
            good = false;
            return;            
        end            
    end
end
