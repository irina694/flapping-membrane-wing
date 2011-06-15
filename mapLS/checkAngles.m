function good=checkAngles(x,y);

    good=true;
    N=length(x);
    
    for i=2:N-1
        e1=[x(i-1)-x(i),y(i-1)-y(i)]/sqrt((x(i-1)-x(i))^2+(y(i-1)-y(i))^2);
        e2=[x(i+1)-x(i),y(i+1)-y(i)]/sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
        theta1=atan2(e1(2),e1(1));
        theta2=atan2(e2(2),e2(1));
        theta=abs(theta1-theta2)*180/pi
    end


end