function nooverlapping=checkwindingNumber(x,y)

    nooverlapping=true;
    epsilon=1e-6;
    N=length(x);
    x1=x(1:N);
    x2=[x(2:N),x(1)];
    y1=y(1:N);
    y2=[y(2:N),y(1)];
    Lx=(x1-x2).^2;
    Ly=(y1-y2).^2;
    L=min(Lx+Ly);
    dL=sqrt(L)*epsilon;
    
    %wn=zeros(2,N);
    for i=1:N
        %computes bissectriz unity vector \\
        if(i==1)
            eB=computebissectriz([x(N),x(i),x(i+1)],[y(N),y(i),y(i+1)]);
        elseif(i==N)
            eB=computebissectriz([x(i-1),x(i),x(1)],[y(i-1),y(i),y(1)]);
        else
            eB=computebissectriz([x(i-1),x(i),x(i+1)],[y(i-1),y(i),y(i+1)]);
        end
        P1=[x(i),y(i)]+dL*eB;
        P2=[x(i),y(i)]-dL*eB; 
        wn1=windingNumber(P1,x,y);
        wn2=windingNumber(P2,x,y);   
        if(wn1+wn2~=1)
            nooverlapping=false;
            return;
        end
    end    

end


function wn=windingNumber(P,x,y)
    smallN=1e-6;
    wn=0;
    N=length(x);
    for i=1:N-1
        e1=[x(i)-P(1),y(i)-P(2)]/sqrt((x(i)-P(1))^2+(y(i)-P(2))^2);
        e2=[x(i+1)-P(1),y(i+1)-P(2)]/sqrt((x(i+1)-P(1))^2+(y(i+1)-P(2))^2);
        theta1=atan2(e1(2),e1(1));
        if(theta1<0)
            theta1=theta1+2*pi;
        end
        theta2=atan2(e2(2),e2(1));;
        if(theta2<0)
            theta2=theta2+2*pi;
        end
        %theta2=theta2-(sign(theta2)-1)*pi;
        if(theta2<theta1-pi)
            theta=theta2-(theta1-2*pi);           
        elseif(theta2>theta1+pi)
            theta=theta2-(theta1+2*pi);           
        else
            theta=theta2-theta1;
        end
        wn=wn+theta;
    end
    i=N;
    e1=[x(i)-P(1),y(i)-P(2)]/sqrt((x(i)-P(1))^2+(y(i)-P(2))^2);
    e2=[x(1)-P(1),y(1)-P(2)]/sqrt((x(1)-P(1))^2+(y(1)-P(2))^2);
        theta1=atan2(e1(2),e1(1));
        if(theta1<0)
            theta1=theta1+2*pi;
        end
        theta2=atan2(e2(2),e2(1));;
        if(theta2<0)
            theta2=theta2+2*pi;
        end
        if(theta2<theta1-pi)
            theta=theta2-(theta1-2*pi);           
        elseif(theta2>theta1+pi)
            theta=theta2-(theta1+2*pi);           
        else
            theta=theta2-theta1;
        end
        wn=wn+theta;
    wn=wn/(2*pi);
    if(abs(wn)<smallN)
        wn=0;
    end
    if(abs(wn-1)<smallN)
        wn=1;
    end
    
end



function wn=windingNumber1(P,x,y)
    wn=0;
    N=length(x);
    for i=1:N-1
        e1=[x(i)-P(1),y(i)-P(2)]/sqrt((x(i)-P(1))^2+(y(i)-P(2))^2);
        e2=[x(i+1)-P(1),y(i+1)-P(2)]/sqrt((x(i+1)-P(1))^2+(y(i+1)-P(2))^2);
        theta1=atan2(e1(2),e1(1));
        theta1=theta1-(sign(theta1)-1)*pi;
        theta2=atan2(e2(2),e2(1));;
        theta2=theta2-(sign(theta2)-1)*pi;
        theta=theta2-theta1;
        wn=wn+theta;
    end
    i=N;
    e1=[x(i)-P(1),y(i)-P(2)]/sqrt((x(i)-P(1))^2+(y(i)-P(2))^2);
    e2=[x(1)-P(1),y(1)-P(2)]/sqrt((x(1)-P(1))^2+(y(1)-P(2))^2);
    theta1=acos(e1(1));
    theta2=acos(e2(1));
    theta=theta2-theta1;
    wn=wn+theta;
    wn=wn/(2*pi);
end


function wn=windingNumber0(P,x,y)

    wn=0;
    N=length(x);
    for i=1:N-1
        e1=[x(i)-P(1),y(i)-P(2)]/sqrt((x(i)-P(1))^2+(y(i)-P(2))^2);
        e2=[x(i+1)-P(1),y(i+1)-P(2)]/sqrt((x(i+1)-P(1))^2+(y(i+1)-P(2))^2);
        wn=wn+acos(e1(1)*e2(1)+e1(2)*e2(2));
    end
    i=N;
    e1=[x(i)-P(1),y(i)-P(2)]/sqrt((x(i)-P(1))^2+(y(i)-P(2))^2);
    e2=[x(1)-P(1),y(1)-P(2)]/sqrt((x(1)-P(1))^2+(y(1)-P(2))^2);
    wn=wn+acos(e1(1)*e2(1)+e1(2)*e2(2));        
end


function e3=computebissectriz(x,y)
    e1=[x(1)-x(2),y(1)-y(2)]/sqrt((x(1)-x(2))^2+(y(1)-y(2))^2);
    e2=[x(3)-x(2),y(3)-y(2)]/sqrt((x(3)-x(2))^2+(y(3)-y(2))^2);
    e3=e1+e2;
    if(sqrt(e3(1)^2+e3(2)^2)==0) %the 2 vectores are colinear the bisectriz is perpendicular
        e3=[e1(2),-e1(1)];
    else
        e3=e3/sqrt(e3(1)^2+e3(2)^2);     
    end
end
