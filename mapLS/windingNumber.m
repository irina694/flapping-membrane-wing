function wn=checkwindingNumber(x,y)

    epsilon=1e-6;
    N=length(x);
    x1=x(1:N-1);
    x2=x(2:N);
    y1=y(1:N-1);
    y2=y(2:N);
    
    L=min(abs([x1-x2,y1-y2]));
    dL=L*epsilon;
    
    for i=2:N-1
        %computes bissectriz unity vector 
        eB=computebissectriz([x(i-1),x(i),x(i+1)],[y(i-1),y(i),y(i+1)]);
        P1=[x(2),y(2)]+dL*eB;
        P2=[x(2),y(2)]-dL*eB; 
        wn1=windingNumber(P1,x,y);
        
    end
    
    

end



function wn=windingNumber(P,x,y)

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
    e3=e3/sqrt(e3(1)^2+e3(2)^2);     
end
