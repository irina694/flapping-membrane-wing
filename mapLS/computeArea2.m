function area = computeArea2(x,y)
      
        aux=0;
        nn = length(x);
        for j=1:nn-1
            aux=aux+(x(j)*y(j+1)-x(j+1)*y(j))/2;
        end
        aux=aux+(x(nn)*y(1)-x(1)*y(nn))/2;
        area=abs(aux);
end
