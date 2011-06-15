function area = computeArea3(x,y)
        aux=0;
        nn = length(x);
        for j=1:nn-1
            aux=aux+(x(j)-x(j+1))*(y(j)+y(j+1))/2;
        end
        aux=aux+(x(nn)-x(1))*(y(nn)+y(1))/2;
        area=abs(aux);
end
