function hex = drawHex(xC,yC,radius)

lSide = radius;
hHex = lSide*sin(pi/3);

x(1) = xC + lSide;
x(2) = xC + lSide/2;
x(3) = xC - lSide/2;
x(4) = xC - lSide;
x(5) = xC - lSide/2;
x(6) = xC + lSide/2;
x(7) = x(1);

y(1) = yC;
y(2) = yC + hHex;
y(3) = yC + hHex;
y(4) = yC;
y(5) = yC - hHex;
y(6) = yC - hHex;
y(7) = y(1);

for i = 1:6
    side{i} = curve2([x(i) x(i+1)],[y(i) y(i+1)]);
end

hex = geomcoerce('curve',side);