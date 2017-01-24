function circle(center, radius, color)
acc = 0.01;
l = 1;
for k = 0:acc:2*pi
    circ(1,l) = radius*cos(k) + real(center);
    circ(2,l) = radius*sin(k) + imag(center);
    l = l + 1;
end

plot(circ(1,:), circ(2,:), color);