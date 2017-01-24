function circle_m(s,r)
a = linspace(0,2*pi, 100);
v = s+r*exp(j*a);
plot(real(v), imag(v),'r');
end