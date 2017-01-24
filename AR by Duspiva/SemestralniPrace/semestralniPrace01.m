%Semestralni prace z AR
Sp0 =1.2620e-04;
S20 = 1.0305e-04;
Qv = 1.5e-04;
%linearizce
syms H1 H2 k1 k2 Sp S2 Q S cp c2 g;
g = 9.81;
S = 25*10e-4;
cp = 0.6;
c2 = 0.6;

%pomocne promenne
k1= (1/S)*cp*Sp*sqrt(2*g);
k2= (1/S)*c2*S2*sqrt(2*g);

%nelinearni rovnice
f = [-((1/S)*cp*Sp*sqrt(2*g))*sqrt(H1-H2)+1/S*Q;
    ((1/S)*cp*Sp*sqrt(2*g))*sqrt(H1-H2)-((1/S)*c2*S2*sqrt(2*g))*sqrt(H2)];
x = [H1;H2];
u = Q;

%linearizace
A = jacobian(f,x);
B = diff(f,u);
C = [1 0;0 1];

Sp = Sp0;
S2 = S20;
Q = Qv;
fnum = eval(f);
[H1, H2] = solve(fnum);
H1 = eval(H1);
H2 = eval(H2);

A = eval(A);
B = eval(B);
P = ss(A,B,C,0);
Sys = tf(P(2))

% nyquist(P(2,1))

