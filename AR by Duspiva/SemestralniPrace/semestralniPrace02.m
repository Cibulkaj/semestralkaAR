%semestralniPrace02

%linearizce2
syms H1 H2 k1 k2 Sp S2 Q S cp c2 g;
g = 9.81;
S = 25*10e-4;
cp = 0.6;
c2 = 0.6;
H10 = 0.5;
H20 = 0.3;


%nelinearni rovnice
f = [-((1/S)*cp*Sp*sqrt(2*g)) + 1/S*Q;
    ((1/S)*cp*Sp*sqrt(2*g))*sqrt(H1-H2)-((1/S)*c2*S2*sqrt(2*g))*sqrt(H2)];
x = [H1;H2];
u = Q;

%linearizace
A = jacobian(f,x);
B = diff(f,u);
C = [1 0;0 1];

H1 = H10;
H2 = H20;
Q = Qv;
fnum = eval(f);
[S2, Sp] = solve(fnum);
Sp = eval(Sp);
S2 = eval(S2);

A = eval(A);
B = eval(B);
P = ss(A,B,C,0);


