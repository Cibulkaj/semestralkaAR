H10 = 0.5;
H20 = 0.3;
syms H1 H2 Sp S2 Q S cp c2 g K1 K2 s;
K1 = cp*Sp/S*sqrt(2*g);
K2 = c2*S2/S*sqrt(2*g);
f = [-K1*sqrt(H1-H2)+1/S*Q;K1*sqrt(H1-H2)-K2*sqrt(H2)];
x=[H1,H2];
u=Q;


A=jacobian(f,x)
A2 = [diff(f,x(1));diff(f,x(2))];
B=diff(f,u)
C=[0 1];
P=C*(inv(s*eye(2)-A))*B
g=9.81;
S=25*10^(-4);
cp=0.6;
c2=0.6;
%%%%%%%%
P0 = simplify(P)
P = eval(P0);
P0 = collect(P0,s);
H1 = H10;
H2 = H20;
Q = 1.5*10^(-4);
fnum = eval(f);


[S2,SP] = solve(fnum);
% % Sp=eval(solve(fnum(2,1)));
% % Sp2 = Sp
% % S2 = sovlve(subs(fnum(2,1),'Sp',Sp2));
% % 
% % S2 = eval(S2);
% % Sp = eval(Sp);
% % P0= eval(P0);
% % P0 = collect(P0,S);
w=10;
subs(P0,s,sqrt(-1)*w);
