Q10 = 1.5e-4;
T = 0.5;
Qcerpadla = tf(Q10, [T 1]);


Q10 = 1.5e-4;
H10 = 0.6;
H20 = 0.4;
S = 25e-4;
Cp = 0.6;
C2 = 0.6;
g = 9.81;

Sp = Q10/(Cp*sqrt(2*g*(H10-H20)));
S2 = (Cp*Sp*sqrt(H10-H20))/(C2*sqrt(H20));


Q = Q10;
H2 = Q^2/(2*g*S2^2*C2^2);
H1 = (Q^2/(2*g*Sp^2*Cp^2))+H2;

A1 = [-((1/S)*cp*Sp*sqrt(2*g))*sqrt(H1-H2)+1/S*Q;
    ((1/S)*cp*Sp*sqrt(2*g))*sqrt(H1-H2)-((1/S)*c2*S2*sqrt(2*g))*sqrt(H2)];

P = minreal(Qcerpadla*Sys*10);

Sys2 = ss(f,B,C,D);
Ppuvodni = tf(Sys2);

W2 = Wa; clear Wb;
W2 = minreal(W2/P0a); %aditivní -> multiplikativní váhová fce
W1 = 1/2;


% %návrh regulátoru pomocí mixsyn.. :)
% %pøeklad do øeèi Janina kmene
% W3 = W2;
% W2 = [];
% G = P0;
% 
% [K,CL,GAM,INFO]=mixsyn(P0,W1,W2,W3);
% 
% %
%   L=G*K;  % loop transfer function
%   S=inv(1+L); % Sensitivity
%   T=1-S;      % complementary sensitivity
% 
% % konverze ze stavovych modelu na prenosove funkce
% tf(K)
% [Knum,Kden]=ss2tf(K.a,K.b,K.c,K.d)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ks = 0.5;
Ki = 0.1;

CPI = tf([Ks Ki], [1 0]);

Pfin = CPI*P;
S = minreal(1/(1+Pfin));
T = minreal(Pfin/(1+Pfin));

figure; hold on;
sup = norm(W1*S0 + W2*T0, inf);
bodemag(W1*S0 + W2*T0);

