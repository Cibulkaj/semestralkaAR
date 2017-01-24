% P0 = minreal(Qu*P0a);
C = tf([152 9.5],[16 0]);
Pc = tf(Q10,[0.5 1]);
P0 = P0a*Pc;
% L0 = P0a*tf(Q10,[0.5 1])*tf([152 9.5],[16 0]);


L0 = C*P0
S0 = minreal(1/(1+L0));
T0 = minreal(L0/(1+L0));
Sc = minreal(C/(1+L0));
Si = minreal(P0/(1+L0));
nyquist(L0)
L1 = C*P0b*Pc;
hold on;
nyquist(L1)
hold off

W2 = minreal(Wa/P0a) %aditivní -> multiplikativní váhová fce
W1 = 1/2;



% rltool(S0);
% title('Citlivostni funkce');
% rltool(T0);
% legend('Komplementarni citlivostni funkce');
% rltool(Sc);
% legend('Citlivostni funkce rizeni');
% rltool(Si);
% legend('Vstupni citlivostni funkce');

