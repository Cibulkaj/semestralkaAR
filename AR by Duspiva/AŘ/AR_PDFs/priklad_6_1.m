% Priklad syntezy robustniho regulatoru metodou smisene citlivostni funkce
%=====================================================

% zadani ulohy
P0=tf([1, 10, 25],[25,75,75,25]);
G = P0;
W1=tf([13.5],[1,3,3,1]);
W2=[];
W3=tf([1,1],[0.2, 20]);    % v prednaskah je tato funkce oznacovana jako W2

% minimalizace Hinf kriteria
%        || W1*S ||
%        || W3*T || Hinf
% neboli minimalizace kompromisni podminky
% || (|W1*S0|^2+|W2*T0|^2)^(1/2) || Hinf
% viz vztah (3.11) v prednaskach

[K,CL,GAM,INFO]=mixsyn(P0,W1,W2,W3);

%
  L=G*K;  % loop transfer function
  S=inv(1+L); % Sensitivity
  T=1-S;      % complementary sensitivity

% konverze ze stavovych modelu na prenosove funkce
tf(K)
[Knum,Kden]=ss2tf(K.a,K.b,K.c,K.d)
figure(1);
bode(K)

figure(2);
[Lnum,Lden]=ss2tf(L.a,L.b,L.c,L.d);
bode(L)

figure(3);
hold on;
bode(S);
bode(T);

