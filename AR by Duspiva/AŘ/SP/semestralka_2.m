clear all; close all; clc;
%%
%Inicializace

Q10 = 1.5e-4;
T = 0.5;
Qu = tf(Q10, [T 1]);
load('tfs.mat');
W2 = Wa; clear Wb;
W2 = minreal(W2/P0a); %aditivn� -> multiplikativn� v�hov� fce
W1 = 1/2;
% W1 = tf(1/2,[1/(2*pi*0.30),1]);
P0 = minreal(Qu*P0a);

% %n�vrh regul�toru pomoc� mixsyn.. :)
% %p�eklad do �e�i Janina kmene
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
% %bohu�el pro takov�hle W1 a W3 dost na pika�u.. :)

%%
%Regul�tor

% C0 = ltiblock.pid('C','pi');  % tunable PI
% % L0 = C0*P0;
% P = augw(P0,W1,[],W2);
% L0=C0*P;
% options = hinfstructOptions('TargetGain', 2);
% T = hinfstruct(L0, options);
% showTunable(T);

K = 4;
Ki = 0.3;
C = tf([K Ki], [1 0]);

%%
%P�enosov� funkce
L0 = C*P0;
S0 = minreal(1/(1+L0));
T0 = minreal(L0/(1+L0));

%%
%Testov�n� robustn� kvalita - 1a), 1b)
figure; hold on;
sup = norm(W1*S0 + W2*T0, inf)
bodemag(W1*S0 + W2*T0);
plot(xlim, [0 0], 'r--');

%Grafick� test robustn� stability
O = nyquistoptions;
O.ShowFullContour = 'off'; 
figure; hold on; set(gca, 'box', 'on');
for k = 1:length(Ps)
    L = C*Qu*Ps{k};
    nyquist(L, O, 'b');
end
nyquist(L0, O, 'r');
title('Nyquistovo kriterium stability')


%%
%�tlum na na 10 rad/s
Tmax = minreal(C*Qu*(P0a+Wa)/(1+C*Qu*(P0a+Wa)))
Tmin = minreal(C*Qu*(P0a-Wa)/(1+C*Qu*(P0a-Wa)))
gain10 = abs(squeeze(freqresp(Tmax,10)))
gain10 = 20*log10(gain10)

figure; hold on;
% bodemag(Tmax);
bodemag(T0);
%%
%Zes�len� energie �umu m��en�
gainEnergy = norm(Tmin, inf);

%%
%Zes�len� �umu m��en� a chyby na v�stupu
w = 2*pi*50;
gain50 = abs(squeeze(freqresp(Tmin, w)))

Smax = minreal(1/(1+C*Qu*(P0a+Wa)));
w = 2*pi*0.1;
gain01 = abs(squeeze(freqresp(Smax, w)))
% bodemag(Smax)
bodemag(S0);
%%
%Zes�len� na v�stupu p�i poru�e s omezenou energi� na vstupu ��zen�ho syst�mu 
Diy = minreal((Qu*(P0a-Wa))/(1+C*Qu*(P0a-Wa)));
gainMax = norm(Diy, 2)