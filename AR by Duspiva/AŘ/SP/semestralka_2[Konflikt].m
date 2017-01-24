clear all; close all; clc;
%%
%Inicializace

Q10 = 1.5e-4;
T = 0.5;
Qu = tf(Q10, [T 1]);
load('tfs.mat');
W2 = Wa; clear Wb;
W2 = minreal(W2/P0a); %aditivní -> multiplikativní váhová fce
W1 = 1/2;
% W1 = tf(1/2,[1/(2*pi*0.30),1]);
P0 = minreal(Qu*P0a);

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
% %bohužel pro takovýhle W1 a W3 dost na pikaèu.. :)

%%
%Regulátor

% C0 = ltiblock.pid('C','pi');  % tunable PI
% % L0 = C0*P0;
% P = augw(P0,W1,[],W2);
% L0=C0*P;
% options = hinfstructOptions('TargetGain', 2);
% T = hinfstruct(L0, options);
% showTunable(T);
% K = 0:0.05:20; Ki = 0:0.05:20;
% par = [];
% for k = 1:length(K)
% for l = 1:length(Ki);
K = 4;
Ki = 0.3;
C = tf([K Ki], [1 0]);

% Ti = 1/(2*pi*0.01);
% alfa = 10;
% Fi = tf([Ti 1], [alfa*Ti 1]);
% Td = 1/(2*pi*0.01);
% alfa = 0.8;
% Fd = tf(0.1*[Td 1], [Td/alfa 1]);


% C = tf([K(k) Ki(l)], [1 0]);
%%
%Pøenosové funkce
L0 = C*P0;
S0 = minreal(1/(1+L0));
T0 = minreal(L0/(1+L0));

%%
%Testování robustní kvalita - 1a), 1b)
figure; hold on;
sup = norm(W1*S0 + W2*T0, inf);
bodemag(W1*S0 + W2*T0);
plot(xlim, [0 0], 'r--');

%Grafický test robustní stability
figure; hold on; set(gca, 'box', 'on');
for k = 1:length(Ps)
    L = C*Qu*Ps{k};
    nyquist(L, 'g');
end
nyquist(L0, 'b');

%%
%Útlum na na 10 rad/s
Tmax = minreal(C*Qu*(P0a+Wa)/(1+C*Qu*(P0a+Wa)));
Tmin = minreal(C*Qu*(P0a-Wa)/(1+C*Qu*(P0a-Wa)));
gain10 = abs(squeeze(freqresp(Tmax,10)));
gain10 = 20*log10(gain10);

figure; hold on;
% bodemag(Tmax);
bodemag(T0);
%%
%Zesílení energie šumu mìøení
gainEnergy = norm(Tmin, inf);

%%
%Zesílení šumu mìøení a chyby na výstupu
w = 2*pi*50;
gain50 = abs(squeeze(freqresp(Tmin, w)));

Smax = minreal(1/(1+C*Qu*(P0a+Wa)));
Smin = minreal(1/(1+C*Qu*(P0a-Wa)));
w = 2*pi*0.1;
gain01 = abs(squeeze(freqresp(Smax, w)));
% bodemag(Smax)
bodemag(S0);
%%
%Zesílení na výstupu pøi poruše s omezenou energií na vstupu øízeného systému 
Diy = minreal((Qu*(P0a-Wa))/(1+C*Qu*(P0a-Wa)));
gainMax = norm(Diy, 2);
% if sup <= 1 && gainEnergy <= 1.5 && gainMax ~= inf && gain01 < 1.02
%     par(:,length(par)+1) = [K(k);Ki(l)];
% end
% end
% end
%%
%Zesílení ve smyslu maximální hodnoty
% maxValueD = norm(Smax,1);
% maxValueN = norm(Tmax,1);

[As,Bs,Cs,Ds]=tf2ss(Tmin.num{1}, Tmin.den{1});
Tss = ss(As,Bs,Cs,Ds);
[L1normT,err,U,L,tol,niter]=l1norm(Tss,1e-12,24);

[At,Bt,Ct,Dt]=tf2ss(Smin.num{1}, Smin.den{1});
Sss = ss(At,Bt,Ct,Dt);
[L1normS,err,U,L,tol,niter]=l1norm(Sss,1e-12,24);

%%
%Návrh PI regulátoru pomocí metody robustních regionù
% w = logspace(-1,3,10000);
% [re,im] = nyquist(P0,w);
% a = squeeze(re)';
% b = squeeze(im)';
% 
% u = -.5;
% v = -.866; %bezpeènost ve fázi 60°
% % P1 = [u,v];
% u1 = -.5;
% v1 = 0; %bezpeènost v zesílení 2
% % P2 = [u1, v1];
% 
% k = (a.*u+b.*v)./(a.*a+b.*b);
% ki = -((a.*v-b.*u).*w)./(a.*a+b.*b); %výpoèet robustních regionù pro bezp. ve fázi
% 
% k1 = (a.*u1+b.*v1)./(a.*a+b.*b);
% ki1 = -((a.*v1-b.*u1).*w)./(a.*a+b.*b); %výpoèet robustních regionù pro bezp. v zesílení
% 
% figure;
% hold on;
% 
% plot(ki,k,'b');
% % axis([-100,100,-5,5]);
% 
% [ki, k] = ginput(1);
% C = tf([k ki], [1 0]);
% L0 = minreal(C*P0);
% figure;
% nyquist(L0);