clear all;
Q10 = 1.5e-4;
T = 0.5;
Qu = tf(Q10, [T 1]);

P01 = tf([18.1818],[1 0.2279 0.0062]);
P02 = tf([20],[1 0.25 0.0075]);
P03 = tf([16.67],[1 0.2083 0.005208]);

Wa=tf([0 0 1.81818181818182 8.88178419700125e-16 -0.0123966942148760],[1 0.477272727272727 0.0705165289256198 0.00325413223140496 4.64876033057851e-05]);

P1 = P01*Qu;
P2 = P02*Qu;
P3 = P03*Qu;

Kp = 10;
Ki = 0.10;
C = tf([Kp Ki], [1 0]);

%GMK

L1 = P1*C;
L2 = P2*C;
L3 = P3*C;

S1 = minreal(1/(1+L1))
S2 = minreal(1/(1+L2));
S3 = minreal(1/(1+L3));

T1 = minreal(L1/(1+L1))
T2 = minreal(L2/(1+L2));
T3 = minreal(L3/(1+L3));

Si1 = minreal(C/(1+L1));
Si2 = minreal(C/(1+L2));
Si3 = minreal(C/(1+L3));

So1 = minreal(P1/(1+L1));
So2 = minreal(P2/(1+L2));
So3 = minreal(P3/(1+L3));

sys1 = 1+L1;
sys2 = 1+L2;
sys3 = 1+L3;

closed = feedback(L1,1)




data = impulse(-closed);
step(closed)
norma = norm(data,1);

y = db2mag(norma)




stabilni = true;

figure
suptitle('Nominalni prenos P1');
hold on;
subplot(2,2,1);
pzmap(S1);
[p1,z] = pzmap(S1);
title('gmk S');
subplot(2,2,2);
pzmap(T1);
[p2,z] =pzmap(T1);
title('gmk T');
subplot(2,2,3);
pzmap(Si1);
[p3,z] =pzmap(Si1);
title('gmk Si');
subplot(2,2,4);
pzmap(So1);
[p4,z] =pzmap(So1);
title('gmk So');
for i=1:length(p1)
    if(real(p1(i))>0)
        disp('poly');
        stabilni = false;
    end
    if(real(p2(i))>0)
        disp('poly');
        stabilni = false;
    end
    if(real(p3(i))>0)
        disp('poly');
        stabilni = false;
    end
    if(real(p4(i))>0)
        disp('poly');
        stabilni = false;
    end
end

figure
suptitle('Max Prenos P2');
hold on;
subplot(2,2,1);
pzmap(S2);
[p1,z] = pzmap(S2);
title('gmk S');
subplot(2,2,2);
pzmap(T2);
[p2,z] =pzmap(T2);
title('gmk T');
subplot(2,2,3);
pzmap(Si2);
[p3,z] =pzmap(Si2);
title('gmk Si');
subplot(2,2,4);
pzmap(So2);
[p4,z] =pzmap(So2);
title('gmk So');
for i=1:length(p1)
    if(real(p1(i))>0)
        disp('poly');
        stabilni = false;
    end
    if(real(p2(i))>0)
        disp('poly');
        stabilni = false;
    end
    if(real(p3(i))>0)
        disp('poly');
        stabilni = false;
    end
    if(real(p4(i))>0)
        disp('poly');
        stabilni = false;
    end
end

figure
suptitle('Min Prenos P3');
hold on;
subplot(2,2,1);
pzmap(S3);
[p1,z] = pzmap(S3);
title('gmk S');
subplot(2,2,2);
pzmap(T3);
[p2,z] =pzmap(T3);
title('gmk T');
subplot(2,2,3);
pzmap(Si3);
[p3,z] =pzmap(Si3);
title('gmk Si');
subplot(2,2,4);
pzmap(So3);
[p4,z] =pzmap(So3);
title('gmk So');
for i=1:length(p1)
    if(real(p1(i))>0)
        disp('poly');
        stabilni = false;
    end
    if(real(p2(i))>0)
        disp('poly');
        stabilni = false;
    end
    if(real(p3(i))>0)
        disp('poly');
        stabilni = false;
    end
    if(real(p4(i))>0)
        disp('poly');
        stabilni = false;
    end
end
stabilni;



%Nyquist
figure
hold on
title('Nyquistùv graf');

%dodelat primku!!
nyquist(L1,L2,L3);
legend('L0','L1','L2');


%robustnost
W1 = 1/2;

W21 = Wa/P01;
W22 = Wa/P02;
W23 = Wa/P03;

W1S1 = W1*S1;
W1S2 = W1*S2; 
W1S3 = W1*S3; 

W2T1 = W21*T1;
W2T2 = W22*T2;
W2T3 = W23*T3;

Oms = logspace(-2,4,1000);

FRW1S1 = freqresp(W1S1,Oms);
FRW1S2 = freqresp(W1S2,Oms);
FRW1S3 = freqresp(W1S3,Oms);


FRW2T1 = freqresp(W2T1,Oms);
FRW2T2 = freqresp(W2T2,Oms);
FRW2T3 = freqresp(W2T3,Oms);

ABSS1 = abs(squeeze(FRW1S1));
ABSS2 = abs(squeeze(FRW1S2));
ABSS3 = abs(squeeze(FRW1S3));

ABST1 = abs(squeeze(FRW2T1));
ABST2 = abs(squeeze(FRW2T2));
ABST3 = abs(squeeze(FRW2T3));

maxWS = max(ABSS1);
maxWT = max(ABST1);

soucet1 = ABSS1 + ABST1;
soucet2 = ABSS2 + ABST2;
soucet3 = ABSS3 + ABST3;

figure;
hold on;
grid on;

title('Robustni kvalita rizeni |W1S|+|W2T|');
plot(log10(Oms),soucet1);
plot(xlim, [1 1], 'r--');

legend('pro P0','hranice');
maximum = zeros(3,1);
maximum(1) = max(soucet1);
maximum(2) = max(soucet2);
maximum(3) = max(soucet3);
maximumAll = max(maximum);

if(maximumAll > 1)
    disp('robustnost');
    stabilni = false;
end




%robustnost - pozadavek na citlivostni funkci
figure
hold on;
grid on;
bode(T1);
bode(T2);
bode(T3);
limit = 20*log10(2);
[mag1,phase] = bode(T1);
[mag2,phase] = bode(T2);
[mag3,phase] = bode(T3);
maximum = zeros(3,1);
maximum(1) = max(mag1);
maximum(2) = max(mag2);
maximum(3) = max(mag3);
legend('T0','T1','T2');
maximumAll = max(maximum);

%Omezená šíøka pásma
[mag1,phase] = bode(T1,10);
[mag2,phase] = bode(T2,10);
[mag3,phase] = bode(T3,10);

magdb1 = 20*log10(mag1);
magdb2 = 20*log10(mag2);
magdb3 = 20*log10(mag3);

if(magdb1 > -70)
    disp('zesileni na fr');
    stabilni = false;
end
if(magdb2 > -70)
    disp('zesileni na fr');
    stabilni = false;
end
if(magdb3 > -70)
    disp('zesileni na fr');
    stabilni = false;
end

%Omezeni energie libovolneho sumu
energie1 = norm(T1,inf);
energie2 = norm(T2,inf);
energie3 = norm(T3,inf);

if(energie1 > 1.5)
    disp('energie');
    stabilni = false;
end
if(energie2 > 1.5)
    disp('energie');
    stabilni = false;
end
if(energie3 > 1.5)
    disp('energie');
    stabilni = false;
end

stabilni


%Dvojka - harmonicka porucha
%dosazujeme do S a T
w1 = 2*3.14*50;
w2 = 2*3.14*0.1;
zesileni50Hz1 = abs(freqresp(-T1,w1));
zesileni50Hz2 =abs(freqresp(-T2,w1));
zesileni50Hz3 =abs(freqresp(-T3,w1));

zesileni01Hz1 = abs(freqresp(S1,w2));
zesileni01Hz2 =abs(freqresp(S2,w2));
zesileni01Hz3 =abs(freqresp(S3,w2));
%bode(S1,S2,S3)

%Trojka - zesileni externich signalu
%bode(So1,So2,So3);
norma1 = norm(So1,2);
norma2 = norm(So2,2);
norma3 = norm(So3,2);

%ètyøka - zesileni s omezenim
figure
bode(S1);
title('Citlivostní funkce');
legend('S0');
w=0.21;
zesileniCitl = 2.63 %[db]


figure
bode(T1);
title('Komplementární citlivostní funkce');
legend('T0');
w=0.097;
zesileniKomp = -0.758 %[db]

sim('modelPoruchaCtyrka.slx')
figure;
hold on;
grid on;
plot(simout.Time,simout.Data(:,1));
plot(simout.Time,simout.Data(:,2));

maximum1 = norm(simout.Data(:,2), Inf)
plot(xlim, [maximum1 maximum1], 'r--');

legend('Porucha mìøení','Výstup systému','Maximum výstupu systému pøi chybì');
title('Pùsobení chyby mìøení systému.');

figure;
hold on;
grid on;
plot(simout.Time,simout1.Data(:,1));
plot(simout.Time,simout1.Data(:,2));

maximum2 = norm(simout1.Data(:,2), Inf)
plot(xlim, [maximum2 maximum2], 'r--');

legend('Porucha výstupu','Výstup systému','Maximum výstupu systému pøi chybì');
title('Pùsobení chyby výstupu systému.');


maxDb1 = 20*log10(maximum1)
maxDb2 = 20*log10(maximum2)







