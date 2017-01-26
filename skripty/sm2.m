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

S1 = minreal(1/(1+L1));
S2 = minreal(1/(1+L2));
S3 = minreal(1/(1+L3));

T1 = minreal(L1/(1+L1));
T2 = minreal(L2/(1+L2));
T3 = minreal(L3/(1+L3));

Si1 = minreal(C/(1+L1));
Si2 = minreal(C/(1+L2));
Si3 = minreal(C/(1+L3));

So1 = minreal(P1/(1+L1));
So2 = minreal(P2/(1+L2));
So3 = minreal(P3/(1+L3));

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
suptitle('Prenos P2');
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
suptitle('Prenos P3');
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
legend('L1','L2','L3');


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

soucet1 = ABSS1 + ABST1;
soucet2 = ABSS2 + ABST2;
soucet3 = ABSS3 + ABST3;

figure;
hold on;
grid on;

title('Robustni kvalita rizeni |W1S|+|W2T|');
plot(log10(Oms),soucet1);
plot(log10(Oms),soucet2);
plot(log10(Oms),soucet3);
plot(xlim, [1 1], 'r--');

legend('pro P1','pro P2','pro P3');
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
bode(S1);
bode(S2);
bode(S3);
limit = 20*log10(2);
[mag1,phase] = bode(S1);
[mag2,phase] = bode(S2);
[mag3,phase] = bode(S3);
maximum = zeros(3,1);
maximum(1) = max(mag1);
maximum(2) = max(mag2);
maximum(3) = max(mag3);
legend('S1','S2','S3');
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

close all;








