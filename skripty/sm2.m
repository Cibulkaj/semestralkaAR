clear all;


Q10 = 1.5e-4;
T = 0.5;
Qu = tf(Q10, [T 1]);

P0 = tf([18.1818],[1 0.2279 0.0062]);
Wa=tf([0 0 1.81818181818182 8.88178419700125e-16 -0.0123966942148760],[1 0.477272727272727 0.0705165289256198 0.00325413223140496 4.64876033057851e-05]);

P = P0*Qu;
Kp = 10;
Ki = 0.10;

C = tf([Kp Ki], [1 0]);



%GMK

L = P*C;
S = minreal(1/(1+L));
T = minreal(L/(1+L));
Si = minreal(C/(1+L));
So = minreal(P/(1+L));

stabilni = true;

figure
hold on;
subplot(2,2,1);
pzmap(S);
[p1,z] = pzmap(S);
title('gmk S');
subplot(2,2,2);
pzmap(T);
[p2,z] =pzmap(T);
title('gmk T');
subplot(2,2,3);
pzmap(Si);
[p3,z] =pzmap(Si);
title('gmk Si');
subplot(2,2,4);
pzmap(So);
[p4,z] =pzmap(So);
title('gmk So');

sys = 1+C*P;

[p,z] = pzmap(sys);

close all;
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
%dodelat primku!!
nyquist(L);
title('nquist');
[p,z] = pzmap(L);
close all;

W1 = 1/2;
W2 = Wa/P0;
L = P0*C; %prenos otevrene smycky
W1S = W1*S; 
W2T = W2*T;
Oms = logspace(-2,4,1000);


FRW1S = freqresp(W1S,Oms);
FRW2T = freqresp(W2T,Oms);
ABSS = abs(squeeze(FRW1S));
ABST = abs(squeeze(FRW2T));

soucet = ABSS + ABST;

figure;
plot(log10(Oms),soucet);
maximum = max(soucet);

if(maximum > 1)
    disp('robustnost');
    stabilni = false;
end

legend('|W1S|+|W2T|');
title('Robustni kvalita rizeni');
grid on
close all;

%


%Omezená šíøka pásma
[mag,phase] = bode(T,10);
magdb = 20*log10(mag);

if(magdb > -70)
    disp('zesileni na fr');
    stabilni = false;
end

%Omezeni energie libovolneho sumu
energie = norm(T,inf);

if(energie > 1.5)
    disp('energie');
    stabilni = false;
end

stabilni








