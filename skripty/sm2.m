clear all;


Q10 = 1.5e-4;
T = 0.5;
Qu = tf(Q10, [T 1]);

P0 = tf([18.1818],[1 0.2279 0.0062]);
Wa=tf([0 0 1.81818181818182 8.88178419700125e-16 -0.0123966942148760],[1 0.477272727272727 0.0705165289256198 0.00325413223140496 4.64876033057851e-05]);

P = P0*Qu
Kp = 0.010;
Ki = 0.01;

C = tf([Kp Ki], [1 0]);



%GMK

L = P*C
S = minreal(1/(1+L))
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

close all;
for i=1:length(p1)
    if(real(p1(i))>0)
        stabilni = false;
    end
    if(real(p2(i))>0)
        stabilni = false;
    end
    if(real(p3(i))>0)
        stabilni = false;
    end
    if(real(p4(i))>0)
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






