%nominalni kvalita rizeni
close all;
J = 0.0475;
m = 1.5;
r = ureal('r',0.25,'Mode','Range','Percentage',20);
g = 10;
d = 0.1;
l = 0.05;

W1 = tf([10],[1 2 10]);
W2 = 0.2;
P = tf([r],[J d m*l*g]); %prenos
P0 = tf(P.Nominal) % nominalni prenos
C = 20*tf([1 25],[1 300]);
T = (P0*C)/(1+P0*C); %komplementarni citlivostni funkce
pole(T) %je to stabilni?
S = 1/(1+P*C);  % citlivostni funkce
W1S = W1*S; 
L = P0*C; %prenos otevrene smycky
Oms = logspace(-2,4,1000);
FRW1S = freqresp(W1S,Oms);
AW1S = abs(squeeze(FRW1S));
FRS = freqresp(S,Oms);
AS = abs(squeeze(FRS));
FRW1 = freqresp(W1,Oms);
AW1 = abs(squeeze(FRW1));
AIW1=1/AW1;

figure
plot(log10(Oms),AW1S);
legend('|W1S|');

%robustni kvalita rizeni
%W1s mame
W2T = W2*T;
FRW1S = freqresp(W1S,Oms);
FRW2T = freqresp(W2T,Oms);
ABSS = abs(squeeze(FRW1S));
ABST = abs(squeeze(FRW2T));

soucet = ABSS + ABST;

figure;
plot(log10(Oms),soucet);
legend('|W1S|+|W2T|');
title('Robustni kvalita rizeni');
grid on

infnormABST = max(ABST);

figure;
plot(log10(Oms),ABST);
title('Robustni stabilita rizeni');
legend('|W2T|');

grid on

%grafický test, zda se kru?nice neprotínají
FRL = squeeze(freqresp(L,Oms));

figure;
plot(real(FRL),imag(FRL));
hold on
grid on
AW2L2 = abs(squeeze(freqresp(W2*L,Oms)));
fr = 408;
v_ind = [-1,FRL(fr)];
v_ind2 = [AW1(fr),AW2L2(fr)];

fi=0:1:360;
% Grade --> Radian conversion.
fi=fi*pi/180;
R = AW1(fr);
px=ones(size(fi))*real(v_ind(1))+v_ind2(1)*cos(fi);
py=ones(size(fi))*imag(v_ind(1))+v_ind2(1)*sin(fi);
plot(px,py,'-k');
hold on
R = AW2L2(fr);
px=ones(size(fi))*real(v_ind(2))+v_ind2(2)*cos(fi);
py=ones(size(fi))*imag(v_ind(2))+v_ind2(2)*sin(fi);
plot(px,py,'-k');
hold on

%smutny - prekryva se to tak musime navrhnout novy regulator :(
P2 = augw(P,W1,[],W2);
[K,CL,GAM] = hinfsyn(P2); % reg, uz smyeka, norma nekonecno

pole(CL)
GAM %norm(lft(P2,K,Inf);
K = minreal(tf(K))

C = K*tf([1 25],[1 300]); % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
T = (P0*C)/(1+P0*C); %komplementarni citlivostni funkce
pole(T) %je to stabilni?
S = 1/(1+P*C);  % citlivostni funkce
W1S = W1*S; 
L = P0*C; %prenos otevrene smycky
Oms = logspace(-2,4,1000);
FRW1S = freqresp(W1S,Oms);
AW1S = abs(squeeze(FRW1S));
FRS = freqresp(S,Oms);
AS = abs(squeeze(FRS));
FRW1 = freqresp(W1,Oms);
AW1 = abs(squeeze(FRW1));
AIW1=1/AW1;

figure
plot(log10(Oms),AW1S);
legend('|W1S|');

%robustni kvalita rizeni
%W1s mame
W2T = W2*T;
FRW1S = freqresp(W1S,Oms);
FRW2T = freqresp(W2T,Oms);
ABSS = abs(squeeze(FRW1S));
ABST = abs(squeeze(FRW2T));

soucet = ABSS + ABST;

figure;
plot(log10(Oms),soucet);
legend('|W1S|+|W2T|');
title('Robustni kvalita rizeni');
grid on

infnormABST = max(ABST);


