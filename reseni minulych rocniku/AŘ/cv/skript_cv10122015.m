%%
%Model dmychadla
close all; clear all; clc;
J = .0475;  %kgm^2
m = 1.5;    %kg
r = .25;    %m
g = 9.81;   %m/s^2
d = .1;     
l = .05;    %m

%syms s
% C = 20*(s+25)/(s+300);
% W1 = 10/(s^2 + 2*s + 10);
% S = 1 / (1+P.Nominal*C);
% 
% P = (r + delta*r) / (J*s^2+d*s+m*g*l);
r = ureal('r', .25, 'Percentage', 20);
P = tf(r, [J d m*g*l]);
C = tf([20 20*25], [1 300]);    %není robustnì kvalitní   
%C = tf([20 50*25], [1 300]);   %je robustnì kvalitní :)
P0 = tf(P.Nominal);

T = P0*C/(1+P0*C);
S = 1/(1+P0*C);
W1 = tf(10, [1 2 10]);
W2 = .2;

omg = logspace(-2,3,1000);       
W1S = W1*S;
FW1S = squeeze(freqresp(W1S, omg));
AFW1S = abs(FW1S);
W2T = W2*T;
FW2T = squeeze(freqresp(W2T, omg));
AFW2T = abs(FW2T);
plot(AFW1S)     %mìlo by ||FW1S||_inf <= 1
figure;
plot(AFW2T)     %taky <= 1
figure
plot(AFW1S+AFW2T)
hold on
plot(ones(1,length(AFW1S)), 'r')

%%
%cviko 10.12.2015
P2 = augw(P,W1,[],W2);
[K,CL,GAM] = hinfsyn(P2);

pole(CL);           %poly uzavøené smyèky
pole(lft(P2,K));    %totéž

norm(lft(P2,K), inf);	%Hinf norma = GAM
K = tf(K);
K = minreal(K);

S2 = 1/(1+P0*K);
T2 = P0*K/(1+P0*K);
W2T2 = W2*T2;
L2 = P0*K;
W1S2 = W1*S2;
FRWS2 = freqresp(W1S2, omg);
AW1S2 = (abs(squeeze(FRWS2)));
FRW2T2 = freqresp(W2T2, omg);
AW2T2 = (abs(squeeze(FRW2T2)));

figure; plot(log10(omg), AW1S2+AW2T2);
hold on; plot(log10(omg), ones(1,length(AW1S2)), 'r');

L3 = P*K;
nyquist(usample(L3,10), 'b', L3.Nominal, 'r');

close all;  
%%
%model reaktoru
b1 = .00035; b2 = .00134; b3 = .00081; b4 = .00341; b5 = .00105; b6 = .00036;
tau = 5.974e-5;
b = b1+b2+b3+b4+b5+b6;
l1 = .01249; l2 = .03175; l3 = .10945; l4 = .31738; l5 = 1.3521; l6 = 8.66835;

A = [-b/tau l1 l2 l3 l4 l5 l6;
    b1/tau -l1 0 0 0 0 0;
    b2/tau 0 -l2 0 0 0 0;
    b3/tau 0 0 -l3 0 0 0;
    b4/tau 0 0 0 -l4 0 0;
    b5/tau 0 0 0 0 -l5 0;
    b6/tau 0 0 0 0 0 -l6];
B = [(1-b) b1 b2 b3 b4 b5 b6]'/tau;
C = zeros(1,7); C(1) = 1;

vlc = eig(A);

sys = ss(A,B,C,[0]);
G1 = tf(sys);
G1 = minreal(G1);
s = tf([1 0], 1);
G1mod = G1*(1/(s+1)*exp(-.2*s));

w = logspace(-1,3,10000);
[re,im] = nyquist(G1mod,w);
a = squeeze(re)';
b = squeeze(im)';

u = -.5;
v = -.866;      %bezpeènost ve fázi 60°
P1 = [u,v];
u1 = -.5;
v1 = 0;
P2 = [u1, v1];

k = (a.*u+b.*v)./(a.*a+b.*b);
ki = -((a.*v-b.*u).*w)./(a.*a+b.*b);

k1 = (a.*u1+b.*v1)./(a.*a+b.*b);
ki1 = -((a.*v1-b.*u1).*w)./(a.*a+b.*b);

figure;
hold on;

plot(ki,k,'b',ki1,k1,'g');
axis([0,0.06,0,0.035]);

[ki, k] = ginput(1);
figure;
nyquist(G1mod*(k+ki/s),w);
hold on
circle(0+0*sqrt(-1), 1); 