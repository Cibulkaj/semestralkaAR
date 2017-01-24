close all; clear all; clc;
%%
%aø cv 12.11.2015

% a = ureal('a', 5);
%a.NominalValue
r = ureal('r', 1, 'Range', [.9,1.1]);
c = ureal('c', 1, 'Range', [.9,1.1]);
P = tf(1, [r*c, 1]);
N = 10; 

P0 = P.NominalValue;
P10 = usample(P,10);

W = tf([1/.21 0],[1 1]);
omg = logspace(-3, 2, 200);
FRP0 = squeeze(freqresp(P0, omg));
reP0 = real(FRP0);
imP0 = imag(FRP0);

FRP10 = squeeze(freqresp(P10, omg));
reP10 = real(FRP10);
imP10 = imag(FRP10);

figure; plot(reP0, imP0, 'b');
hold on; plot(reP10, imP10, 'r'); 
%nyquist();

indexy = [80 199];
indReP10 = reP10(indexy,:);
indImP10 = imP10(indexy,:);
plot(indReP10, indImP10, 'kx');

FRW = freqresp(W, omg);
AW = abs(squeeze(FRW));
circle(FRP0(indexy(1)),AW(indexy(1)));
%%
%zajímavìjší system
close all;
K = ureal('K', 10, 'Percentage', 15);
tau = ureal('tau', 1, 'Percentage', 15);
lambda = ureal('lambda', .5, 'Percentage', 15);

P = tf(K, [tau 1], 'InputDelay', lambda.Nominal);
N = 10;
P10 = usample(P,N);
P0 = tf(P.Nominal);
omg = logspace(-2,3,100);
W2 = tf([1.15*tau.Nominal, 1.15], [.85*tau.Nominal, 1], 'InputDelay', .15*lambda.Nominal)-1;

FRP0 = squeeze(freqresp(P0, omg));
reP0 = real(FRP0);
imP0 = imag(FRP0);

FRP10 = (squeeze(freqresp(P10, omg)));
reP10 = real(FRP10);
imP10 = imag(FRP10);

figure; plot(reP0, imP0, 'b');
hold on;
% plot(reP10, imP10, 'rx'); 
%nyquist();

indexy = [5 40];
indReP10 = reP10(indexy,:);
indImP10 = imP10(indexy,:);
plot(indReP10, indImP10, 'kx');

FRW2 = freqresp(W2, omg);
AW2 = abs(squeeze(FRW2));
OAW2 = [AW2(1:73); AW2(73)*ones(length(AW2)-73,1)];
FRW2P0 = freqresp(W2*P0, omg);
AW2P0 = abs(squeeze(FRW2P0));
circle(FRP0(indexy(1)),AW2P0(indexy(1)));
circle(FRP0(indexy(2)),AW2P0(indexy(2)));
%%
%cvièení 19.11.2015
figure; 
for i=1:N
    plot(abs((FRP10(:,i)-FRP0))/abs(FRP0), 'k');
    hold on;
end
plot(OAW2, 'b');

figure;
plot(abs(FRP10),'g'); hold on; %|P(jw)|
plot(abs(FRP0),'r'); %|Pnom(jw)|

figure;
plot(reP0, imP0, 'b'); hold on;
indexy = [10 20 30 40 50 60 70 80 90];
indReP0 = reP10(indexy,:);
indImP0 = imP10(indexy,:);
plot(indReP0, indImP0, 'gx'); hold on;
FRW2P0 = freqresp(W2*P0, omg);
AW2P0 = abs(squeeze(FRW2P0));

for i = 1:length(indexy)
    circle(FRP0(indexy(i)), AW2P0(indexy(i))); 
    hold on;
end
close all;
%robustní stabilita
C = tf([.15*tau.Nominal, .15], [1,0]);
L0 = P0*C;
L0 = minreal(L0);
S0 = inv(1+L0);
T0 = L0/(1+L0);

FRS0 = freqresp(S0,omg);
AS0 = abs(squeeze(FRS0));
FRT0 = freqresp(T0,omg);
AT0 = abs(squeeze(FRT0));

W2T0 = squeeze(freqresp(W2*T0,omg));    %W(jw)T0(jw)
AW2T = abs(W2T0);

figure;
plot(AW2T, 'r');

W1 = tf([1,2.5],[5,0]);
FRW1 = freqresp(W1, omg);   %W1(jw)
AW1 = abs(squeeze(FRW1));   %|W1(jw)|

W1S = squeeze(freqresp(W1*S0, omg));
AW1S = abs(W1S);

figure; plot(AW1S, 'r');
%robustní kvalita øízení
figure; plot(AW2T+AW1S, 'g');