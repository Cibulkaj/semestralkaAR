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
%ale mam to na hovno
close all;
K = ureal('K', 10, 'Percentage', 15);
tau = ureal('tau', 1, 'Percentage', 15);
lambda = ureal('lambda', .15, 'Percentage', 15);

P = tf(K, [tau 1], 'InputDelay', lambda.Nominal);
N = 10;
P10 = usample(P,N);
P0 = tf(P.Nominal);
omg = logspace(-2,3,100);
W2 = tf([1.15*tau.Nominal, 1.15], [.85*tau.Nominal, 1], 'InputDelay', .15*lambda.Nominal)-1;

FRP0 = squeeze(freqresp(P0, omg));
reP0 = real(FRP0);
imP0 = imag(FRP0);

FRP10 = squeeze(freqresp(P10, omg));
reP10 = real(FRP10);
imP10 = imag(FRP10);

figure; plot(reP0, imP0, 'b');
hold on; plot(reP10, imP10, 'r'); 
%nyquist();

indexy = [20 99];
indReP10 = reP10(indexy,:);
indImP10 = imP10(indexy,:);
plot(indReP10, indImP10, 'kx');

FRW2 = freqresp(W2, omg);
AW2 = abs(squeeze(FRW2));
OAW2 = [AW2(1:73); AW2(73)*ones(length(AW2)-73,1)];
%circle(FRP0(indexy(1)),AW(indexy(1)));

figure; 
for i=1:N
    plot(abs((FRP10(:,i)-FRP0))/abs(FRP0), 'k');
    hold on;
end

plot(OAW2, 'b');