clear all; close all; clc;
%%
%Inicializace konstant
Q10 = 1.5e-4;
H10 = 0.8;
H20 = 0.2;
S = 25e-4;
Cp = 0.6;
C2 = 0.6;
g = 9.81;

%%
%Lineariozvaný model A
Sp = Q10/(Cp*sqrt(2*g*(H10-H20)));
S2 = (Cp*Sp*sqrt(H10-H20))/(C2*sqrt(H20));

A = [-(Cp*Sp*sqrt(2*g))/(2*S*sqrt(H10-H20)) (Cp*Sp*sqrt(2*g))/(2*S*sqrt(H10-H20));
    (Cp*Sp*sqrt(2*g))/(2*S*sqrt(H10-H20)) (-Cp*Sp*sqrt(2*g))/(2*S*sqrt(H10-H20))-(C2*S2*g)/(S*sqrt(2*g*H20))];

B = [1/S;0];
C = eye(2); D = zeros(2,1);

%%
%Linearizovaný model - A Q = 1.2*Q, stejné nastavení pøep. ventilù
Q = 1.2*Q10;
H2 = Q^2/(2*g*S2^2*C2^2);
H1 = (Q^2/(2*g*Sp^2*Cp^2))+H2;

A1 = [-(Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2)) (Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2));
    (Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2)) (-Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2))-(C2*S2*g)/(S*sqrt(2*g*H2))];

%%
%Linearizovaný model - B Q = 1.2*Q, stejná výška hladin
Sp_1 = Q/(Cp*sqrt(2*g*(H10-H20)));
S2_1 = (Cp*Sp_1*sqrt(H10-H20))/(C2*sqrt(H20));

A2 = [-(Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20)) (Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20));
    (Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20)) (-Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20))-(C2*S2_1*g)/(S*sqrt(2*g*H20))];
%%
%Pøenosy Q1->H2
sys = ss(A,B,C,D);
P = tf(sys);
P = P(2)

sys = ss(A1,B,C,D);
P1 = tf(sys);
P1 = P1(2)

sys = ss(A2,B,C,D);
P2 = tf(sys);
P2 = P2(2)


O = nyquistoptions;
O.ShowFullContour = 'off'; 
% %%Nquist
% figure
% hold on
% nyquist(P);
% nyquist(P2);
% legend('Pùvodní pøítok Q1','Zvýšený pøítok Q = 1.2*Q');

%%
%Neurèitosti
%Pøípad a)
Q = 1.1*Q10;
H2 = Q^2/(2*g*S2^2*C2^2);
H1 = (Q^2/(2*g*Sp^2*Cp^2))+H2;

A1 = [-(Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2)) (Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2));
    (Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2)) (-Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2))-(C2*S2*g)/(S*sqrt(2*g*H2))];

P0a = tf(ss(A1,B,C,D));
P0a = P0a(2);
P = P
P0a = P0a
Wa = P-P0a



Ps = cell(0);
figure; hold on; set(gca, 'box', 'on');
Qmod = 1:0.03:1.2;
for k = 1:length(Qmod)
    Q = Qmod(k)*Q10;
    H2 = Q^2/(2*g*S2^2*C2^2);
    H1 = (Q^2/(2*g*Sp^2*Cp^2))+H2;

    A1 = [-(Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2)) (Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2));
        (Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2)) (-Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2))-(C2*S2*g)/(S*sqrt(2*g*H2))];
    Ps{length(Ps)+1} = tf(ss(A1,B,C,D));
    Ps{length(Ps)} = Ps{length(Ps)}(2);
    nyquist(Ps{length(Ps)}, O, 'g');
end
nyquist(P0a, 'b');
title('Neurèitost pro variantu A')
omega = logspace(-2,1e-100,9);
omega = [0 omega];
FRP0a = squeeze(freqresp(P0a,omega));
FRWa = squeeze(freqresp(Wa,omega));
for k = 1:length(omega)
    circle(FRP0a(k), abs(FRWa(k)), 'm');
%     for i=1:1
%        plot(real(FRP0a(k)) +  real(FRWa(i)),imag(FRWa(i)) -imag(FRP0a(k)),'r*'); 
%     end
    
end

%Pøípad b)
Q = 1.1*Q10;
Sp_1 = Q/(Cp*sqrt(2*g*(H10-H20)));
S2_1 = (Cp*Sp_1*sqrt(H10-H20))/(C2*sqrt(H20));

A2 = [-(Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20)) (Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20));
    (Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20)) (-Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20))-(C2*S2_1*g)/(S*sqrt(2*g*H20))];

P0b = tf(ss(A2,B,C,D));
P=P
P0b = P0b(2)

Wb = P-P0b

figure; hold on; set(gca, 'box', 'on');
Qmod = 1:0.03:1.2;
for k = 1:length(Qmod)
    Q = Qmod(k)*Q10;
    Sp_1 = Q/(Cp*sqrt(2*g*(H10-H20)));
    S2_1 = (Cp*Sp_1*sqrt(H10-H20))/(C2*sqrt(H20));
    A2 = [-(Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20)) (Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20));
    (Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20)) (-Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20))-(C2*S2_1*g)/(S*sqrt(2*g*H20))];
    
    F = tf(ss(A2,B,C,D));
    F = F(2);
    nyquist(F, O, 'g');
end
nyquist(P0b, 'b');

title('Neurèitost pro variantu B')
FRP0b = squeeze(freqresp(P0b,omega));
FRWb = squeeze(freqresp(Wb,omega));
for k = 1:length(omega)
    circle(FRP0b(k), abs(FRWb(k)), 'r');
end

%Porovnání
figure; 
hold on;
grid on;
bodemag(Wa, 'b');
bodemag(Wb, 'g');
legend('varianta A', 'varianta B');
title('Porovnání obou variant v Bodeho charakteristice.');
grid on;


% %Porovnání
% Ra = zeros(1,length(omega));
% for k = 1:length(Ra)
%     Ra(k) = abs(FRWa(k));
% end
% 
% Rb = zeros(1,length(omega));
% for k = 1:length(Rb)
%     Rb(k) = abs(FRWb(k));
% end
% save('tfs.mat', 'P0a', 'Wa', 'Ps');