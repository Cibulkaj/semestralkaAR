%linearizce
syms H1 H2 k1 k2 Sp S2 Q S cp c2 g;
g = 9.81;
S = 25*10e-4;
cp = 0.6;
c2 = 0.6;
H10 = 0.5;
H20 = 0.3;
Qv = 1.5e-04;

%%
%Lineariozvaný model
Sp0 = 1.2620e-04;
S20 = 1.0305e-04;

A = [-(cp*sp*sqrt(2*g))/(2*S*sqrt(H10-H20)) (cp*sp*sqrt(2*g))/(2*S*sqrt(H10-H20));
    (cp*sp*sqrt(2*g))/(2*S*sqrt(H10-H20)) (-cp*sp*sqrt(2*g))/(2*S*sqrt(H10-H20))-(c2*s2*g)/(S*sqrt(2*g*H20))];

B = [1/S;0];
C = eye(2); D = zeros(2,1);

%%
%Linearizovaný model - Q = 1.2*Q, stejné nastavení pøep. ventilù
Q = 1.2*Q10;
H2 = Q^2/(2*g*s2^2*C2^2);
H1 = (Q^2/(2*g*Sp^2*Cp^2))+H2;

A1 = [-(Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2)) (Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2));
    (Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2)) (-Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2))-(C2*S2*g)/(S*sqrt(2*g*H2))];

%%
%Linearizovaný model - Q = 1.2*Q, stejná výška hladin
Sp_1 = Q/(Cp*sqrt(2*g*(H10-H20)));
S2_1 = (Cp*Sp_1*sqrt(H10-H20))/(C2*sqrt(H20));

A2 = [-(Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20)) (Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20));
    (Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20)) (-Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20))-(C2*S2_1*g)/(S*sqrt(2*g*H20))];
A1 = [-(Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2)) (Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2));
    (Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2)) (-Cp*Sp*sqrt(2*g))/(2*S*sqrt(H1-H2))-(C2*S2*g)/(S*sqrt(2*g*H2))];

%%
%Pøenosy Q1->H2
sys = ss(A,B,C,D);
P = tf(sys);
P = P(2);

sys = ss(A1,B,C,D);
P1 = tf(sys);
P1 = P1(2);

sys = ss(A2,B,C,D);
P2 = tf(sys);
P2 = P2(2);

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
Wa = P-P0a;

O = nyquistoptions;
O.ShowFullContour = 'off'; 

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

omega = logspace(-2,1e-100,9);
omega = [0 omega];
FRP0a = squeeze(freqresp(P0a,omega));
FRWa = squeeze(freqresp(Wa,omega));
for k = 1:length(omega)
    circle(FRP0a(k), abs(FRWa(k)), 'r');
end

%Pøípad b)
Q = 1.1*Q10;
Sp_1 = Q/(Cp*sqrt(2*g*(H10-H20)));
S2_1 = (Cp*Sp_1*sqrt(H10-H20))/(C2*sqrt(H20));

A2 = [-(Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20)) (Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20));
    (Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20)) (-Cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H10-H20))-(C2*S2_1*g)/(S*sqrt(2*g*H20))];

P0b = tf(ss(A2,B,C,D));
P0b = P0b(2);

Wb = P-P0b;

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

FRP0b = squeeze(freqresp(P0b,omega));
FRWb = squeeze(freqresp(Wb,omega));
for k = 1:length(omega)
    circle(FRP0b(k), abs(FRWb(k)), 'r');
end

%Porovnání
figure; hold on;
bodemag(Wa, 'b');
bodemag(Wb, 'r');
legend('Pripad a', 'Pripad b');

%Porovnání
Ra = zeros(1,length(omega));
for k = 1:length(Ra)
    Ra(k) = abs(FRWa(k));
end

Rb = zeros(1,length(omega));
for k = 1:length(Rb)
    Rb(k) = abs(FRWb(k));
end
save('tfs.mat', 'P0a', 'Wa', 'Ps');