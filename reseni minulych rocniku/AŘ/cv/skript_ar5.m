clear; %clc; 
close all;
%%
%cv5
w=50;
xi=.1;
sys=tf(w^2, [1, 2*xi*w, w^2]);
[hinf, fpeak] = norm(sys,inf);

% figure;
% bode(sys);

figure;
for xi=.1:.1:1
    sys=tf(w^2, [1, 2*xi*w, w^2]);
    bode(sys);
    hold on;
end

clear; %clc; 
close all;
%%
%druha cast
maxtime = 500;
time = 0:0.01:maxtime;
w = 0.7071;
sinwave = sin(w*time);
h2sin = norm(sinwave,2)

squarewave = square(w*time);
h2square = norm(squarewave,2)

clear; clc; close all
%%
%treti cast
s = tf('s');
%K=56/1.5;
K = 1;
G = (K*(s+3)*(s+5)/(s*(s+7)*(s+8)));
C = 1/s;
sys = feedback(C*G,1);
t = 0:0.1:150;
%u = ones(1,251);
u = t;
[y,t,x] = lsim(sys,u,t);
plot(t,y,'b',t,u,'m');

clear; clc; close all
%%
%ctvrta cast
warning off;
s = tf('s');
K = 9;
P = K*1/(10*s+1);
sim('model.slx');