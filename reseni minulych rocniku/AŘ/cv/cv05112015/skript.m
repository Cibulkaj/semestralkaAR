close all; clear all;
warning off;
%%
%cvi�en� z a� s kolegou �echem

Ts = .01;
P = 1;
I = 1;
D = 0;
w = 1;
Td = .2*pi;
zpozdeni = tf(1,1, 'ioDelay', Td);
sys = tf(1, [1 1]);
sim('model.slx');