% communication systems, Armin Panjehpour - 98101288
% questions are all in this file and they`re seperated in sections
clear; clc; close all;

%% question 1 - multi-path channels
clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% part 1.1
fs = 1000; % sample rate
t = 0:1/fs:10; % time vector
f0 = 2; % in Hz, x(t) frequency

x1 = cos(2*pi*f0.*t);

figure;
plot(t,x1);
xlabel('Time(s)');
ylabel('x1');
title('x1(t)');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% part 1.2
% multi-path channel 

% ai are the impulses coeffiecients which follow a guassian distribution

mu = 0;
sigma = 0.5;
Tm = 10;
n = [1 10]; % length of a
rng(0);
a = normrnd(mu,sigma,n); % coeffiecients of impulses
tau = unifrnd(0,Tm,n); % shiftings of impulses

% syms t;
% hc = ; % time domain of the channel
syms f;
Hc = sum(a.*exp(-1i*2*pi*f.*tau)); % frequency domain of the channel
fplot(f,abs(Hc)); % magnitude(frequency domain of the channel)









%%