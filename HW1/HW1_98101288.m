% communication systems, Armin Panjehpour - 98101288
% questions are all in this file and they`re seperated in sections
clear; clc; close all;

%% question 1 - multi-path channels
clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% part 1.1
fs = 1000; % sample rate
t = 0:1/fs:10; % time vector
f0 = 2; % in Hz, x(t) frequency



x = cos(2*pi*f0.*t);

figure;
plot(t,x);
xlabel('Time(s)');
ylabel('x');
title('x(t)');
grid on;
grid minor;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% part 1.2
% multi-path channel 

% ai are the impulses coeffiecients which follow a guassian distribution

mu = 0;
sigma = 0.5;
Tm = 0.01;
n = [1 10]; % length of a
a = normrnd(mu,sigma,n); % coeffiecients of impulses
tau = unifrnd(0,Tm,n); % shiftings of impulses


% frequency domain of the channel
syms f;
syms HcF(f);
HcF(f) = sum(a.*exp(-1i*2*pi*f.*tau)); 
figure;
fplot(abs(HcF),[-1000 1000]); % magnitude(frequency domain of the channel)
grid on;
grid minor;
title('magnitude(Hc(f))');
xlabel('f(Hz)');
legend('Tm = 0.01')

figure;
fplot(angle(HcF),[-1000 1000]); % phase(frequency domain of the channel)
grid on;
grid minor;
title('phase(Hc(f))');
xlabel('f(Hz)');
legend('Tm = 0.01')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% part 1.4
tau = (1/fs).*randi([0 floor(Tm*fs)],n); % samples
HcSampled(f) = sum(a.*exp(-1i*2*pi*f.*tau)); % frequency domain of the channel

N = 5000;

t = 0:1/fs:100-1/fs; % time vector
output = zeros(1,length(t));
for i =1:(length(t))/N
    tprime = ((i-1)*N+1)/fs:1/fs:((i)*N)/fs;
    a = normrnd(mu,sigma,n); % coeffiecients of impulses
    tau = (1/fs).*randi([0 floor(Tm*fs)],n); % samples
    for j=1:max(n)
        output(int64(tprime*fs)) = output(int64(tprime*fs)) + a(j).*cos(2*pi*f0.*((tprime-1/fs)-tau(j)));
    end
end

figure;
plot(t,output)
title('output signal')
xlabel('t');
legend('N = 5000');
grid on;
grid minor;

%% question 2 - signal reconstruction in multipath channels
% m-tapped-delay line equalizer
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% part 2.3
fs = 1000; % sample rate
t = 0:1/fs:10; % time vector
f0 = 2; % in Hz, x(t) frequency
x = cos(2*pi*f0.*t).*(stepfun(t,0)-stepfun(t,2));

figure;
plot(t,x);
grid on;
grid minor;
title('input signal')
xlabel('Time(s)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% part 2.4
syms f;
syms HcF(f);
HcF(f) = exp(-1i*10*pi*f) + 0.4*exp(-1i*10.02*pi*f); 
figure;
fplot(abs(HcF),[-100 100]); % magnitude(frequency domain of the channel)
grid on;
grid minor;
title('magnitude(Hc(f))');
xlabel('f(Hz)');


figure;
fplot(angle(HcF),[-1 1]); % phase(frequency domain of the channel)
grid on;
grid minor;
title('phase(Hc(f))');
xlabel('f(Hz)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% part 2.5 
y = cos(2*pi*f0.*(t-5)).*(stepfun(t,5)-stepfun(t,7))+...
    0.4*cos(2*pi*f0.*(t-5.01)).*(stepfun(t,5.01)-stepfun(t,7.01));

figure;
subplot(1,2,1);
plot(t,x);
grid on;
grid minor;
title('input signal');
xlabel('Time(t)');

subplot(1,2,2);
plot(t,y);
grid on;
grid minor;
title('output signal');
xlabel('Time(t)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% part 2.6 and 2.7
a = [1 0.4];

k = a(1);
k2 = a(2)/k;

mu = 0;
sigma = 0.5;
Tm = 0.01;
n = [1 10]; % length of a
tau = [5 5.01];

t0 = tau(1);
t2 = tau(2) - t0;
Heq = 1./(1+k2*exp(-1i*2*pi*f*t2)); % fourier transform of Heq

% taylor series of Heq - Ci and Ti s
syms HeqTaylor(f)
HeqTaylor = 1;
m = 10;
C = zeros(1,m+1);
T = zeros(1,m);
C(1) = 1;

for i=1:m
    HeqTaylor = HeqTaylor + (-1)^i*(k2*exp(-1i*2*pi*f*t2))^i;
    C(i+1) = ((-1)^i)*((k2)^i);
    T(i) = t2*i;
end



% outputs
yprime = cos(2*pi*f0.*(t-t0)).*(stepfun(t,t0)-stepfun(t,2+t0));% output of the ideal channel

xm = y; % reconstructed signal
RMSerror = zeros(1,10); % RMS empty vector
for i=1:m
    xm = xm + C(i+1).*(cos(2*pi*f0.*(t-5-T(i))).*(stepfun(t,5+T(i))-stepfun(t,7+T(i)))+...
    0.4*cos(2*pi*f0.*(t-5.01-T(i))).*(stepfun(t,5.01+T(i))-stepfun(t,7.01+T(i))));
    RMSerror(i) = sum(abs(yprime-xm).^2); % RMS calculation
end

figure;
plot(xm) % output of the m-tapped
title('output of the m-tapped xm')
grid on;
grid minor;

% RMS
figure;
plot(RMSerror)
title('RMS')
xlabel('m')
ylabel('error')
grid on;
grid minor;

figure;
semilogy(1:10,RMSerror);
title('Logarithmic-Yaxis RMS')
xlabel('m')
ylabel('error')
grid on;
grid minor;

figure; 
subplot(2,2,1);
plot(t,x);
grid on;
grid minor;
title('x(t)','Interpreter','latex')
xlabel('Time(s)')

subplot(2,2,2);
plot(t,y);
grid on;
grid minor;
title('y(t)','Interpreter','latex')
xlabel('Time(s)')


subplot(2,2,3);
plot(t,yprime);
grid on;
grid minor;
title('$\hat{y}(t)$','Interpreter','latex')
xlabel('Time(s)')


subplot(2,2,4);
plot(t,xm);
grid on;
grid minor;
title('$\tilde{x}(t)$','Interpreter','latex')
xlabel('Time(s)')
legend('m=10')
