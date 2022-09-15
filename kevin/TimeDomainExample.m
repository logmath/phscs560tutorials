clear; close all; clc;
%% make signal 
fs=44000; % sample rate
dt=1/fs; % sample period
T=5;  % total sample length
t=0:dt:(T-dt); % time array
N = length(t); % number of samples
[B,A]=butter(2,[500,5000]/(fs/2)); % make a band-pass filter 
noise = filter(B,A,randn(1,length(t))); % make band-limited gaussian white noise
x=sin(2*pi*4000*t) + 3*noise; % add the noise to a sine wave at 4khz

%% play sound
soundsc(x,fs) % play the sound
pause(5)

%% plot
figure(1)
plot(t,x) % plot the signal
xlabel('time (s)')
ylabel('amplitude (arbitrary)')

%% look at histogram
figure(2)
subplot 131
nb = 100; % number of bins for my histogram
h = histogram(x,nb);% make a histogram of the data 
xlabel('amplitude (arbitrary)')
ylabel('number of counts')

subplot 132
h = histogram(x,nb,'normalization','pdf');% make a pdf estimate of data
xlabel('amplitude (arbitrary)')
ylabel('estimate of probability density')

subplot 133
clipped = x;
clipped(x>3) =3;
clipped(x<-3) =-3;
h = histogram(clipped,nb,'normalization','pdf');% make a pdf estimate of data
xlabel('amplitude (arbitrary)')
ylabel('estimate of probability density')
title('what if the signal was clipped?')


%% autocorrelation
[Rxx,lags] = xcorr(x);
lags = lags/fs;
figure(3);
plot(lags,Rxx)
xlabel('\tau')
ylabel('R_{xx}')

%% delay and cross correlation
y = circshift(x,fs);
[Ryx,lags] = xcorr(x,y);
lags = lags/fs;
figure(4)
plot(lags,Ryx)
xlabel('\tau')
ylabel('R_{xy}')


%% impulse response
[hello,hellofs] = audioread('C:\Users\kevin\Documents\Sound recordings\hello hello hello.m4a');
soundsc(hello,hellofs)
pause(5)
x = hello(:,1); % just grab one channel from the stereo recording
N = length(x);
t = (0:N-1)*1/hellofs;% make t array
h = zeros(size(x)); % initialize impulse response
h(1) =1; % make it a delta function at delay 0, meaning the system does nothing
delay = 0.5;
h(hellofs*delay) = 1; % add an echo after a delay 
y = conv(x,h);
soundsc(y,hellofs)

figure(5)
subplot(3,1,1)
plot(t,x)
ylabel('x')

subplot(3,1,2)
stem(t,h)
ylabel('h')

y= y(1:length(t));
subplot(3,1,3)
plot(t,y)
xlabel('t (seconds)')
ylabel('y')

%% check correlations
[Rxx,lagsxx] = xcorr(x,x,floor(N/2)); % the third input is just to truncate the correlation at a max lag of half the duration
[Ryx,lagsxy] = xcorr(y,x,floor(N/2));
[Ryy,lagsyy] = xcorr(y,y,floor(N/2));

Ryxtest = conv(Rxx,h);
Ryxtest = Ryxtest(1:N);
Ryytest = conv(Ryxtest,h);
Ryytest = Ryytest(1:N);

figure(6)
subplot 311
plot(lagsxx/hellofs,Rxx)
title('R_{xx}')

subplot 312
plot(lagsxy/hellofs,Ryx,'o','displayname','Rxy')
hold on
plot(lagsxy/hellofs,Ryxtest,'displayname','conv(Rxx,h)')
hold off
title('R_{xy}')
legend

subplot 313
plot(lagsyy/hellofs,Ryy,'o','displayname','Ryy')
hold on
plot(lagsyy/hellofs,circshift(Ryytest,-hellofs*delay),'displayname','conv(conv(Rxx,h),h)') % I dont't know why i had to shift the stuff here to make it match, i dont understand super well how the arrays exactly line up
hold off
title('R_{yy}')
legend




