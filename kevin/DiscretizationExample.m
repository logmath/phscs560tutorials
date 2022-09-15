%% discretization example
clear;close all;clc;
%% make waveform
T = 1.024;
dt = 0.0001;
t = (0:dt:T).';
ap = [2 3 1.5 2.5];
phis = [0 pi/4 2*pi/3 pi];
effs = [31.25 125 187.5 375];
a = sum(real(ap.*exp(1j*(phis+2*pi*effs.*t))),2);
% a = a/4;
%% sample it
fs = 1000;
dtnew = 1/fs;
skip = dtnew/dt;
N = 1024;
tsampled =t(1:skip:end); 
asampled = a(1:skip:end);
plot(t,a,'-',tsampled,asampled,'o')
xlim([0 .031]);

%% simulate discretization
for M =[4 8 16 32 64] % for several bit depths
Vpp = 18;
Q = Vpp/2^M; % this is the "quantization level"

adis = Q*round(asampled./Q);
error = adis-asampled;

if M==4 % only plot one example
    hold on
    bits = -Vpp/2:Q:Vpp/2;
for kk = 1:length(bits)
    plot([0 .031],[bits(kk) bits(kk)],'g--')
end
    stem(tsampled,error)
end
%% SNR
errorvar = 1/12*(Vpp/2^M)^2;
varsinusoidfullbitdepth = Vpp^2/8;
SNRest =10*log10(varsinusoidfullbitdepth/errorvar)
altSNRest = 1.76+20*M*log10(2);
SNR = 10*log10(var(a)/var(error))

% %% 
% spectra = abs(fftshift(fft(adis))).^2;
% f  = 0:1/T:length(t)/2*1/T;
pause(1)
end