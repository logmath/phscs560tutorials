%% get a transfer function using a chirp
% Measuring transfer functions and impulse responses  section 2.3
% non-periodic ssweeps
clear;close all;clc;
fs = 44100;
fmaxc = 10000;
[Blow,Alow]=butter(8,fmaxc/(fs/2)); % lowpass filter w/corner at 10 K
dt = 1/fs;
T = 10;
N =fs*T; 
t = (0:N-1).'*dt;
xc = chirp(t,10,T,fmaxc);
xn = filter(Blow,Alow,randn(size(xc)));
pref = 20e-6;

%% make impulse response
h = randn(size(t)).*exp(-10/T*t);
h(1) = 1;
h = filter(Blow,Alow,1/10*h); % low pass filter
% h = 50*h/sum(h.^2); % make H small so the output isnt way louder than the input
figure(1)
subplot(2,1,1)
plot(t,h,'displayname','Actual IR')
legend
%% output
yc = conv(xc,h);
yn = conv(xn,h);
ty = (0:(2*N-2))*dt;
SNR = 80; % SNR in dB
noisec = rms(yc)/(10^(SNR/20))*randn(size(yc));
noisen = rms(yn)/(10^(SNR/20))*randn(size(yn)); %add noise on the output 
recordedyc =filter(Blow,Alow,yc + noisec);

figure(4)

yyaxis left
plot(t,xc,'displayname','input')
yyaxis right
plot(ty,recordedyc,'displayname','output with noise')

legend
%% calculate h from chirip
padx = [xc; zeros(length(yc)-length(xc),1)];
padt = (0:length(padx)-1)*dt;
dfpad = fs/length(yc);
Npad = length(yc);
fmaxpad = floor(Npad/2)*dfpad;
fpad = 0:dfpad:fmaxpad;

output = 1/Npad*fft(recordedyc);
outlevel = 20*log10(abs(output(1:floor(Npad/2)+1))/sqrt(dfpad)/pref); % make them spectral densities
input = 1/Npad*fft(padx);
inlevel = 20*log10(abs(input(1:floor(Npad/2)+1))/sqrt(dfpad)/pref);
RTF = output./input;
RTFlevel = 20*log10(abs(RTF(1:floor(Npad/2)+1))); % RTF should not be reference to pref, it should be reference to 1
RIR = ifft(RTF,'symmetric');

figure(2)
subplot(3,1,1)
semilogx(fpad,outlevel,'displayname','output spectrum')
hold on
semilogx(fpad,inlevel,'displayname','input spectrum')
legend
ylim([0 100])

figure(2)
subplot(3,1,2)
semilogx(fpad,RTFlevel,'displayname','RTF from chirp')
legend

figure(1)
subplot(2,1,1)
hold on
plot(padt,RIR,'--','displayname','IR from fft of chirp')
xlim([0 1])
legend
%% calculate h from Gxy/Gxx for chirp
 [f,Gxx,Gyy,Gxy,coh,OASPL,epsilon]= spec(fs,padx,recordedyc,'ns',2^15);
 RTFautochirp = Gxy./Gxx;
 
 figure(2)
 subplot(3,1,1)
 hold on
semilogx(f,10*log10(Gxx/pref^2),'displayname','Gxx')
semilogx(f,10*log10(Gyy/pref^2),'displayname','Gyy')
 xlim([10 fs/2])
 
 figure(2)
 subplot(3,1,2)
 hold on
 semilogx(f,20*log10(abs(RTFautochirp)),'displayname','RTF from autospec of chirp')
xlim([10 fs/2])

IRautochirp = ifft([RTFautochirp; flipud(RTFautochirp)],'symmetric');

figure(1)
tauto = (0:length(IRautochirp)-1)*1/fs;
subplot(2,1,2)
plot(tauto,IRautochirp,'displayname','IR from autospectra of chirp')
legend

%% calculate h from white noise using simple fourier transform
padx = [xn; zeros(length(yn)-length(xn),1)];
recordedyn =filter(Blow,Alow,yn + noisen); 
RTF = fft(recordedyn)./fft(padx);
RTFlevel = 20*log10(abs(RTF(1:floor(Npad/2)+1)));
RIR = ifft(RTF);

figure(2)
subplot(3,1,2)
semilogx(fpad,RTFlevel,'displayname','RTF from white noise fft')

figure(1)
subplot(2,1,1)
hold on
plot(padt,RIR,'displayname','IR from white noise fft')
legend
xlim([0 T])

%% calculate h from Gxy/Gxx for white noise
 [f,Gxx,Gyy,Gxy,coh,OASPL,epsilon]= spec(fs,padx,recordedyn,'ns',2^15);
 RTFautonoise = Gxy./Gxx;
 
 figure(2)
 subplot(3,1,2)
 hold on
 semilogx(f,20*log10(abs(RTFautonoise)),'displayname','RTF from autospec of white noise')
xlim([10 fs/2])
legend
IRautonoise = ifft([flipud(RTFautonoise(2:end)); RTFautonoise(1:end-1)],'symmetric');
ylim([-50 50])

subplot(3,1,3)
semilogx(f,10*log10(Gxx/pref^2),'displayname','Gxx')
hold on
semilogx(f,10*log10(Gyy/pref^2),'displayname','Gyy')
legend
ylim([0 100])
xlim([10 fs/2])

figure(1)
tauto = (0:length(IRautonoise)-1)*1/fs;
subplot(2,1,1)
hold on
plot(tauto,IRautonoise,'displayname','IR from autospectra of white noise')
legend

% %% periodic chirp
% [B,A]=butter(2,[100,1000]/(fs/2)); % make a band-pass filter for the impulse response of the mic
% IRmic = impz(B,A);
% output =conv(xc,IRmic);
% pc =[];
% for numreps = 1:20
%     pc = [pc; chir];
%     ypc = [ypc; output + rms(output)/(10^(SNR/20))*randn(size(output))];
% end
% pc = pc.';
% ypc = ypc.';
% df = 1/T;
% fmaxperiodic = floor(length(ypc)/2)*df;
% f = 0:df:fmaxperiodic-df;
% 
% padpc = [pc; zeros(length(ypc)-length(pc),numreps)];
% 
% 
% %% straight fourier transform
% Hmicstriaght= fft(reshape(ypc,[],1))./fft(reshape(padpc,[],1));
% Nbig = length(Hmicstriaght);
% dfbig = fs/Nbig;
% fbig = (0:floor(Nbig/2)-1)*dfbig;
% figure(3)
% semilogx(fbig,20*log10(sqrt(2)*abs(Hmicstriaght(1:floor(Nbig/2)))/pref),'displayname','fft of whole thing')
% title('periodic chirp')
% legend
