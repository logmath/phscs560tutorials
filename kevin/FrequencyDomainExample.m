clear; close all; clc;

%% make signal 
fs=44000; % sample rate
dt=1/fs; % sample period
T=5; % total sample length
t=0:dt:(T-dt); % time array
N = length(t); % number of samples
[B,A]=butter(2,[500,5000]/(fs/2)); % make a band-pass filter 
noise = 3*filter(B,A,randn(1,length(t))); % make band-limited gaussian white noise
tonefreq = 4000;
tone = sin(2*pi*tonefreq*t);
x= tone+noise; % add the noise to a sine wave at 4khz
%% play sound
% soundsc(x,fs) % play the sound

%% plot
figure(1)
plot(t,x) % plot the signal
xlabel('time (s)')
ylabel('amplitude (arbitrary)')

%% fourier transform
figure(2)
subplot(1,2,1)
xhat = 1/N*fft(x); % MATLAB's fft algorithm doesn't multiply by 1/N for whatever reason
plot(abs(xhat))
xlabel('sample')
ylabel('|xhat|')
title('output of fft(x)')

subplot(1,2,2)
xhatshift = fftshift(xhat);
plot(abs(xhatshift))
xlabel('sample')
ylabel('|xhat|')
title('fftshift(fft(x))')

%% define double sided frequency array
df = 1/T;
fmax = floor(N/2)*df; % for even N, this is the Nyquist rate. for odd N, this is 0.5*df less than the Nyquist rate

% The DFT is N-periodic. The N+1th sample (if you went there) is the
% same as the first sample (the Zero frequency bin).

if mod(N,2)==0 % is N is even
    f = -fmax:df:fmax-df;
% for even N, fmax is the nyquist rate. 
% For an even numbered N, the first sample is the 0 frequency and the (N/2
% +1)th sample is the sample at the Nyquist rate (by whatever standard that
% defines the FFT, the Nyquist rate is placed at the end of the negative 
% side of the double sided array, so the positive side only goes up one bin
% less than the nyquist.
% The Nth sample is then the sample right before the 0 frequency 
% (or equivalently the sampling frequency) again. because of periodicity 
% the Nth sample is the same as the first negative frequency, the N-1th 
% sample is the same as the second negative frequency, and so forth. 
else % N is odd
    f = -fmax:df:fmax;
% for odd N fmax is not the Nyquist, its half a bin below the Nyquist
% for odd N the first sample is the zero frequency and the floor(N/2)th
% sample represents the frequcny half a bin below the Nyquist. The
% floor(N/2)+1 to Nth frequencies then, by periodicity, are the same as the
% -floor(N/2)th ascending to the sample right before the zero frequency
% again.
end
unshiftedf = (0:N-1)*df;
fss = 0:df:fmax;

%% plot fourier transform with appropriate frequencies
figure(3)
subplot(1,2,1)
plot(unshiftedf/1000,abs(xhat))
xlabel('frequency(kHz)')
ylabel('|xhat|')
title('raw output of fft')

subplot(1,2,2)
plot(f/1000,abs(xhatshift))
xlabel('frequency(kHz)')
ylabel('|xhat|')
title('shifting to be zero-centered')

%% plot single sided on a loglog scale
xhatss = sqrt(2)*xhat(1:floor(N/2)+1); % multiply the complex values by sqrt(2), so that when its squared its the same as being multiplied by 2, because we lost half the energy by only going single sided
figure(4)
loglog(fss,abs(xhatss))
title('single sided fft on a loglog frequency scale')

%% plot single sided spectrum on semilogx scale with the amplitues converted to level
figure(5)
pref=20e-6;
spl = 20*log10(abs(xhatss)/pref);
semilogx(fss,spl)
xlim([20 20000])
ylabel('SPL (dB re 20\muPa)')
xlabel('Frequency (Hz)')

%% calculate OASPL
OASPLtime= 20*log10(rms(x)/pref)
OASPLfreq = 10*log10(sum(abs(xhat).^2)./pref^2)
OASPLfreqss = 10*log10(sum(abs(xhatss).^2)./pref^2) 
OASPLtonetime =  20*log10(rms(tone)/pref) % calculate the level of just the tone
OASPLnoisetime = 20*log10(rms(noise)/pref) % calculate level of just the noise

%% window
w = hann(N);
xwindowed = w.'.*x;
figure(6)
plot(xwindowed)
title('windowed x')
OASPLwindowed = 20*log10(rms(xwindowed)/pref)
energylost2window = OASPLtime - OASPLwindowed

Xwindowed = abs(1/N*fft(xwindowed)).^2;
figure(6)
semilogx(fss,10*log10(2*Xwindowed(1:floor(N/2)+1)/pref^2))
title('fft of windowed x')
xlim([20 20000])


%% Autospec
for ns = [2^15 2^13 2^10]
figure(7)
[f,Gxx,OASPLauto]= autospec(x,fs,ns);
semilogx(f,10*log10(Gxx/pref^2),'displayname',['Block size = ' num2str(ns)])
hold on
xlim([20 20000])
title('autospectral density ')
xlabel('f (hz)')
ylabel('SPL (dB re 20\muPa)')
legend('location','south')

figure(8)
[f,Gxx,OASPLauto]= autospec(x,fs,ns,N,1);
semilogx(f,10*log10(Gxx/pref^2),'displayname',['Block size = ' num2str(ns)])
xlim([20 20000])
hold on
title('autospectrum ')
xlabel('f (hz)')
ylabel('SPL (dB re 20\muPa)')
legend('location','south')
end

%% autospec from correlation
[Rxx,lags] = xcorr(x,'biased');% need biased to get the level right, it basically just divides by N
Nrxx = length(Rxx);
df = fs/Nrxx;
fssrxx = (0:floor(Nrxx/2)-1)*df;
Sxx = abs(1/Nrxx*fft(Rxx));
Gauto = 2*Sxx(1:floor(Nrxx/2));

figure(8)
hold on
h = semilogx(fssrxx,10*log10(Gauto/pref^2),'displayname','calculated from Rxx');
uistack(h,'bottom')

figure(7)
hold on
h= semilogx(fssrxx,10*log10(Gauto/df/pref^2),'displayname','calculated from Rxx');
uistack(h,'bottom')

%% get amplitude of sine wave from autospec
window = hann(2^15); % make a window for the 2^15 bin width
W = mean(abs(window).^2); % mean squared value of the window
w = mean(w); % mean value of the window
[fme,Gxxme,~,~,~,OASPLme] = spec(fs,x/w); 
 [~,ind]= min(abs(fme-tonefreq)); % find the index where the tone frequency is
 fme(ind) % this is the frequency bin closest to the tone frequency
sqrt(Gxxme(ind)) % this is the amplitude estimation of the sine wave created in line 12
 