clear;close all;clc;
% Lab 5 long chirp vs 20s white noise
dataloc = 'C:\Users\kevin\Box\Lab 5 560 Data\';
pref = 20e-6;
%% windscreen on
onref = binfileload(dataloc,'ID',6,0);
on = binfileload(dataloc,'ID',6,1);
fs = 102400;
[f,Gxx,Gyy,Gxy,coh,OASPL,epsilon] = spec(fs,onref,on);
Hon = Gxy./Gxx;
figure
semilogx(f,coh)
figure
%% windscreen off
offref = binfileload(dataloc,'ID',7,0);
off = binfileload(dataloc,'ID',7,1);
[f,Gxx,Gyy,Gxy,coh,OASPL,epsilon] = spec(fs,offref,off);
Hoff = Gxy./Gxx;
%% Influence of windscreen
H = Hon./Hoff;
semilogx(f,20*log10(abs(Hoff)),f,10*log10(abs(Hon).^2),f,20*log10(abs(H)));
legend({'off','on','H of windscreen'})
xlim([10 50000])

%% only grab the single chirp
onstart = 135600;
offstart = 129800;
clsoeon = on((0:((3*fs)-1))+onstart);
onref = onref((0:((3*fs)-1))+onstart);
off = off((0:((3*fs)-1))+offstart);
offref = offref((0:((3*fs)-1))+offstart);
%% do the ffts
Hon1 = fft(on)./fft(onref);
Hoff1 = fft(off)./fft(offref);
H1 = Hon1./Hoff1;
df = 1/3;
f1 = (0:((3*fs)-1))*df;

semilogx(f1(1:3*fs),10*log10(abs(Hon1(1:3*fs))))
hold on
semilogx(f1(1:3*fs),10*log10(abs(Hoff1(1:3*fs))))
semilogx(f1(1:3*fs),10*log10(abs(H1(1:3*fs))))
legend({'on','off','windscreen'})