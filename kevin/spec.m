function [f,Gxx,Gyy,Gxy,coh,OASPL,epsilon] = spec(fs,x,y,options)
% call  [f,Gxx,Gyy,Gxy,coh,epsilon] = spec(fs,x,y,options)
% This program calulates the autospectra of signals x and y as well as the
% crossspectrum and coherence between them. 
% A Hanning windowing is used with 50% overlap. Per Bendat and Piersol, 
% Section 11.6.3
% they are scaled by the mean-square value of the window for overall
% amplitude scaling purposes.
%% Inputs:
% x and y   the one-dimentional signals 
% fs        sampling frequency 
%% Options
% ns        number of samples per block. default is 2^15 
% density   a boolean to decide if you want to divide the spectra by the 
%           bin width to make them auto/cross sepectral densities
% N         how many samples of the signals to take, default is the whole
%           signal, choose a power of 2 for faster ffts
%% Outputs
% f = output frequency array
% Gxx and Gyy are the autospectrum (or autospectral densities)
% Gxy = cross spectrum (or density)
% coh = the coherence
% OASPL(1) = the overall SPL calculated for x in frequency domain
% OASPL(2) = the overall SPL calculated for x in time domain
% OASPL(3) = the overall SPL calculated for y in frequency domain
% OASPL(4) = the overall SPL calculated for y in time domain
% epsilon = the normalized random error of the coherence calculation
% 
% By Kevin Leete, taken heavily from Kent Gee and Alan Wall's autospec and
% crosspec functions

arguments
    fs (1,1)
    x  (:,1)
    y (:,1) = x;
    options.ns (1,1) = 2^15;
    options.density (1,1) = 1;
    options.N =length(x); %choose number of samples to be a power of 2 for speed
end
x = x-mean(x);% zero mean the data
y = y-mean(y);
x = x(1:options.N); %truncate recording to be right length
y = y(1:options.N);
df = fs/options.ns; % frequency resolution is dependant on the sampling rate and
            % the length of the blocks
f = df * (0:(floor(options.ns/2)-1)); % single sided frequency array
ww = hann(options.ns); % generates hanning window
W = mean(ww.*conj(ww)); % used to scale the PSD
numblocks = floor(2*options.N/options.ns-1); % number of blocks (with 50% overlap)
epsilon = 1/sqrt(numblocks);
if epsilon >.1
    disp(sprintf(['for ' num2str(numblocks) ' blocks, ' char(cellstr(char(949))) ' = ' num2str(round(epsilon,3)) '.\n Consider decreasing ns or increasing N']))
    %bendat and piersol give \e = 1/sqrt(numblocks) as the normalized rms
    %error
end
xblocks=zeros(options.ns,numblocks);
yblocks=zeros(options.ns,numblocks); 
for k = 1:numblocks
    if k==1
        start(k) = 1;
        fin(k) = options.ns;
    else
        start(k) = (k-1)*options.ns/2;
        fin(k) = (k+1)*options.ns/2-1;
    end    
    xblocks(:,k) = ww.*x(start(k):fin(k));% divide signal into blocks %window each block
    yblocks(:,k) = ww.*y(start(k):fin(k));
end

X = sqrt(1/W)*1/options.ns*fft(xblocks,options.ns); % do fourier transform and scale for energy conservation
Y = sqrt(1/W)*1/options.ns*fft(yblocks,options.ns);
Xss = sqrt(2) * X(1:options.ns/2,:); % make single sided
Yss = sqrt(2) * Y(1:options.ns/2,:);
Gxx = mean(Xss.*conj(Xss),2);  % make it a autospectrum (Pa^2) by averaging squared values over blocks
Gyy = mean(Yss.*conj(Yss),2); 
Gxy = mean(Xss.*conj(Yss),2);  %calculate the cross spectrum
coh = abs(Gxy).^2./(Gxx.*Gyy); % calculate coherence
OASPL(1) = 10*log10(sum(Gxx)./(20e-6)^2); % OASPL of x from frequency domain
OASPL(2) = 20*log10(rms(x)/20e-6); % OASPL of x from time domain
OASPL(3) = 10*log10(sum(Gyy)./(20e-6)^2); %OASPL of y from frequency domain
OASPL(4) = 20*log10(rms(y)/20e-6); % OASPL of y from time domain
if options.density %make them auto and cross spectral densities (Pa^2/Hz)
    Gxx = Gxx/df;
    Gyy = Gyy/df;
    Gxy = Gxy/df;
end

end









