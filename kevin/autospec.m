function [f,Gxx,OASPL] = autospec(x,fs,ns,N,unitflag)
arguments
    x (:,1) % signal
    fs (1,1) % sample rate
    ns (1,1) = 2^15; % block size
    N (1,1) = 2^floor(log2(length(x))); % length of array to grab from x
    unitflag = 0 % flag to make it autospectrum, default is density
end
% This program calulates the autospectral density or autospectrum and the OASPL of a signal.
% Hanning windowing is used, with 50% overlap. Per Bendat and Piersol, Gxx 
% is scaled by the mean-square value of the window to recover the correct OASPL.
%
%   call [f,Gxx,OASPL] = autospec(x,fs,ns,N);
% 
%   Outputs: 
%   f = frequency array for plotting
%   Gxx = Single-sided autospectrum or autospectral density, depending on unitflag
%   OASPL = Overall sound pressure level
%
%   Inputs:
%   x = time series data.
%   fs = sampling frequency
%   ns = number of samples per block.  Default is 2^15 if not specified.
%   N = total number of samples.  If N is not an integer multiple of ns, 
%       the samples less than ns in the last block are discarded.  Default   
%       is nearest lower power of 2 if not specified.
%   unitflag = 1 for autospectrum, 0 for autospectral density.  Default is
%   autospectral density
%
%   Authors: Kent Gee and Alan Wall

x = x(1:N);

%single-sided FREQUENCY ARRAY
f = fs*(0:ns/2-1)/ns;
df = f(3) - f(2);   %Width of frequency bins.

%Enforce zero-mean
x = x-mean(x);

%HANNING WINDOW
ww = hann(ns);
W = mean(ww.*conj(ww)); %Used to scale the ASD for energy conservation

%SPLITS DATA INTO BLOCKS
% Divides total data set into blocks of length ns with 50% overlap, and
% windowed.

blocks(1,:) = ww.*x(1:ns);

for k = 2:floor(2*N/ns-1)  
    blocks(k,:) = ww.*x((k-1)*ns/2:(k+1)*ns/2-1); 
end                                                 

%COMPLEX PRESSURE AMPLITUDE
X = sqrt(2*df/ns/fs/W)*fft(blocks,ns,2);    %Correct values for a single-
                                            %sided fft.   
Xss = X(:,1:ns/2);  %Takes first ns/2 points to make it single-sided. 
                    %Units are Pa.

%AUTOSPECTRAL DENSITY
    Gxx =(1/df)*mean(conj(Xss).*Xss,1); %Units are Pa^2/Hz
    

%OVERALL SOUND PRESSURE LEVEL
if nargout > 2
    OASPL = 20*log10(sqrt(sum(Gxx*df))/2e-5); 
end

%AUTOSPECTRUM SCALING?
    Gxx=Gxx*df^unitflag; % units are Pa^2


end






