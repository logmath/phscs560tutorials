function H = TransferFunction(f,B,A)
arguments
    f (1,:) 
    B (1,:) = 1;  % numerator   
    A (1,:) = 1;  % denominator
end
m = (1:length(B)).'; %column vector of powers
n = (1:length(A));
H = (B*exp(1j*2*pi*f).^-(m-1))./((A*exp(1j*2*pi*f).^-(n-1)));  
end