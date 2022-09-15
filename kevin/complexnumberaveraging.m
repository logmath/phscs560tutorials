clear;close all;clc;
z = 10*(1/sqrt(2) +1/sqrt(2)*1j);
N =100;
er = abs(z)*(1-.1*randn(N,1)).*exp(1j*(angle(z)+0.5.*pi/2*randn(N,1))); %number with random amplitude and phase noise
guess = mean(abs(er))*exp(1j*mean(angle(er)));
errors = [zeros(N,1) er];

plot([0 z],'r-','linewidth',5)
hold on
for n = 1:N
    plot(errors(n,:),'k')
end
plot([0 guess],'g--','linewidth',5)
plot(mean(errors),'b:','linewidth',5)
legend('actual number','guess','number with random error added')
grid on
axis image

