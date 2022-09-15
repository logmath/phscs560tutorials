% BeamPatterns
clear;close all;clc
N =10;
theta = -90:.1:90;
for kd =0.1:.1:10
    H = 1/N*abs(sin(N/2*kd*sind(theta))./sin(1/2*kd*sind(theta)));
    figure(1)
    subplot 211
    plot(theta,H)
    xlabel('angle (deg)')
    ylabel('|H|')
    title(['line array N = ' num2str(N) ', kd = ' num2str(kd)])
    
    subplot 212
    polarplot((theta)*pi/180,20*log10(H))
    rlim([-20 1])
    thetalim([-90 90])
    set(gca,'ThetaZeroLocation','top')
    rticks(-20:10:0)
    
    figure(2) % baffled circular piston with radius d
    Hbcp = abs(2*besselj(1,kd*sind(theta))./(kd*sind(theta)));
    subplot 211
    plot(theta,Hbcp)
    xlabel('angle (deg)')
    ylabel('|H|')
    title(['Baffled Circular Piston ka = ' num2str(kd)])
    ylim([0 1])
    
    subplot 212
    polarplot((theta)*pi/180,20*log10(Hbcp))
    rlim([-20 1])
    thetalim([-90 90])
    set(gca,'ThetaZeroLocation','top')
    rticks(-20:10:0)
    pause(.1)
    
    
end