function HinfContEx1a

% Continuous time H-infinity filtering example, theta = 7/16

P0 = 1;
c = (9*P0+4)/(40*P0-160);
t = 0 : 0.01 : 10;
P = (4 + 160 * c * exp(5*t/2)) ./ (-9 + 40 * c * exp(5*t/2));

close all;
plot(t,P);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time'); ylabel('Riccati equation solution P');
axis([0 10 1 5]);
