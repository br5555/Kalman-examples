function DiscreteKFEx2Plot

% Create plots for example in Discrete Kalman Filter chapter.

[x1Array, xhat100, K100] = DiscreteKFEx2(1.00);
[x1Array, xhat010, K010] = DiscreteKFEx2(0.10);
[x1Array, xhat001, K001] = DiscreteKFEx2(0.01);
[x1Array, xhat000, K000] = DiscreteKFEx2(0.00);

t = 0 : 50;
close all;
figure; hold on;
plot(t, x1Array, 'k');
plot(t, xhat100, 'b*--');
plot(t, xhat010, 'rd-');
plot(t, xhat001, 'mo:');
plot(t, xhat000, 'k+');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('time');
legend('true state', 'xhat (Q = 1)', 'xhat (Q = 0.1)', 'xhat (Q = 0.01)', 'xhat (Q = 0)');
set(gca,'box','on');

figure; hold on;
plot(t, K100, 'b*--');
plot(t, K010, 'rd-');
plot(t, K001, 'mo:');
plot(t, K000, 'k+');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('time'); 
legend('Q = 1', 'Q = 0.1', 'Q = 0.01', 'Q = 0');
set(gca,'box','on');
