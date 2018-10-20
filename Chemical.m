function Chemical

% Optimal State Estimation, by Dan Simon

N = 51;
x = [10; 5];
xhat = [8; 7];
P = eye(2);
R = 0.01;

xhatArr = xhat;
PArr(:,:,1) = P;
for k = 1 : N-1
   H = [1 0.99^(k-1)];
   y = H * x + sqrt(R) * randn;
   K = P * H' * inv(H * P * H' + R);
   xhat = xhat + K * (y - H * xhat);
   P = P - K * H * P;
   % Save data in arrays.
   xhatArr = [xhatArr xhat];
   PArr(:,:,k+1) = P;
end

% Plot the results.
close all; % close all figures
figure; % open a new figure
set(gcf,'Color','White');
k = 0 : N-1;
subplot(2,1,1); hold on;
plot(k, xhatArr(1,:), 'r-', 'LineWidth', 2);
plot(k, xhatArr(2,:), 'b:', 'LineWidth', 2);
set(gca,'FontSize',12);
ylabel('estimates');
legend('x_1', 'x_2');

subplot(2,1,2); hold on;
plot(k, squeeze(PArr(1,1,:)), 'r-', 'LineWidth', 2);
plot(k, squeeze(PArr(2,2,:)), 'b:', 'LineWidth', 2);
set(gca,'FontSize',12);
ylabel('variances');
legend('P(1,1)', 'P(2,2');
xlabel('time step');
