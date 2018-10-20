function DiscreteKFAlt

% Simulate a discrete-time scalar Kalman filter.

tf = 10; % final time

F = 1; % state transition matrix
H = 1; % measurement matrix
Q = 1; % process noise covariance
R = 1; % measurement noise covariance

x = 0; % initial true state
xhat = 0; % initial estimate of x
Pplus = 1; % initial estimation error covariance

xArray = [];
xhatArray = [];
KArray = [];
PArray = [];

randn('state', 2);

for t = 0 : tf
   % Simulate the system
   x = F * x + sqrt(Q) * randn;
   y = H * x + sqrt(R) * randn;
   % Kalman filter
   Pminus = F * Pplus * F' + Q;
   K = Pminus * H' * inv(H * Pminus * H' + R);
   K = (1 + sqrt(5)) / (3 + sqrt(5)); % steady state Kalman gain
   xhat = F * xhat;
   xhat = xhat + K * (y - H * xhat);
   Pplus = (1 - K * H) * Pminus * (1 - K * H)' + K * R * K';
   % Save data for later
   xArray = [xArray; x];
   xhatArray = [xhatArray; xhat];
   KArray = [KArray; K];
   PArray = [PArray; Pminus];
end

% Plot results
close all;
t = 0 : tf;

figure;
plot(t, xArray, 'r-', t, xhatArray, 'b--');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('time');
legend('true state', 'estimated state');

figure;
plot(t, KArray, 'r-', t, PArray, 'b--');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('time');
legend('Kalman Gain', 'Estimation Covariance');

disp(['RMS estimation error = ', num2str(std(xArray - xhatArray))]);