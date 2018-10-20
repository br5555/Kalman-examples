function Correlated(M, MFilter)

% Kalman filter simulation using correlated process and measurement noise.
% This illustrates the improvement in filter results that can be attained
% when the correlation is taken into account.
% M = correlation between process noise and measurement noise.
% MFilter = value of M used in Kalman filter.

kf = 50; % number of time steps in simulation
phi = 0.8; % system matrix
H = 1; % measurement matrix
Q = 1; % process noise covariance
R = 0.1; % measurement noise covariance

% Compute the eigendata of the covariance matrix so that we can simulate correlated noise.
Q1 = [Q M; M' R];
[d, lambda] = eig(Q1);
if (lambda(1) < 0) || (lambda(2) < 0)
    disp('Q1 is not positive semidefinite');
    return;
end
ddT = d * d';
if (norm(eye(size(ddT)) - ddT) > eps)
   disp(['d is not orthogonal. d = ',d]);
   return;
end
   
x = 0; % initial state
xhatplus = x; % initial state estimate
Pplus = 0; % initial uncertainty in state estimate

xArray = [];
xhatArray = [];
KArray = [];
PArray = [];
zArray = [];

randn('state', 0); % initialize random number generator

for k = 1 : kf
   % Generate correlated process noise (w) and measurement noise (n)
   v = [sqrt(lambda(1,1))*randn; sqrt(lambda(2,2))*randn];
   Dv = d * v;
   w = Dv(1);
   n = Dv(2);
   % Simulate the system dynamics and the measurement
   x = phi * x + w;
   z = H * x + n;
   % Simulate the Kalman filter
   Pminus = phi * Pplus * phi' + Q;
   K = inv(H * Pminus * H' + H * MFilter + MFilter' * H' + R);
   K = (Pminus * H' + MFilter) * K;
   xhatminus = phi * xhatplus;
   xhatplus = xhatminus + K * (z - H * xhatminus);
   Pplus = Pminus - K * (H * Pminus + MFilter');
   % Save data for plotting
   xArray = [xArray x];
   xhatArray = [xhatArray xhatplus];
   KArray = [KArray K];
   PArray = [PArray Pplus];
   zArray = [zArray z];
end

% Plot the data
k = 1 : kf;
close all;

figure;
plot(k, zArray - xArray, 'r-', k, xhatArray - xArray, 'b:');
set(gca,'FontSize',12); set(gcf,'Color','White');
title(['M = ',num2str(M),', MFilter = ',num2str(MFilter)]);
xlabel('Time');
legend('Measurement Error', 'Estimation Error');

figure;
plot(k, KArray);
set(gca,'FontSize',12); set(gcf,'Color','White');
title(['Kalman Gain for M = ',num2str(M),', MFilter = ',num2str(MFilter)]);
xlabel('Time');

figure;
plot(k, PArray);
set(gca,'FontSize',12); set(gcf,'Color','White');
title(['Estimation Error Covariance for M = ',num2str(M),', MFilter = ',num2str(MFilter)]);
xlabel('Time');

% Compute statistics
err = zArray - xArray;
err = norm(err)^2 / kf;
disp(['Measurement Error Variance = ',num2str(err)]);
err = xhatArray - xArray;
err = norm(err)^2 / kf;
disp(['Estimation Error Variance = ',num2str(err)]);
disp(['Analytical Variance = ',num2str(Pplus)]);
