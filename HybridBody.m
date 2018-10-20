function [AltErr, VelErr, BallErr] = HybridBody

% Hybrid extended Kalman filter example.
% Track a body falling through the atmosphere.
% Outputs are:
%   AltErr = RMS altitude estimation error
%   VelErr = RMS velocity estimation error
%   BallErr = RMS ballistic coefficient estimation error

rho0 = 0.0034; % lb-sec^2/ft^4
g = 32.2; % ft/sec^2
k = 22000; % ft
R = 100; % measurement variance (ft^2)

x = [100000; -6000; 2000]; % initial state
xhat = [100010; -6100; 2500]; % initial state estimate
H = [1 0 0]; % measurement matrix

P = [500 0 0; 0 20000 0; 0 0 250000];
T = 0.5; % measurement time step
tf = 16; % simulation length
dt = tf / 40000; % time step for integration
xArray = x;
xhatArray = xhat;
for t = T : T : tf
   % Simulate the system.
   for tau = dt : dt : T
      xdot(1,1) = x(2);
      xdot(2,1) = rho0 * exp(-x(1)/k) * x(2)^2 / 2 / x(3) - g;
      xdot(3,1) = 0;
      x = x + xdot * dt;
   end
   % Simulate the measurement.
   z = H * x + sqrt(R) * randn;
   % Simulate the continuous-time part of the filter.
   for tau = dt : dt : T
      xhatdot(1,1) = xhat(2);
      xhatdot(2,1) = rho0 * exp(-xhat(1)/k) * xhat(2)^2 / 2 / xhat(3) - g;
      xhatdot(3,1) = 0;
      xhat = xhat + xhatdot * dt;
      F = [0 1 0; -rho0 * exp(-xhat(1)/k) * xhat(2)^2 / 2 / k / xhat(3) ...
            rho0 * exp(-xhat(1)/k) * xhat(2) / xhat(3) ...
            -rho0 * exp(-xhat(1)/k) * xhat(2)^2 / 2 / xhat(3)^2; ...
            0 0 0];
      Pdot = F * P + P * F';
      P = P + Pdot * dt;
   end
   % Simulate the discrete-time part of the filter.
   K = P * H' * inv(H * P * H' + R);
   xhat = xhat + K * (z - H * xhat);
   P = (eye(3) - K * H) * P * (eye(3) - K * H)' + K * R * K';
   % Save data for plotting.
   xArray = [xArray x];
   xhatArray = [xhatArray xhat];
end

% Plot data
close all;
t = 0 : T : tf;

figure;
plot(t, xArray(1,:) - xhatArray(1,:));
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time (seconds)');
ylabel('Altitude Estimation Error (feet)');

AltErr = std(xArray(1,:) - xhatArray(1,:));
disp(['Hybrid EKF RMS altitude estimation error = ', num2str(AltErr)]);

figure;
plot(t, xArray(2,:) - xhatArray(2,:));
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time (seconds)');
ylabel('Velocity Estimation Error (feet/sec)');

VelErr = std(xArray(2,:) - xhatArray(2,:));
disp(['Hybrid EKF RMS velocity estimation error = ', num2str(VelErr)]);

figure;
plot(t, xArray(3,:) - xhatArray(3,:));
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time (seconds)');
ylabel('Ballistic Coefficient Estimation Error');

BallErr = std(xArray(3,:) - xhatArray(3,:));
disp(['Hybrid EKF RMS ballistic coefficient estimation error = ', num2str(BallErr)]);

figure;
plot(t, xArray(1,:)/1000);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time (seconds)');
ylabel('True Altitude (thousands of feet)');

figure;
plot(t, xArray(2,:));
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time (seconds)');
ylabel('True Velocity (feet/sec)');