function HybridUKF

% Hybrid extended Kalman filter example.
% Track a body falling through the atmosphere.
% This example is taken from [Jul00], which was based on [Ath68].

rho0 = 2; % lb-sec^2/ft^4
g = 32.2; % ft/sec^2
k = 2e4; % ft
R = 10^4; % measurement noise variance (ft^2)
Q = diag([0 0 0]); % process noise covariance
M = 10^5; % horizontal range of position sensor
a = 10^5; % altitude of position sensor

x = [3e5; -2e4; 1e-3];
xhat = [3e5; -2e4; 1e-3];
xhatukf = xhat;

P = diag([1e6 4e6 10]);
Pukf = P;

T = 0.5; % measurement time step
randn('state',sum(100*clock)); % random number generator seed

tf = 30; % simulation length (seconds)
dt = 0.001; % time step for integration (seconds)
xArray = x;
xhatArray = xhat;
xhatukfArray = xhatukf;
Parray = diag(P);
Pukfarray = diag(Pukf);

W = ones(6,1) / 6; % UKF weights

for t = T : T : tf
   % Simulate the system.
   for tau = dt : dt : T
      xdot(1,1) = x(2);
      xdot(2,1) = rho0 * exp(-x(1)/k) * x(2)^2 / 2 * x(3) - g;
      xdot(3,1) = 0;
      xdot = xdot + sqrt(dt * Q) * [randn; randn; randn];
      x = x + xdot * dt;
   end
   % Simulate the noisy measurement.
   z = sqrt(M^2 + (x(1)-a)^2) + sqrt(R) * randn;
   % Simulate the continuous-time part of the Kalman filter (time update).
   for tau = dt : dt : T
      xhatdot(1,1) = xhat(2);
      xhatdot(2,1) = rho0 * exp(-xhat(1)/k) * xhat(2)^2 / 2 * xhat(3) - g;
      xhatdot(3,1) = 0;
      xhat = xhat + xhatdot * dt;
      F = [0 1 0; -rho0 * exp(-xhat(1)/k) * xhat(2)^2 / 2 / k * xhat(3) ...
            rho0 * exp(-xhat(1)/k) * xhat(2) * xhat(3) ...
            rho0 * exp(-xhat(1)/k) * xhat(2)^2 / 2; ...
            0 0 0];
      H = [(xhat(1) - a) / sqrt(M^2 + (xhat(1)-a)^2) 0 0];
      Pdot = F * P + P * F' + Q * dt - P * H' * inv(R) * H * P;
      P = P + Pdot * dt;
   end
   % Simulate the discrete-time part of the Kalman filter (measurement update).
   H = [(xhat(1) - a) / sqrt(M^2 + (xhat(1)-a)^2) 0 0];
   K = P * H' * inv(H * P * H' + R);
   zhat = sqrt(M^2 + (xhat(1)-a)^2);
   xhat = xhat + K * (z - zhat);
   P = (eye(3) - K * H) * P * (eye(3) - K * H)' + K * R * K';
   % Generate the UKF sigma points.
   [root,p] = chol(3*Pukf);
   
   for i = 1 : 3
       sigma(:,i) = xhatukf + root(i,:)';
       sigma(:,i+3) = xhatukf - root(i,:)';
   end
   for i = 1 : 6
       xbreve(:,i) = sigma(:,i);
   end
   % UKF time update
   for i = 1 : 6
       for tau = dt : dt : T
          xbrevedot(1,1) = xbreve(2,i);
          xbrevedot(2,1) = rho0 * exp(-xbreve(1,i)/k) * xbreve(2,i)^2 / 2 * xbreve(3,i) - g;
          xbrevedot(3,1) = 0;
          xbreve(:,i) = xbreve(:,i) + xbrevedot * dt;
      end
  end
  xhatukf = zeros(3,1);
  for i = 1 : 6
      xhatukf = xhatukf + W(i) * xbreve(:,i);
  end
  Pukf = zeros(3,3);
  for i = 1 : 6
      Pukf = Pukf + W(i) * (xbreve(:,i) - xhatukf) * (xbreve(:,i) - xhatukf)';
  end
  Pukf = Pukf + Q;
  % UKF measurement update
  for i = 1 : 6
      zukf(:,i) = sqrt(M^2 + (xbreve(1,i)-a)^2);
  end
  zhat = 0;
  for i = 1 : 6
      zhat = zhat + W(i) * zukf(:,i);
  end
  Py = 0;
  Pxy = zeros(3,1);
  for i = 1 : 6
      Py = Py + W(i) * (zukf(:,i) - zhat) * (zukf(:,i) - zhat)';
      Pxy = Pxy + W(i) * (xbreve(:,i) - xhat) * (zukf(:,i) - zhat)';
  end
  Py = Py + R;
  Kukf = Pxy * inv(Py);
  xhatukf = xhatukf + Kukf * (z - zhat);
  Pukf = Pukf - Kukf * Py * Kukf';      
   
   % Save data for plotting.
   xArray = [xArray x];
   xhatArray = [xhatArray xhat];
   xhatukfArray = [xhatukfArray xhatukf];
   Parray = [Parray diag(P)];
   Pukfarray = [Pukfarray diag(Pukf)];
end

close all;
t = 0 : T : tf;
figure; 
semilogy(t, abs(xArray(1,:) - xhatArray(1,:)), 'b'); hold;
%plot(t, sqrt(Parray(1,:)), 'b--');
semilogy(t, abs(xArray(1,:) - xhatukfArray(1,:)), 'r:');
%plot(t, sqrt(Pukfarray(1,:)), 'r--');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('Position Estimation Error');
legend('Kalman filter', 'Unscented filter');

figure; 
semilogy(t, abs(xArray(2,:) - xhatArray(2,:)), 'b'); hold;
%plot(t, sqrt(Parray(2,:)), 'b--');
semilogy(t, abs(xArray(2,:) - xhatukfArray(2,:)), 'r:');
%plot(t, sqrt(Pukfarray(2,:)), 'r--');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('Velocity Estimation Error');
legend('Kalman filter', 'Unscented filter');

figure; 
semilogy(t, abs(xArray(3,:) - xhatArray(3,:)), 'b'); hold;
%plot(t, sqrt(Parray(3,:)), 'b--');
semilogy(t, abs(xArray(3,:) - xhatukfArray(3,:)), 'r:');
%plot(t, sqrt(Pukfarray(3,:)), 'r--');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('Ballistic Coefficient Estimation Error');
legend('Kalman filter', 'Unscented filter');

figure;
plot(t, xArray(1,:));
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('True Position');

figure;
plot(t, xArray(2,:));
title('Falling Body Simulation', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('True Velocity');
