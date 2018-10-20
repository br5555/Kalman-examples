function Multiple

% Multiple model Kalman filtering.
% Second order system.
% Corrected March 17, 2009, due to typos in Equations 10.35 and 10.37.

zeta = 0.1; % damping ratio
wn = [sqrt(4); sqrt(4.4); sqrt(4.8)]; % possible wn values
N = size(wn, 1); % number of parameter sets
pr = [0.1; 0.6; 0.3]; % a priori probabilities
% Compute the initial estimate of wn
wnhat = 0;
for i = 1 : N
   wnhat = wnhat + wn(i) * pr(i);
end
T = 0.1; % sample period
Qc = 1000; % continuous time process noise variance
R = diag([10 10]); % discrete time measurement noise covariance
H = eye(2); % measurement matrix
q = size(H, 1); % number of measurements
x = [0; 0]; % initial state

% Compute the alternative Lambda and phi matrices.
for i = 1 : N
   Ai = [0 1; -wn(i)^2 -2*zeta*wn(i)];
   Bi = [0; wn(i)^2];
   Fi = expm(Ai*T);
   F(:,:,i) = Fi;
   %Lambda(:,:,i) = (phii - eye(size(phii))) * inv(Fi) * Li;
   Pplus(:,:,i) = zeros(size(Fi));
   xhat(:,i) = x;
end

B = [0; wn(1)^2];
Q = B * Qc * B' * T; % discrete time process noise covariance

tf = 60; % Length of simulation
% Create arrays for later plotting
wnhatArray = [wnhat];
prArray = [pr];
for t = T : T : tf
   % Simulate the system.
   % The first parameter set is the true parameter set.
   w = sqrt(Q) * randn(2, 1);
   x = F(:,:,1) * x + w;
   y = H * x + sqrt(R) * randn(2, 1);
   % Run a separate Kalman filter for each parameter set.
   for i = 1 : N
      Pminus(:,:,i) = F(:,:,i) * Pplus(:,:,i) * F(:,:,i)';
      Pminus(:,:,i) = Pminus(:,:,i) + Q;
      K = Pminus(:,:,i) * H' * inv(H * Pminus(:,:,i) * H' + R);
      xhat(:,i) = F(:,:,i) * xhat(:,i);
      r = y - H * xhat(:,i); % measurment residual
      S = H * Pminus(:,:,i) * H' + R; % covariance of measurement residual
      pdf(i) = exp(-r'*inv(S)*r/2) / ((2*pi)^(q/2)) / sqrt(det(S));
      xhat(:,i) = xhat(:,i) + K * (y - H * xhat(:,i));
      Pplus(:,:,i) = (eye(2) - K * H) * Pminus(:,:,i) * (eye(2) - K * H)' + K * R * K';
   end
   % Compute the sum that appears in the denominator of the probability expression.
   Prsum = 0;
   for i = 1 : N
      Prsum = Prsum + pdf(i) * pr(i);
   end
   % Update the probability of each parameter set.
   for i = 1 : N
      pr(i) = pdf(i) * pr(i) / Prsum;
   end
   % Compute the best state estimate and the best parameter estimate.
   xhatbest = 0;
   wnhat = 0;
   for i = 1 : N
      xhatbest = xhatbest + pr(i) * xhat(:,i);
      wnhat = wnhat + pr(i) * wn(i);
   end
   % Save data for plotting.
   wnhatArray = [wnhatArray wnhat];
   prArray = [prArray pr];
end

close all;

t = 0 : T : tf;
figure;
plot(t, wnhatArray.^2);
title('Estimate of square of natural frequency', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');

figure;
plot(t, prArray(1,:), 'b-', t, prArray(2,:), 'k--', t, prArray(3,:), 'r:');
title('Probabilities of square of natural frequency', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
legend('Probability that \omega_n^2 = 4', 'Probability that \omega_n^2 = 4.4', 'Probability that \omega_n^2 = 4.8');
