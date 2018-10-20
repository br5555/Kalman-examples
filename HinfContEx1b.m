function HinfContEx1b

% Continuous time H-infinity filtering example
% INPUTS:
%   theta = H-infinity performance bound

close all

[t, z0] = HinfFilter(0, 7/16);
[t, z1] = HinfFilter(1, 7/16);
[t, zk] = HinfFilter(0, 0);

figure;
plot(t,z0,'r-', t,z1,'b:');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time'); ylabel('Estimation Error');
legend('time varying filter', 'steady state filter');

figure;
plot(t,z0,'r-', t,zk,'b:');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time'); ylabel('Estimation Error');
legend('\theta = 7/16', '\theta = 0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, ztilde] = HinfFilter(SteadyState, theta)

% Initialize the random number generator seed.
% Try 0 and 1 for qualitatively different results.
randn('state', 0); 
% Define system matrices.
A = 1;
C = 1;
L = 1;
Q = 1;
R = 1;
S = 1;
P = 1;
tf = 5; % simulation time
dt = 0.1; % simulation step size
x = 0; % initial state
xhat = 0; % initial state estimate
enorm = 0;
ztildenorm = 0;
% Create arrays for plotting data
zArray = []; 
zhatArray = [];
for t = 0 : dt : tf
    % Generate noise
    w = 10 * randn(1);
    v = 10 * randn(1);
    v = 10 + v;
    % Simulate system dynamics
    xdot = A * x + w;
    x = x + xdot * dt;
    y = C * x + v;
    z = L * x;
    % Compute filter gain
    if (SteadyState == 1)
        K = 4;
    else
        K = P * C' * inv(R);
        Pdot = A * P + P * A' + Q - K * C * P + theta * P * L' * S * L * P;
        P = P + Pdot * dt;
    end
    % Filter dynamics
    xhatdot = A * xhat + K * (y - C * xhat);
    xhat = xhat + xhatdot * dt;
    zhat = L * xhat;
    % Save date for later plotting
    zArray = [zArray z];
    zhatArray = [zhatArray zhat];
    % Compute integral of noise squared, and integral of estimation error
    % squared.
    enorm = enorm + (w^2 + v^2) * dt;
    ztildenorm = ztildenorm + (z - zhat)^2 * dt;
end
% Plot data
t = 0 : dt : tf;
figure;
plot(t,zArray,'r', t,zhatArray,'b:');
if (SteadyState == 1)
    TitleStr = 'H\infty filter simulation results: steady state filter';
else
    TitleStr = 'H\infty filter simulation results: time varying filter';
end
title(TitleStr, 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time'); ylabel('');
legend('true output', 'estimated output');

figure;
ztilde = zArray - zhatArray;
plot(t,ztilde);
if (SteadyState == 1)
    TitleStr = 'H\infty estimation error: steady state filter (\theta = ';
else
    TitleStr = 'H\infty estimation error: time varying filter (\theta = ';
end
TitleStr = [TitleStr, num2str(theta), ')'];
title(TitleStr, 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time'); ylabel('');

disp(['infinity norm = ', num2str(sqrt(ztildenorm/enorm))]);