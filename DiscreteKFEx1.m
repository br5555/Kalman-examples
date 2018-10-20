function DiscreteKFEx1(N)

% Discrete time Kalman filter for position estimation of a Newtonian system.
% This example illustrates the effectiveness of the Kalman filter for state
% estimation. It also shows how the variance of the estimation error 
% propagates between time steps and decreases as each measurement is processed.
% INPUTS: N = number of time steps.

if ~exist('N', 'var')
    N = 6;
end

T = 5; % time between measurements
sigma = 30; % position measurement standard deviation
R = sigma^2;
P0 = [100 0 0; 0 10 0; 0 0 1]; % initial state estimate uncertainty
A = [0 1 0; 0 0 1; 0 0 0];
H = [1 0 0];
F = [1 T T*T/2; 0 1 T; 0 0 1]; % state transition matrix
x = [1; 1; 1]; % initial state
xhat = x; % initial state estimate

posArray = [];
xhatArray = [];
yArray = [];
Pplus = P0;
Varminus = [];
Varplus = [P0(1,1)];
KArray = [];

for k = 1 : N
    % Simulate the system and measurement
    x = F * x;
    y = H * x + sigma * randn;
    % Estimate the state
    Pminus = F * Pplus * F';
    K = Pminus * H' * inv(H * Pminus * H' + R);
    xhat = F * xhat;
    xhat = xhat + K * (y - H * xhat);
    Pplus = (eye(3) - K * H) * Pminus * (eye(3) - K * H)' + K * R * K';
    % Save data for plotting
    posArray = [posArray x(1)];
    xhatArray = [xhatArray xhat];
    yArray = [yArray y];
    Varminus = [Varminus Pminus(1,1)];
    Varplus = [Varplus Pplus(1,1)];
    KArray = [KArray K];
end

% Plot the results
close all;
k = 1 : N;
plot(k, yArray-posArray, 'r:');
hold;
plot(k, xhatArray(1,:)-posArray, 'b-');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('time step'); ylabel('position');
legend('measurement error', 'estimation error');

figure; hold;
for k = 1 : N-1
    plot([k-1 k], [Varplus(k) Varminus(k+1)]);
    plot([k k], [Varminus(k+1) Varplus(k+1)]);
end
set(gca,'FontSize',12); set(gcf,'Color','White'); set(gca,'Box','on');
xlabel('time step');
ylabel('position estimation error variance');
