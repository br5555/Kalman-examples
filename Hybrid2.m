function [EKFErr, EKF2Err, EKFiErr] = Hybrid2

% Optimal State Estimation, by Dan Simon
%
% Hybrid second order Kalman filter example.
% Track a body falling through the atmosphere.
% This example is taken from [Ath68].
% Outputs:
%   EKFErr = 1st order EKF estimation error (RMS) (altitude, velocity, ballistic coefficient)
%   EKF2Err = 2nd order EKF estimation error (RMS) (altitude, velocity, ballistic coefficient)
%   EKFiErr = iterated EKF estimation error (RMS) (altitude, velocity, ballistic coefficient)

% Revision March 3, 2009 - The Q terms were being incorrectly multiplied by dt in the RungeKutta routine,
% which was incorrect (although it does not affect the results in this case since Q = 0).

global rho0 g k dt

rho0 = 2; % lb-sec^2/ft^4
g = 32.2; % ft/sec^2
k = 2e4; % ft
R = 10^4; % measurement noise variance (ft^2)
Q = diag([0 0 0]); % process noise covariance
M = 10^5; % horizontal range of position sensor
a = 10^5; % altitude of position sensor

x = [3e5; -2e4; 1e-3]; % true state
xhat = [3e5; -2e4; 1e-3]; % first order EKF estimate
xhat2 = xhat; % second order EKF estimate
xhati = xhat; % iterated EKF estimate
N = 7; % max number of iterations to execute in the iterated EKF

P = diag([1e6 4e6 10]);
P2 = P;
Pi = P;

T = 1; % measurement time step
randn('state',sum(100*clock)); % random number generator seed

tf = 30; % simulation length (seconds)
dt = 0.02; % time step for integration (seconds)
xArray = x;
xhatArray = xhat;
xhat2Array = xhat2;
xhatiArray = xhati;
Parray = diag(P);
P2array = diag(P2);
Piarray = diag(Pi);

for t = T : T : tf
    % Simulate the system.
    for tau = dt : dt : T
        % Fourth order Runge Kutta ingegration
        [dx1, dx2, dx3, dx4] = RungeKutta(x);
        x = x + (dx1 + 2 * dx2 + 2 * dx3 + dx4) / 6;
        x = x + sqrt(dt * Q) * [randn; randn; randn] * dt;
    end
    % Simulate the noisy measurement.
    z = sqrt(M^2 + (x(1)-a)^2) + sqrt(R) * randn;
    
    % First order Kalman filter.
    % Simulate the continuous-time part of the first order Kalman filter (time update).
    for tau = dt : dt : T
        [dx1, dx2, dx3, dx4] = RungeKutta(xhat);
        xhat = xhat + (dx1 + 2 * dx2 + 2 * dx3 + dx4) / 6;
        F = [0 1 0; -rho0 * exp(-xhat(1)/k) * xhat(2)^2 / 2 / k * xhat(3) ...
                rho0 * exp(-xhat(1)/k) * xhat(2) * xhat(3) ...
                rho0 * exp(-xhat(1)/k) * xhat(2)^2 / 2; ...
                0 0 0];
        H = [(xhat(1) - a) / sqrt(M^2 + (xhat(1)-a)^2) 0 0];
        [dP1, dP2, dP3, dP4] = RungeKuttaP(P, F, Q, dt, H, R);
        P = P + (dP1 + 2 * dP2 + 2 * dP3 + dP4) / 6;
    end
    % Force the ballistic coefficient estimate to be non-negative.
    xhat(3) = max(xhat(3), 0);
    % Simulate the discrete-time part of the first order Kalman filter (measurement update).
    H = [(xhat(1) - a) / sqrt(M^2 + (xhat(1)-a)^2) 0 0];
    K = P * H' * inv(H * P * H' + R);
    zhat = sqrt(M^2 + (xhat(1)-a)^2);
    xhat = xhat + K * (z - zhat);
    % Force the ballistic coefficient estimate to be non-negative.
    xhat(3) = max(xhat(3), 0);
    P = (eye(3) - K * H) * P;
    
    % Iterated Kalman filter.
    % Simulate the continuous-time part of the iterated Kalman filter (time update).
    for tau = dt : dt : T
        [dx1, dx2, dx3, dx4] = RungeKutta(xhati);
        xhati = xhati + (dx1 + 2 * dx2 + 2 * dx3 + dx4) / 6;
        F = [0 1 0; -rho0 * exp(-xhati(1)/k) * xhati(2)^2 / 2 / k * xhati(3) ...
                rho0 * exp(-xhati(1)/k) * xhati(2) * xhati(3) ...
                rho0 * exp(-xhati(1)/k) * xhati(2)^2 / 2; ...
                0 0 0];
        H = [(xhati(1) - a) / sqrt(M^2 + (xhati(1)-a)^2) 0 0];
        [dP1, dP2, dP3, dP4] = RungeKuttaP(Pi, F, Q, dt, H, R);
        Pi = Pi + (dP1 + 2 * dP2 + 2 * dP3 + dP4) / 6;
    end
    % Force the ballistic coefficient estimate to be non-negative.
    xhati(3) = max(xhati(3), 0);
    % Simulate the discrete time part of the iterated Kalman filter (measurement update);
    xhatminus = xhati;
    Pminus = Pi;
    for i = 1 : N
        H = [(xhati(1) - a) / sqrt(M^2 + (xhati(1)-a)^2) 0 0];
        K = Pminus * H' * inv(H * Pminus * H' + R);
        zhat = sqrt(M^2 + (xhati(1)-a)^2);
        xhati = xhatminus + K * ((z - zhat) - H * (xhatminus - xhati));
        Pi = (eye(3) - K * H) * Pminus;
        % Force the ballistic coefficient estimate to be non-negative.
        xhati(3) = max(xhati(3), 0);
    end
    
    % Second order Kalman filter.
    % Simulate the continuous-time part of the second order Kalman filter (time update).
    for tau = dt : dt : T
        [dx, dP] = RungeKutta2(xhat2, P2, Q, R);
        xhat2 = xhat2 + (dx(:,1) + 2 * dx(:,2) + 2 * dx(:,3) + dx(:,4)) / 6;
        P2 = P2 + (dP(:,:,1) + 2 * dP(:,:,2) + 2 * dP(:,:,3) + dP(:,:,4)) / 6;
    end
    % Force the ballistic coefficient estimate to be non-negative.
    xhat2(3) = max(xhat2(3), 0);
    % Simulate the discrete-time part of the second order Kalman filter (measurement update).
    H = [(xhat2(1) - a) / sqrt(M^2 + (xhat2(1)-a)^2) 0 0];
    zhat = sqrt(M^2 + (xhat2(1)-a)^2);
    D = zeros(3,3);
    D(1,1) = 1/zhat * (1 - (xhat2(1) - a)^2 / zhat / zhat);
    L = 1/2 * trace(D * P2 * D * P2);
    K = P2 * H' * inv(H * P2 * H' + R + L);
    pie = 1/2 * K * [1] * trace(D * P2);
    xhat2 = xhat2 + K * (z - zhat) - pie;
    % Force the ballistic coefficient estimate to be non-negative.
    xhat2(3) = max(xhat2(3), 0);
    P2 = P2 - P2 * H' * (H * P2 * H' + R + L)^(-1) * H * P2;
    
    % Save data for plotting.
    xArray = [xArray x];
    xhatArray = [xhatArray xhat];
    xhat2Array = [xhat2Array xhat2];
    xhatiArray = [xhatiArray xhati];
    Parray = [Parray diag(P)];
    P2array = [P2array diag(P2)];
    Piarray = [Piarray diag(Pi)];
end

close all;
t = 0 : T : tf;
figure; 
semilogy(t, abs(xArray(1,:) - xhatArray(1,:)), 'b', 'LineWidth', 2); hold;
semilogy(t, abs(xArray(1,:) - xhat2Array(1,:)), 'r:', 'LineWidth', 2);
semilogy(t, abs(xArray(1,:) - xhatiArray(1,:)), 'k--', 'LineWidth', 2);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('Altitude Estimation Error');
legend('Kalman filter', 'Second order filter', 'Iterated Kalman filter');

AltErr = std(xArray(1,:) - xhatArray(1,:));
disp(['1st Order EKF RMS altitude estimation error = ', num2str(AltErr)]);
Alt2Err = std(xArray(1,:) - xhat2Array(1,:));
disp(['2nd Order EKF RMS altitude estimation error = ', num2str(Alt2Err)]);
AltiErr = std(xArray(1,:) - xhatiArray(1,:));
disp(['Iterated EKF RMS altitude estimation error = ', num2str(AltiErr)]);

figure; 
semilogy(t, abs(xArray(2,:) - xhatArray(2,:)), 'b', 'LineWidth', 2); hold;
semilogy(t, abs(xArray(2,:) - xhat2Array(2,:)), 'r:', 'LineWidth', 2);
semilogy(t, abs(xArray(2,:) - xhatiArray(2,:)), 'k--', 'LineWidth', 2);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('Velocity Estimation Error');
legend('Kalman filter', 'Second order filter', 'Iterated Kalman filter');

VelErr = std(xArray(2,:) - xhatArray(2,:));
disp(['1st Order EKF RMS velocity estimation error = ', num2str(VelErr)]);
Vel2Err = std(xArray(2,:) - xhat2Array(2,:));
disp(['2nd Order EKF RMS velocity estimation error = ', num2str(Vel2Err)]);
VeliErr = std(xArray(2,:) - xhatiArray(2,:));
disp(['Iterated EKF RMS velocity estimation error = ', num2str(VeliErr)]);

figure; 
semilogy(t, abs(xArray(3,:) - xhatArray(3,:)), 'b', 'LineWidth', 2); hold;
semilogy(t, abs(xArray(3,:) - xhat2Array(3,:)), 'r:', 'LineWidth', 2);
semilogy(t, abs(xArray(3,:) - xhatiArray(3,:)), 'k--', 'LineWidth', 2);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('Ballistic Coefficient Estimation Error');
legend('Kalman filter', 'Second order filter', 'Iterated Kalman filter');

BallErr = std(xArray(3,:) - xhatArray(3,:));
disp(['1st Order EKF RMS ballistic coefficient estimation error = ', num2str(BallErr)]);
Ball2Err = std(xArray(3,:) - xhat2Array(3,:));
disp(['2nd Order EKF RMS ballistic coefficient estimation error = ', num2str(Ball2Err)]);
BalliErr = std(xArray(3,:) - xhatiArray(3,:));
disp(['Iterated EKF RMS ballistic coefficient estimation error = ', num2str(BalliErr)]);

figure;
plot(t, xArray(1,:), 'LineWidth', 2);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('True Altitude');
grid;

figure;
plot(t, xArray(2,:), 'LineWidth', 2);
title('Falling Body Simulation', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('True Velocity');
grid;

EKFErr = [AltErr; VelErr; BallErr];
EKF2Err = [Alt2Err; Vel2Err; Ball2Err];
EKFiErr = [AltiErr; VeliErr; BalliErr];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx1, dx2, dx3, dx4] = RungeKutta(x)
% Fourth order Runge Kutta integration for the falling body system.
global rho0 g k dt
dx1(1,1) = x(2);
dx1(2,1) = rho0 * exp(-x(1)/k) * x(2)^2 / 2 * x(3) - g;
dx1(3,1) = 0;
dx1 = dx1 * dt;
xtemp = x + dx1 / 2;
dx2(1,1) = xtemp(2);
dx2(2,1) = rho0 * exp(-xtemp(1)/k) * xtemp(2)^2 / 2 * xtemp(3) - g;
dx2(3,1) = 0;
dx2 = dx2 * dt;
xtemp = x + dx2 / 2;
dx3(1,1) = xtemp(2);
dx3(2,1) = rho0 * exp(-xtemp(1)/k) * xtemp(2)^2 / 2 * xtemp(3) - g;
dx3(3,1) = 0;
dx3 = dx3 * dt;
xtemp = x + dx3;
dx4(1,1) = xtemp(2);
dx4(2,1) = rho0 * exp(-xtemp(1)/k) * xtemp(2)^2 / 2 * xtemp(3) - g;
dx4(3,1) = 0;
dx4 = dx4 * dt;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx, dP] = RungeKutta2(x, P, Q, R)
% Fourth order Runge Kutta integration for the falling body system for the second order Kalman filter estimate.
global rho0 g k dt
% First Runge Kutta increment for x and P.
xtemp = x;
Ptemp = P;
exp1 = exp(-xtemp(1) / k);
F = [0 1 0; -rho0 * exp1 * xtemp(2)^2 / 2 / k * xtemp(3) ...
        rho0 * exp1 * xtemp(2) * xtemp(3) ...
    rho0 * exp1 * xtemp(2)^2 / 2; ...
        0 0 0];
F1 = zeros(3,3);
F2 = [rho0 / k / k * exp1 * xtemp(2)^2 / 2 * xtemp(3)  -rho0 / k * exp1 * xtemp(2) * xtemp(3)  -rho0 / k * exp1 * xtemp(2)^2 / 2 ;
    -rho0 / k * exp1 * xtemp(2) * x(3)  rho0 * exp1 * xtemp(3)  rho0 * exp(-xtemp(1) / k) * xtemp(2) ;
    -rho0 / k * exp1 * xtemp(2)^2 / 2  rho0 * exp1 * xtemp(2)  0 ];
F3 = zeros(3,3);
ubreve = 1/2 * ([1 ; 0 ; 0] * trace(F1 * Ptemp) + [0 ; 1 ; 0] * trace(F2 * Ptemp) + [0 ; 0 ; 1] * trace(F3 * Ptemp));
dx(1,1) = xtemp(2);
dx(2,1) = rho0 * exp(-xtemp(1)/k) * xtemp(2)^2 / 2 * xtemp(3) - g;
dx(3,1) = 0;
dx(:,1) = (dx(:,1) + ubreve) * dt;
xtemp = x + dx(:,1) / 2;
dP(:,:,1) = F * Ptemp + Ptemp * F' + Q * dt;
dP(:,:,1) = dP(:,:,1) * dt;
Ptemp = P + dP(:,:,1) / 2;
% Second Runge Kutta increment for x and P.
exp1 = exp(-xtemp(1) / k);
F = [0 1 0; -rho0 * exp1 * xtemp(2)^2 / 2 / k * xtemp(3) ...
        rho0 * exp1 * xtemp(2) * xtemp(3) ...
    rho0 * exp1 * xtemp(2)^2 / 2; ...
        0 0 0];
F1 = zeros(3,3);
F2 = [rho0 / k / k * exp1 * xtemp(2)^2 / 2 * xtemp(3)  -rho0 / k * exp1 * xtemp(2) * xtemp(3)  -rho0 / k * exp1 * xtemp(2)^2 / 2 ;
    -rho0 / k * exp1 * xtemp(2) * x(3)  rho0 * exp1 * xtemp(3)  rho0 * exp(-xtemp(1) / k) * xtemp(2) ;
    -rho0 / k * exp1 * xtemp(2)^2 / 2  rho0 * exp1 * xtemp(2)  0 ];
F3 = zeros(3,3);
ubreve = 1/2 * ([1 ; 0 ; 0] * trace(F1 * Ptemp) + [0 ; 1 ; 0] * trace(F2 * Ptemp) + [0 ; 0 ; 1] * trace(F3 * Ptemp));
dx(1,2) = xtemp(2);
dx(2,2) = rho0 * exp(-xtemp(1)/k) * xtemp(2)^2 / 2 * xtemp(3) - g;
dx(3,2) = 0;
dx(:,2) = (dx(:,2) + ubreve) * dt;
xtemp = x + dx(:,2) / 2;
dP(:,:,2) = F * Ptemp + Ptemp * F' + Q * dt;
dP(:,:,2) = dP(:,:,2) * dt;
Ptemp = P + dP(:,:,2) / 2;
% Third Runge Kutta increment for x and P.
exp1 = exp(-xtemp(1) / k);
F = [0 1 0; -rho0 * exp1 * xtemp(2)^2 / 2 / k * xtemp(3) ...
        rho0 * exp1 * xtemp(2) * xtemp(3) ...
    rho0 * exp1 * xtemp(2)^2 / 2; ...
        0 0 0];
F1 = zeros(3,3);
F2 = [rho0 / k / k * exp1 * xtemp(2)^2 / 2 * xtemp(3)  -rho0 / k * exp1 * xtemp(2) * xtemp(3)  -rho0 / k * exp1 * xtemp(2)^2 / 2 ;
    -rho0 / k * exp1 * xtemp(2) * x(3)  rho0 * exp1 * xtemp(3)  rho0 * exp(-xtemp(1) / k) * xtemp(2) ;
    -rho0 / k * exp1 * xtemp(2)^2 / 2  rho0 * exp1 * xtemp(2)  0 ];
F3 = zeros(3,3);
ubreve = 1/2 * ([1 ; 0 ; 0] * trace(F1 * Ptemp) + [0 ; 1 ; 0] * trace(F2 * Ptemp) + [0 ; 0 ; 1] * trace(F3 * Ptemp));
dx(1,3) = xtemp(2);
dx(2,3) = rho0 * exp(-xtemp(1)/k) * xtemp(2)^2 / 2 * xtemp(3) - g;
dx(3,3) = 0;
dx(:,3) = (dx(:,3) + ubreve) * dt;
xtemp = x + dx(:,3);
dP(:,:,3) = F * Ptemp + Ptemp * F' + Q * dt;
dP(:,:,3) = dP(:,:,3) * dt;
Ptemp = P + dP(:,:,3);
% Fourth Runge Kutta increment for x and P.
exp1 = exp(-xtemp(1) / k);
F = [0 1 0; -rho0 * exp1 * xtemp(2)^2 / 2 / k * xtemp(3) ...
        rho0 * exp1 * xtemp(2) * xtemp(3) ...
    rho0 * exp1 * xtemp(2)^2 / 2; ...
        0 0 0];
F1 = zeros(3,3);
F2 = [rho0 / k / k * exp1 * xtemp(2)^2 / 2 * xtemp(3)  -rho0 / k * exp1 * xtemp(2) * xtemp(3)  -rho0 / k * exp1 * xtemp(2)^2 / 2 ;
    -rho0 / k * exp1 * xtemp(2) * x(3)  rho0 * exp1 * xtemp(3)  rho0 * exp(-xtemp(1) / k) * xtemp(2) ;
    -rho0 / k * exp1 * xtemp(2)^2 / 2  rho0 * exp1 * xtemp(2)  0 ];
F3 = zeros(3,3);
ubreve = 1/2 * ([1 ; 0 ; 0] * trace(F1 * Ptemp) + [0 ; 1 ; 0] * trace(F2 * Ptemp) + [0 ; 0 ; 1] * trace(F3 * Ptemp));
dx(1,4) = xtemp(2);
dx(2,4) = rho0 * exp(-xtemp(1)/k) * xtemp(2)^2 / 2 * xtemp(3) - g;
dx(3,4) = 0;
dx(:,4) = (dx(:,4) + ubreve) * dt;
dP(:,:,4) = F * Ptemp + Ptemp * F' + Q * dt;
dP(:,:,4) = dP(:,:,4) * dt;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dP1, dP2, dP3, dP4] = RungeKuttaP(P, F, Q, dt, H, R)
% Fourth order Runge Kutta integration for the covariance of the Kalman
% filter for the falling body system.
dP1 = F * P + P * F' + Q - P * H' * inv(R) * H * P;
dP1 = dP1 * dt;
Ptemp = P + dP1 / 2;
dP2 = F * Ptemp + Ptemp * F' + Q - Ptemp * H' * inv(R) * H * Ptemp;
dP2 = dP2 * dt;
Ptemp = P + dP2 / 2;
dP3 = F * Ptemp + Ptemp * F' + Q - Ptemp * H' * inv(R) * H * Ptemp;
dP3 = dP3 * dt;
Ptemp = P + dP3;
dP4 = F * Ptemp + Ptemp * F' + Q - Ptemp * H' * inv(R) * H * Ptemp;
dP4 = dP4 * dt;
return