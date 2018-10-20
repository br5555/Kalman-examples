function [PctImprovement] = FixLagSmooth(duration, dt, N, measnoise)

% function FixLagSmooth(duration, dt)
% Updated September 2016
%
% Fixed lag smoother simulation for a vehicle traveling along a road.
% INPUTS
%   duration = length of simulation (seconds)
%   dt = step size (seconds)
%   N = one less than smoothing lag
%   measnoise
% OUTPUTS
%   PctImprovement = array containing the percent improvement due to
%                    smoothing with a lag of 1, 2, ..., N+1

accelnoise = 10; % acceleration noise (feet/sec^2)

if ~exist('duration', 'var')
    duration = 10;
end
if ~exist('dt', 'var')
    dt = 0.1;
end
if ~exist('N', 'var')
    N = 30;
end
if ~exist('measnoise', 'var')
    measnoise = 10; % position measurement noise (feet)
end

a = [1 dt; 0 1]; % transition matrix
b = [dt^2/2; dt]; % input matrix
c = [1 0]; % measurement matrix
x = [0; 0]; % initial state vector
xhat = x; % initial state estimate

Sw = accelnoise^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2]; % process noise covariance
P = [1 0; 0 1]; % initial estimation covariance
Sz = measnoise^2; % measurement error covariance

PCol = zeros(2,2,N+1);
PColOld = zeros(2,2,N+1);
PSmooth = zeros(2,2,N+1);
PSmoothOld = zeros(2,2,N+1);

% Initialize arrays for later plotting.
pos = []; % true position array
poshat = []; % estimated position array
posmeas = []; % measured position array
vel = []; % true velocity array
velhat = []; % estimated velocity array

kk=0;
for t = 0 : dt: duration,
    kk=kk+1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TIME INDEX
    % Use a constant commanded acceleration of 1 foot/sec^2.
    u = 1;
    % Simulate the linear system.
    ProcessNoise = accelnoise * [(dt^2/2)*randn; dt*randn];
    x = a * x + b * u + ProcessNoise;
    % Simulate the noisy measurement
    MeasNoise = measnoise * randn;
    y = c * x + MeasNoise;
    % Extrapolate the most recent state estimate to the present time.
    xhat = a * xhat + b * u;
    % Smoothed estimated at time (kk) given measurements up to and including time (kk-1) - this is the same as the prior xhat calculated on the previous line
    xhatSmoothArr(:,kk+1,kk) = xhat; 
    % Form the Innovation vector.
    Inn = y - c * xhat;
    % Compute the covariance of the Innovation.
    s = c * P * c' + Sz;
    % Form the Kalman Gain matrix.
    K = a * P * c' / s;
    KSmooth = K;
    % Update the state estimate.
    xhat = xhat + K * Inn;
    
    % Compute the covariance of the estimation error.
    PColOld(:,:,1) = P;
    PSmoothOld(:,:,1) = P;
    P = a * P * a' - a * P * c' / s * c * P * a' + Sw;
    % Save some parameters for plotting later.
    pos = [pos; x(1)];
    posmeas = [posmeas; y];
    poshat = [poshat; xhat(1)];
    vel = [vel; x(2)];
    velhat = [velhat; xhat(2)];
    for i = 1 : N+1
        KSmooth = PColOld(:,:,i) * c' / s;
        PSmooth(:,:,i+1) = PSmoothOld(:,:,i) - PColOld(:,:,i) * c' * KSmooth' * a';
        PCol(:,:,i+1) = PColOld(:,:,i) * (a - K * c)';
        if kk < i, break, end
        % Calculate smoothed estimated at time (kk+1-i) given measurements up to and including time kk
        xhatSmoothArr(:,kk+2-i,kk+1) = xhatSmoothArr(:,kk+2-i,kk) + KSmooth * Inn;
    end
    PSmoothOld = PSmooth;
    PColOld = PCol;
    
end

for i = 1 : N
    PctImprovement(i) = 100 * trace(P - PSmooth(:,:,i+1)) / trace(P);
end

close all
figure
plot(PctImprovement);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Number of lag intervals'); ylabel('Percent improvement');
title(['Measurement Noise = ',num2str(measnoise)], 'FontSize', 12);

figure
hold on
plot(0:kk-1, poshat,'b-')
plot(0:kk-1, posmeas,'r.')

Temp = NaN(kk, 1);
for i = 1 : kk-1
    Temp(i) = xhatSmoothArr(1, i, i+1);
end
plot(0:kk-1, Temp, 'k-')    % Smoothed first state with a lag of 1 time step

Temp = NaN(kk, 1);
for i = 1 : kk-5
    Temp(i) = xhatSmoothArr(1, i, i+5);
end
plot(0:kk-1, Temp, 'g')     % Smoothed first state with a lag of 5 time steps

Temp = NaN(kk, 1);
for i = 1 : kk-30
    Temp(i) = xhatSmoothArr(1, i, i+30);
end
plot(0:kk-1, Temp, 'm--')   % Smoothed first state with a lag of 30 time steps

set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time step')
legend('Posteriori estimation (time lag = 0)', 'Measurement', 'Smoothed estimate (time lag = 1)', 'Smoothed estimate (time lag = 5)', 'Smoothed estimate (time lag = 30)')
hold off
