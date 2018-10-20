function [xsmooth] = FixPtSmooth(duration, dt, measnoise, PlotFlag)

% function FixPtSmooth(duration, dt)
%
% Fixed point smoother simulation for a two state navigation problem.
% INPUTS
%   duration = length of simulation (seconds)
%   dt = step size (seconds)
%   measnoise = standard deviation of position measurement noise
%   PlotFlag = true/false flag saying whether or not to generate plots
% OUTPUT
%   xsmooth = estimate of the initial state after all measurements have been processed

if ~exist('duration', 'var')
    duration = 10;
end
if ~exist('dt', 'var')
    dt = 0.1;
end
if ~exist('measnoise', 'var')
    measnoise = 1; 
end
if ~exist('PlotFlag', 'var')
    PlotFlag = true;
end

accelnoise = 0.2; % acceleration noise (feet/sec^2)

a = [1 dt; 0 1]; % transition matrix
b = [dt^2/2; dt]; % input matrix
c = [1 0]; % measurement matrix

Sw = accelnoise^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2]; % process noise covariance
P = [1 0; 0 1]; % initial the a posteriori estimation covariance
Sz = measnoise^2; % measurement error covariance

x = [0; 0]; % initial state vector
xhat = x + sqrt(P) * 2 * ones(size(x)); % initial the a posteriori state estimate (in error by 2 sigma)

u = 10; % Use a constant commanded acceleration

xhat = a * xhat + b * u; % initial the a priori state estimate
P = a * P * a' + Sw; % initial the a priori covariance

% Initialize the fixed point smoother
Sigma = P;
Pi = P;
xsmooth = xhat;

% Initialize arrays for later plotting.
pos = []; % true position array
poshat = []; % estimated position array
possmooth = []; % smoothed position array
posmeas = []; % measured position array
vel = []; % true velocity array
velhat = []; % estimated velocity array
velsmooth = []; % smoothed velocity array
TrPArray = []; % Trace of standard estimation covariance
TrPiArray = []; % Trace of smoothed estimation covariance

randn('state', sum(100*clock)); % initialize the random number generator

% Simulate the linear system and the noisy measurement.
ProcessNoise = accelnoise * [(dt^2/2)*randn; dt*randn];
x = a * x + b * u + ProcessNoise;
y = c * x + measnoise * randn;
% Form the Innovation vector.
Inn = y - c * xhat;

for t = 0 : dt: duration
    % Save some parameters for plotting later.
    pos = [pos; x(1)];
    posmeas = [posmeas; y];
    poshat = [poshat; xhat(1)];
    possmooth = [possmooth; xsmooth(1)];
    vel = [vel; x(2)];
    velhat = [velhat; xhat(2)];
    velsmooth = [velsmooth; xsmooth(2)];
    TrPArray = [TrPArray trace(P)];
    TrPiArray = [TrPiArray trace(Pi)];
    
    % Compute the covariance of the Innovation.
    s = c * P * c' + Sz;
    % Form the Kalman Gain matrix.
    K = a * P * c' * inv(s);
    
    % Compute the smoothed estimate.
    Lambda = Sigma * c' * inv(s); % smoother gain
    Pi = Pi - Sigma * c' * Lambda'; % covariance of smoothed estimate
    Sigma = Sigma * (a - K * c)'; % cross covariance of standard and smoothed estimates
    xsmooth = xsmooth + Lambda * Inn; % smoothed estimate of the initial state
    
    % Update the state estimate.
    xhat = a * xhat + K * Inn + b * u;
    % Compute the covariance of the estimation error.
    P = a * P * a' - a * P * c' * K' + Sw;
    
    % Simulate the linear system and the noisy measurement.
    ProcessNoise = accelnoise * [(dt^2/2)*randn; dt*randn];
    x = a * x + b * u + ProcessNoise;
    y = c * x + measnoise * randn;
    % Form the Innovation vector.
    Inn = y - c * xhat;
    
end

if ~PlotFlag
    return;
end

disp(['Pi(1,1) = ', num2str(Pi(1,1)), ', Pi(2,2) = ', num2str(Pi(2,2))]);
Improve = 100 * (TrPArray(1) - trace(Pi)) / TrPArray(1);
disp(['Smoothing Improvement = ',num2str(Improve),' %']);

% Plot the results
close all;
t = 0 : dt : duration;

figure;
plot(t,TrPArray,'r-', t,TrPiArray,'b:');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time');
ylabel('Analytical Estimation Covariance', 'FontSize', 12);
legend('Standard Filter', 'Smoothed Filter'); 

figure;
plot(t,possmooth-pos(1),'r', t,velsmooth-vel(1),'b:');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time');
ylabel('Smoothing Error at Initial Time', 'FontSize', 12);
legend('Position Error', 'Velocity Error'); 

figure;
plot(t,pos,'r-', t,posmeas,'b:', t,poshat,'k--');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time');
ylabel('Position');
legend('True Position', 'Measured Position', 'Estimated Position');

figure;
plot(t,pos-posmeas,'r-', t,pos-poshat,'b:');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time (sec)');
ylabel('Position Error');
legend('Position Measurement Error', 'Position Estimation Error');

figure;
plot(t,vel,'r-', t,velhat,'b:');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time');
ylabel('Velocity');
legend('True Velocity', 'Estimated Velocity');

figure;
plot(t,vel-velhat);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time');
ylabel('Velocity Estimation Error');
