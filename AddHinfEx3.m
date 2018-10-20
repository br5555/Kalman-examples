function AddHinfEx3

% Robust Kalman filtering with model uncertainties.
% Based on paper by Hung and Fang.

J = 10; % motor moment of inertia

JFilter = 100 * J; % assumed moment of inertia
%JFilter = J;

Friction = 40; % viscous friction coefficient
c = 5; % voltage-to-torque constant
alpha = Friction / J;
alphaFilter = Friction / JFilter;
k = c / J;
u = 0; % input voltage
sigmaw = 2; % std dev of torque disturbance
sigmav = 0.2; % std dev of angle measurement noise (degrees)
sigmav = sigmav * pi / 180;

tf = 10; % simulation time (seconds)
%tf = 40;

A = [0 1; 0 -alpha]; % continuous time system matrix
AFilter = [0 1; 0 -alphaFilter]; % assumed system matrix
B = [0; c/J]; % continuous time input matrix
Bw = [0; 1/JFilter]; % continuous time noise input matrix
Q = Bw * sigmaw * Bw'; % process noise covariance
H = [1 0]; % measurement matrix

x0 = [0; 10]; % initial state

% Continuous time simulation
dt = 0.01;
x = x0;
yArr = [];
for t = 0 : dt : tf
    xdot = A * x + B * u;
    x = x + xdot * dt;
    y = H * x;
    yArr = [yArr y];
end

% Discrete time simulation
T = 0.1;
Ad = expm(A * T); % discretized system matrix
AdFilter = expm(AFilter * T); % assumed discretized system matrix
Bd = c / Friction * [T - 1/alpha + exp(-alpha*T)/alpha ; 1 - exp(-alpha*T)]; % discretized input matrix
Bwd = 1 / Friction * [T - 1/alpha + exp(-alpha*T)/alpha ; 1 - exp(-alpha*T)]; % discretized noise input matrix
x = x0;
xhat = x;
Pplus = Q;
ydArr = [];
xtildeArr = [];
% Robust filter parameters
N = .00001*eye(2);
M1 = eye(2);
M2 = [0 0];
D1 = [Bwd'; 0 0]';
D2 = [0 1];
theta = .6;
alpha = .10;
eps = 1e-8;
S1 = 1e4*eye(2);
S2 = 1e4*eye(2);
Q1 = S1;
Q2 = S2;
R11 = D1 * D1' + alpha * M1 * M1';
R12 = D1 * D2' + alpha * M1 * M2';
R22 = D2 * D2' + alpha * M2 * M2';
L = .001*eye(2);
xhatr = x;
xrobustArr = [];
for t = 0 : T : tf;
    % System and measurement simulation
    y = H * x + sigmav * randn;
    x = Ad * x + Bd * u + Bwd * sigmaw * randn;
    % Kalman filter
    Pminus = AdFilter * Pplus * AdFilter' + Q;
    K = Pminus * H' * inv(H * Pminus * H' + sigmav^2);
    xhat = AdFilter * xhat + Bd * u;
    xtildeArr = [xtildeArr x-xhat];
    xhat = xhat + K * (y - H * xhat);
    Pplus = (eye(2) - K * H) * Pminus * (eye(2) - K * H)' + K * sigmav^2 * K';
    % Save data for plotting
    ydArr = [ydArr y];
    % Robust filter equations
    A = AdFilter;
    R1 = inv(inv(Q2) - N' * N / alpha) * A';
    R2 = inv(R1) * inv(inv(Q2) - N' * N / alpha) * inv(R1');
    A1 = A + R11 * inv(R1);
    C1 = H + R12' * inv(R1);
    temp1 = inv(inv(Q1) - theta^2 * L' * L);
    R = C1 * temp1 * C1' + R12' * R2 * R12 + R22;
    % Robust state estimate 
    G = (A1 * temp1 * C1' + R11 * R2 * R12 + R12) * inv(R);
    F = A1 - G * C1;
    xhatr = F * xhatr + G * y;
    xrobustArr = [xrobustArr x-xhatr];
    % Riccati equations for robust filter
    Q1 = A1 * temp1 * A1' + R11 + R11 * R2 * R11' - ...
        (A1 * temp1 * C1' + R11 * R2 * R12 + R12) * inv(R) * (A1 * temp1 * C1' + R11 * R2 * R12 + R12)' + eps * eye(2);
    Q2 = A * Q2 * A' + A * Q2 * N' * inv(alpha * eye(2) - N * Q2 * N') * N * Q2 * A' + R11 + eps * eye(2);
    % Check conditions for validity of robust state estimate
    lambda = eig(Q1);
    for i = 1 : 2
        if ~isreal(lambda(i)) || (lambda(i) <= 0)
            disp(['Q1 is not positive definite - t = ', num2str(t)]);
            return;
        end
    end
    lambda = eig(eye(2) / theta^2 - L * Q1 * L');
    for i = 1 : 2
        if ~isreal(lambda(i)) || (lambda(i) <= 0)
            disp(['Q1 condition not satisfied - t = ', num2str(t)]);
            return;
        end
    end
    lambda = eig(Q2);
    for i = 1 : 2
        if ~isreal(lambda(i)) || (lambda(i) <= 0)
            disp(['Q2 is not positive definite - t = ', num2str(t)]);
            return;
        end
    end
    lambda = eig(alpha * eye(2) - N * Q2 * N');
    for i = 1 : 2
        if ~isreal(lambda(i)) || (lambda(i) <= 0)
            disp(['Q2 condition not satisfied - t = ', num2str(t)]);
            return;
        end
    end
end

close all;
t = 0 : dt : tf;
td = 0 : T : tf;

%figure;
%plot(t, yArr, td, ydArr);
%legend('Continuous time', 'Discrete time');

xtildeArr = xtildeArr * 180 / pi;
xrobustArr = xrobustArr * 180 / pi;

figure;
plot(td, xtildeArr(1,:), 'r:', td, xrobustArr(1,:), 'b-');
set(gca,'FontSize',12); set(gcf,'Color','White');
legend('Kalman filter', 'Robust filter');
xlabel('seconds'); ylabel('deg');

figure;
plot(td, xtildeArr(2,:), 'r:', td, xrobustArr(2,:), 'b-');
set(gca,'FontSize',12); set(gcf,'Color','White');
legend('Kalman filter', 'Robust filter');
xlabel('seconds'); ylabel('deg/sec');

% Compute RMS estimation errors
iStart = -1;
for i = 1 : size(xtildeArr,2)
    if (abs(xtildeArr(1,i)) <= 1) && (iStart < 0)
        iStart = i;
    elseif abs(xtildeArr(1,i)) > 1 
        iStart = -1;
    end
end
len = size(xtildeArr,2) - iStart + 1;
KalmanRMSPos = sqrt(norm(xtildeArr(1,iStart:end))^2 / len);
KalmanRMSVel = sqrt(norm(xtildeArr(2,iStart:end))^2 / len);

iStart = -1;
for i = 1 : size(xrobustArr,2)
    if (abs(xrobustArr(1,i)) <= 1) && (iStart < 0)
        iStart = i;
    elseif abs(xrobustArr(1,i)) > 1
        iStart = -1;
    end
end
len = size(xrobustArr,2) - iStart + 1;
RobustRMSPos = sqrt(norm(xrobustArr(1,iStart:end))^2 / len);
RobustRMSVel = sqrt(norm(xrobustArr(2,iStart:end))^2 / len);

disp(['Kalman filter RMS estimation errors = ', num2str(KalmanRMSPos), ', ', num2str(KalmanRMSVel)]);
disp(['Robust filter RMS estimation errors = ', num2str(RobustRMSPos), ', ', num2str(RobustRMSVel)]);
disp(['Q1 = diag(', num2str(Q1(1,1)), ', ', num2str(Q1(2,2)), ')']);