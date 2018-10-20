function [EstErrK, EstErrUs, EstErrU] = HybridSimplex(SphericalFlag)

% Optimal State Estimation, by Dan Simon
%
% Hybrid unscented Kalman filter example with simplex sigma points or spherical sigma points.
% Track a body falling through the atmosphere.
% This example is taken from [Jul00], which was based on [Ath68].
% INPUTS:
%    SphericalFlag = flag indicating whether or not to use spherical sigma
%    points (if false, then use simplex sigma points)
% OUTPUTS:
%    EstErrK = Kalman filter RMS est error for position, velocity, ballistic coefficient
%    EstErrUs = Standard UKF RMS est error for position, velocity, ballistic coefficient
%    EstErrU = Minimal (spherical or simplex) UKF RMS est error for position, velocity, ballistic coefficient

if ~exist('SphericalFlag', 'var')
    SphericalFlag = true;
end

rho0 = 2; % lb-sec^2/ft^4
g = 32.2; % ft/sec^2
k = 2e4; % ft
R = 10^4; % measurement noise variance (ft^2)
Q = diag([0 0 0]); % process noise covariance
M = 10^5; % horizontal range of position sensor
a = 10^5; % altitude of position sensor

x = [3e5; -2e4; 1000];
xhat = 1.01 * [3e5; -2e4; 1000]; % EKF estimate
xhatukf = xhat; % simplex or spherical UKF estimate
xhatukfs = xhat; % standard UKF estimate

P = diag(0.1*abs(x));
Pukf = P;
Pukfs = P;

T = 0.5; % measurement time step
randn('state',sum(100*clock)); % random number generator seed

tf = 30; % simulation length (seconds)
dt = 0.001; % time step for integration (seconds)
xArray = x;
xhatArray = xhat;
xhatukfArray = xhatukf;
xhatukfsArray = xhatukfs;

N = length(x);

% Generate the weights for the simplex or spherical UKF.
if SphericalFlag
    W(1) = 0;
    W(2) = (1 - W(1)) / (N + 1);
    W(3) = W(2);
    W(4) = W(3);
    W(5) = W(4);
else
    W(1) = 0;
    W(2) = (1 - W(1)) / 2^N;
    W(3) = W(2);
    W(4) = 2 * W(3);
    W(5) = 2 * W(4);
end
% Generate the sigma points for the simplex or spherical UKF.
% Initialization
sigma = zeros(3,5);
sigma(1,1) = 0;
sigma(1,2) = -1 / sqrt(2 * W(2));
sigma(1,3) = 1 / sqrt(2 * W(2));
% j = 2 iteration
j = 2;
if SphericalFlag
    sigma(1:j,1) = [sigma(1:j-1,1) ; 0];
    sigma(1:j,2) = [sigma(1:j-1,2) ; -1/sqrt(j*(j+1)*W(2))];
    sigma(1:j,3) = [sigma(1:j-1,3) ; -1/sqrt(j*(j+1)*W(2))];
    sigma(1:j,4) = [zeros(j-1,1) ; j/sqrt(j*(j+1)*W(2))];
else
    sigma(1:j,1) = [sigma(1:j-1,1) ; 0];
    sigma(1:j,2) = [sigma(1:j-1,2) ; -1/sqrt(2*W(j+2))];
    sigma(1:j,3) = [sigma(1:j-1,3) ; -1/sqrt(2*W(j+2))];
    sigma(1:j,4) = [zeros(j-1,1) ; 1/sqrt(2*W(j+2))];
end
% j = 3 iteration
j = 3;
if SphericalFlag
    sigma(1:j,1) = [sigma(1:j-1,1) ; 0];
    sigma(1:j,2) = [sigma(1:j-1,2) ; -1/sqrt(j*(j+1)*W(2))];
    sigma(1:j,3) = [sigma(1:j-1,3) ; -1/sqrt(j*(j+1)*W(2))];
    sigma(1:j,4) = [sigma(1:j-1,4) ; -1/sqrt(j*(j+1)*W(2))];
    sigma(1:j,5) = [zeros(j-1,1) ; j/sqrt(j*(j+1)*W(2))];
else
    sigma(1:j,1) = [sigma(1:j-1,1) ; 0];
    sigma(1:j,2) = [sigma(1:j-1,2) ; -1/sqrt(2*W(j+2))];
    sigma(1:j,3) = [sigma(1:j-1,3) ; -1/sqrt(2*W(j+2))];
    sigma(1:j,4) = [sigma(1:j-1,4) ; -1/sqrt(2*W(j+2))];
    sigma(1:j,5) = [zeros(j-1,1) ; 1/sqrt(2*W(j+2))];
end

Ws = ones(2*N,1) / 2 / N; % standard UKF weights

for t = T : T : tf
    % Simulate the system.
    for tau = dt : dt : T
        xdot(1,1) = x(2);
        xdot(2,1) = rho0 * exp(-x(1)/k) * x(2)^2 / 2 / x(3) - g;
        xdot(3,1) = 0;
        xdot = xdot + sqrt(dt * Q) * [randn; randn; randn];
        x = x + xdot * dt;
    end
    % Simulate the noisy measurement.
    z = sqrt(M^2 + (x(1)-a)^2) + sqrt(R) * randn;
    % Simulate the continuous-time part of the Kalman filter (time update).
    for tau = dt : dt : T
        xhatdot(1,1) = xhat(2);
        xhatdot(2,1) = rho0 * exp(-xhat(1)/k) * xhat(2)^2 / 2 / xhat(3) - g;
        xhatdot(3,1) = 0;
        xhat = xhat + xhatdot * dt;
        temp = rho0 * exp(-xhat(1)/k) * xhat(2) / xhat(3);
        F = [0 1 0; -temp * xhat(2) / 2 / k temp ...
            -temp * xhat(2) / 2 / xhat(3); ...
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
    P = (eye(N) - K * H) * P * (eye(N) - K * H)' + K * R * K';
    % Start of simplex / spherical UKF - matrix square root calculation
    [root,p] = chol(Pukf);
    for i = 1 : N+2
        xbreve(:,i) = xhatukf + root * sigma(:,i);
    end
    % Simplex / spherical UKF time update
    for i = 1 : N+2
        for tau = dt : dt : T
            xbrevedot(1,1) = xbreve(2,i);
            xbrevedot(2,1) = rho0 * exp(-xbreve(1,i)/k) * xbreve(2,i)^2 / 2 / xbreve(3,i) - g;
            xbrevedot(3,1) = 0;
            xbreve(:,i) = xbreve(:,i) + xbrevedot * dt;
        end
    end
    xhatukf = zeros(N,1);
    for i = 1 : N+2
        xhatukf = xhatukf + W(i) * xbreve(:,i);
    end
    Pukf = Q;
    for i = 1 : N+2
        Pukf = Pukf + W(i) * (xbreve(:,i) - xhatukf) * (xbreve(:,i) - xhatukf)';
    end
    % Simplex / spherical UKF measurement update
    for i = 1 : N+2
        zukf(:,i) = sqrt(M^2 + (xbreve(1,i)-a)^2);
    end
    zhat = 0;
    for i = 1 : N+2
        zhat = zhat + W(i) * zukf(:,i);
    end
    Py = R;
    Pxy = zeros(N,1);
    for i = 1 : N+2
        Py = Py + W(i) * (zukf(:,i) - zhat) * (zukf(:,i) - zhat)';
        Pxy = Pxy + W(i) * (xbreve(:,i) - xhatukf) * (zukf(:,i) - zhat)';
    end
    Kukf = Pxy * inv(Py);
    xhatukf = xhatukf + Kukf * (z - zhat);
    Pukf = Pukf - Kukf * Py * Kukf';
    % Start of standard UKF - generate the UKF sigma points.
    [root,p] = chol(N*Pukfs);
    for i = 1 : N
        sigmas(:,i) = xhatukfs + root(i,:)';
        sigmas(:,i+N) = xhatukfs - root(i,:)';
    end
    for i = 1 : 2*N
        xbreve(:,i) = sigmas(:,i);
    end
    % Standard UKF time update
    for i = 1 : 2*N
        for tau = dt : dt : T
            xbrevedot(1,1) = xbreve(2,i);
            xbrevedot(2,1) = rho0 * exp(-xbreve(1,i)/k) * xbreve(2,i)^2 / 2 / xbreve(3,i) - g;
            xbrevedot(3,1) = 0;
            xbreve(:,i) = xbreve(:,i) + xbrevedot * dt;
        end
    end
    xhatukfs = zeros(N,1);
    for i = 1 : 2*N
        xhatukfs = xhatukfs + Ws(i) * xbreve(:,i);
    end
    Pukfs = Q;
    for i = 1 : 2*N
        Pukfs = Pukfs + Ws(i) * (xbreve(:,i) - xhatukfs) * (xbreve(:,i) - xhatukfs)';
    end
    % Standard UKF measurement update
    for i = 1 : 2*N
        zukf(:,i) = sqrt(M^2 + (xbreve(1,i)-a)^2);
    end
    zhat = 0;
    for i = 1 : 2*N
        zhat = zhat + Ws(i) * zukf(:,i);
    end
    Py = R;
    Pxy = zeros(N,1);
    for i = 1 : 2*N
        Py = Py + Ws(i) * (zukf(:,i) - zhat) * (zukf(:,i) - zhat)';
        Pxy = Pxy + Ws(i) * (xbreve(:,i) - xhatukfs) * (zukf(:,i) - zhat)';
    end
    Kukf = Pxy * inv(Py);
    xhatukfs = xhatukfs + Kukf * (z - zhat);
    Pukfs = Pukfs - Kukf * Py * Kukf';
    % Save data for plotting.
    xArray = [xArray x];
    xhatArray = [xhatArray xhat];
    xhatukfArray = [xhatukfArray xhatukf];
    xhatukfsArray = [xhatukfsArray xhatukfs];
end

EstErrK = std(xArray' - xhatArray');
EstErrU = std(xArray' - xhatukfArray');
EstErrUs = std(xArray' - xhatukfsArray');
disp(['RMS altitude est err = ', num2str(EstErrK(1)), ' (EKF), ', num2str(EstErrUs(1)), ...
    ' (Std UKF), ', num2str(EstErrU(1)), ' (Minimal UKF)']);
disp(['RMS velocity est err = ', num2str(EstErrK(2)), ' (EKF), ', num2str(EstErrUs(2)), ...
    ' (Std UKF), ', num2str(EstErrU(2)), ' (Minimal UKF)']);
disp(['RMS ball coeff est err = ', num2str(EstErrK(3)), ' (EKF), ', num2str(EstErrUs(3)), ...
    ' (Std UKF), ', num2str(EstErrU(3)), ' (Minimal UKF)']);
