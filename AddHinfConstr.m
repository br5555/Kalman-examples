function [ErrKarray, ErrKCarray, ErrHinfarray, ErrHinfCarray] = AddHinfConstr(g, T, tf)

% function AddHinfConstr
% This m-file simulates a vehicle tracking problem.
% The vehicle state is estimated with a minimax filter.
% In addition, with the a priori knowledge that the vehicle is on
% a particular road, the vehicle state is estimated with a 
% constrained minimax filter.
% This m-file also simulates a Kalman filter and constrained
% Kalman filter so you can compare results.
% The state consists of the north and east position, and the
% north and east velocity of the vehicle.
% The measurement consists of north and east position.
% For further details see the web site 
% http://www.csuohio.edu/simond/minimaxconstrained/.
% INPUTS
%   g = gamma (I suggest 40)
%   T = time step in seconds (I suggest 1)
%   tf = final time in seconds (I suggest 120)
% OUTPUTS
%   ErrKarray = time varying array of error of Kalman unconstrained state estimate
%   ErrKCarray = time varying array of error of Kalman constrained state estimate
%   ErrHinfarray = time varying array of error of Minimax unconstrained state estimate
%   ErrHinfCarray = time varying array of error of Minimax constrained state estimate

if ~exist('g', 'var')
    g = 40;
end
if ~exist('T', 'var')
    T = 1;
end
if ~exist('tf', 'var')
    tf = 120;
end

Q = diag([4, 4, 1, 1]); % Process noise covariance (m, m, m/sec, m/sec)
Qsqrt = sqrt(Q);

R = diag([900, 900]); % Measurement noise covariance (m, m)
Rsqrt = sqrt(R);

theta = pi / 3; % heading angle (measured CCW from east)
tantheta = tan(theta);

% Define the initial state x, initial unconstrained Kalman filter estimate xhat,
% and initial constrained Kalman filter estimate xtilde.
x = [0; 0; tantheta; 1] * 100;
xhat = x;
xtilde = x;
P = diag([R(1,1), R(2,2), Q(1,1), Q(2,2)]); % Initial estimation error covariance

% AccelDecelFlag is used to simulate the vehicle alternately accelerating and
% decelerating, as if in traffic.
AccelDecelFlag = 1;

% System matrix.
A = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];

% Input matrix.
B = [0; 0; T*sin(theta); T*cos(theta)];

% Normalized measurement matrix.
C = inv(Rsqrt) * [1 0 0 0; 0 1 0 0];

% State constraint matrices.
D = [1 -tantheta 0 0; 0 0 1 -tantheta];
% Normalize D so that D*D'=I.
D = D / sqrt(1 + tantheta^2);
V = D' *  D;
d = [0; 0];

% Initialize arrays for saving data for plotting.
xarray = [];
xhatarray = [];
xtildearray = [];
randn('state', sum(100*clock));

% Minimax initialization.
% Make sure that xtildeinf satisfies the state constraint.
Qbar = P; 
Qtilde = P;
xhatinf = x;
xtildeinf = x;
xhatinfarray = [];
xtildeinfarray = [];

for t = T : T : tf
    
    % Get the noise-corrupted measurement z.
    z = C * x;
    MeasErr = randn(size(z));
    z = z + MeasErr;
    
    % Set the known input u.
    if AccelDecelFlag == 1
        if (x(3) > 30) | (x(4) > 30)
            AccelDecelFlag = -1;
        end
    else
        if (x(3) < 5) | (x(4) < 5)
            AccelDecelFlag = 1;
        end
    end
    u = 1 * AccelDecelFlag;
    
    % Run the unconstrained minimax filter.
    Pinf = inv(eye(4) - Qbar / g / g + Qbar * C' * C) * Qbar;
    Qbar = A * Pinf * A' + Q;
    K = A * Pinf * C';
    headinghat = atan2(xhatinf(3), xhatinf(4));  
    Bhat = [0; 0; T*sin(headinghat); T*cos(headinghat)];
    xhatinf = A * xhatinf + Bhat * u + K * (z - C * xhatinf);
    xhatinfarray = [xhatinfarray xhatinf];
    
    % Run the constrained minimax filter.
    Ptilde = inv(eye(4) - Qtilde / g / g + Qtilde * C' * C) * Qtilde;
    Qtilde = (eye(4) - V) * A * Ptilde * A' * (eye(4) - V) + Q;
    lambda = eig(eye(4) - Qtilde / g / g);
    for i = 1 : 4
        if ~isreal(lambda(i)) || lambda(i) <= 0
            disp(['positive semidefinite condition fails - t = ', num2str(t)]);
            return;
        end
    end
    K = (eye(4) - V) * A * Pinf * C';
    headingtilde = atan2(xtildeinf(3), xtildeinf(4));
    Btilde = [0; 0; T*sin(headingtilde); T*cos(headingtilde)];
    xtildeinf = A * xtildeinf + Btilde * u + K * (z - C * xtildeinf);
    xtildeinfarray = [xtildeinfarray xtildeinf];
    
    % Run the unconstrained Kalman filter.
    K = A * P * C' * inv(C * P * C' + eye(size(R)));
    % Update the state estimation error covariance.
    P = (A * P - K * C * P) * A' + Q;   
    % Estimate the heading on the basis of the state estimate.
    headinghat = atan2(xhat(3), xhat(4));  
    Bhat = [0; 0; T*sin(headinghat); T*cos(headinghat)];
    xhat = A * xhat + Bhat * u + K * (z - C * xhat);
    xhatarray = [xhatarray xhat];
    % Find the constrained Kalman filter estimate.
    xtilde = xhat - D' * inv(D*D') * (D * xhat - d);
    xtildearray = [xtildearray xtilde];
    % Simulate the system dynamics.
    x = A * x + B * u + Qsqrt * randn(size(x));
    
    % Uncomment the following line to add unmodeled process noise.
    x = x + [0; 0; 1; 1];
    
    % Constrain the vehicle (i.e., the true state) to the straight road.
    if abs(x(1) - tantheta * x(2)) > 2
        x(2) = (x(2) + x(1) * tantheta) / (1 + tantheta^2);
        x(1) = x(2) * tantheta;
    end
    if abs(x(3) - tantheta * x(4)) > 0.2
        x(4) = (x(4) + x(3) * tantheta) / (1 + tantheta^2);
        x(3) = x(4) * tantheta;
    end
    xarray = [xarray x];
end

% Compute average position estimation errors.
EstError = xarray - xhatarray;
EstError = sqrt(EstError(1,:).^2 + EstError(2,:).^2);
EstError = mean(EstError);
disp(['Average Kalman Unconstrained Position Estimation Error = ', num2str(EstError)]);

EstErrorConstr = xarray - xtildearray;
EstErrorConstr = sqrt(EstErrorConstr(1,:).^2 + EstErrorConstr(2,:).^2);
EstErrorConstr = mean(EstErrorConstr);
disp(['Average Kalman Constrained Position Estimation Error (W=I) = ', num2str(EstErrorConstr)]);

EstError = xarray - xhatinfarray;
EstError = sqrt(EstError(1,:).^2 + EstError(2,:).^2);
EstError = mean(EstError);
disp(['Average Minimax Unconstrained Position Estimation Error = ', num2str(EstError)]);

EstErrorConstr = xarray - xtildeinfarray;
EstErrorConstr = sqrt(EstErrorConstr(1,:).^2 + EstErrorConstr(2,:).^2);
EstErrorConstr = mean(EstErrorConstr);
disp(['Average Minimax Constrained Position Estimation Error = ', num2str(EstErrorConstr)]);
disp(' ');

% Compute average velocity estimation errors.
EstError = xarray - xhatarray;
EstError = sqrt(EstError(3,:).^2 + EstError(4,:).^2);
EstError = mean(EstError);
disp(['Average Kalman Unconstrained Velocity Estimation Error = ', num2str(EstError)]);

EstErrorConstr = xarray - xtildearray;
EstErrorConstr = sqrt(EstErrorConstr(3,:).^2 + EstErrorConstr(4,:).^2);
EstErrorConstr = mean(EstErrorConstr);
disp(['Average Kalman Constrained Velocity Estimation Error (W=I) = ', num2str(EstErrorConstr)]);

EstError = xarray - xhatinfarray;
EstError = sqrt(EstError(3,:).^2 + EstError(4,:).^2);
EstError = mean(EstError);
disp(['Average Minimax Unconstrained Velocity Estimation Error = ', num2str(EstError)]);

EstErrorConstr = xarray - xtildeinfarray;
EstErrorConstr = sqrt(EstErrorConstr(3,:).^2 + EstErrorConstr(4,:).^2);
EstErrorConstr = mean(EstErrorConstr);
disp(['Average Minimax Constrained Velocity Estimation Error = ', num2str(EstErrorConstr)]);

close all;
t = 0 : T : tf-T;

% Plot the position errors.
figure;
plot(t, xarray(1, :), ':', t, xarray(2, :), '-');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('True Position');
xlabel('seconds');
ylabel('meters');

figure;
plot(t, xarray(1, :) - xhatarray(1, :), ':', ...
    t, xarray(2, :) - xhatarray(2, :), '-');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('Kalman Position Estimation Error (Unconstrained)');
xlabel('seconds');
ylabel('meters');

figure;
plot(t, xarray(1, :) - xtildearray(1, :), ':', ...
    t, xarray(2, :) - xtildearray(2, :), '-');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('Kalman Position Estimation Error (Constrained)');
xlabel('seconds');
ylabel('meters');

figure;
plot(t, xarray(1, :) - xhatinfarray(1, :), ':', ...
    t, xarray(2, :) - xhatinfarray(2, :), '-');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('Minimax Position Estimation Error (Unconstrained)');
xlabel('seconds');
ylabel('meters');

figure;
plot(t, xarray(1, :) - xtildeinfarray(1, :), ':', ...
    t, xarray(2, :) - xtildeinfarray(2, :), '-');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('Minimax Position Estimation Error (Constrained)');
xlabel('seconds');
ylabel('meters');

% Plot the velocity errors.
figure;
plot(t, xarray(3, :), ':', t, xarray(4, :), '-');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('True Velocity');
xlabel('seconds');
ylabel('meters/sec');

figure;
plot(t, xarray(3, :) - xhatarray(3, :), ':', ...
    t, xarray(4, :) - xhatarray(4, :), '-');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('Kalman Velocity Estimation Error (Unconstrained)');
xlabel('seconds');
ylabel('meters/sec');

figure;
plot(t, xarray(3, :) - xtildearray(3, :), ':', ...
    t, xarray(4, :) - xtildearray(4, :), '-');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('Kalman Velocity Estimation Error (Constrained)');
xlabel('seconds');
ylabel('meters/sec');

figure;
plot(t, xarray(3, :) - xhatinfarray(3, :), ':', ...
    t, xarray(4, :) - xhatinfarray(4, :), '-');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('Minimax Velocity Estimation Error (Unconstrained)');
xlabel('seconds');
ylabel('meters/sec');

figure;
plot(t, xarray(3, :) - xtildeinfarray(3, :), ':', ...
    t, xarray(4, :) - xtildeinfarray(4, :), '-');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('Minimax Velocity Estimation Error (Constrained)');
xlabel('seconds');
ylabel('meters/sec');

% Compute estimation errors
ErrKarray = xarray - xhatarray; % Kalman unconstrained
ErrKCarray = xarray - xtildearray; % Kalman constrained
ErrHinfarray = xarray - xhatinfarray; % Minimax unconstrained
ErrHinfCarray = xarray - xtildeinfarray; % Minimax constrained