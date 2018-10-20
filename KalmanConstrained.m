function KalmanConstrained

% function KalmanConstrained
% This m-file simulates a vehicle tracking problem.
% The vehicle state is estimated with a Kalman filter.
% In addition, with the a priori knowledge that the vehicle is on
% a particular road, the vehicle state is estimated with a
% constrained Kalman filter.
% The state consists of the north and east position, and the
% north and east velocity of the vehicle.
% The measurement consists of north and east positions.

tf = 300; % final time (seconds)
T = 3; % time step (seconds)

Q = diag([4, 4, 1, 1]); % Process noise covariance (m, m, m/sec, m/sec)
Qsqrt = sqrt(Q);

R = diag([900, 900]); % Measurement noise covariance (m, m)
Rsqrt = sqrt(R);

% Measurement noise covariance for perfect measurement formulation.
R1 = diag([900, 900, 0, 0]);
R1sqrt = sqrt(R1);

theta = 0.9 * pi; % heading angle (measured CCW from east)
tantheta = tan(theta);

% Define the initial state x, initial unconstrained filter estimate xhat,
% and initial constrained Kalman filter estimate xtilde.
x = [0; 0; tantheta; 1] * 100;
xhat = x; % Unconstrained Kalman filter
xhat1 = x; % Kalman filter with perfect measurements
xtilde = x; % Constrained Kalman filter (W=I)
xtildeP = x; % Constrained Kalman filter (W=inv(P))
u = 1; % control input (alternates between +1 and -1)

% Initial estimation error covariance
P = diag([R(1,1), R(2,2), Q(1,1), Q(2,2)]);
% Initial estimation error covariance for perfect measurement formulation
P1 = P;

D = [1 -tantheta 0 0; 0 0 1 -tantheta]; % State constraint matrix.
F = [1 0 T 0 ; 0 1 0 T ; 0 0 1 0 ; 0 0 0 1]; % system matrix
B = [0; 0; T*sin(theta); T*cos(theta)]; % input matrix
H = [1 0 0 0 ; 0 1 0 0]; % measurement matrix
H1 = [H ; D]; % measurement matrix for perfect measurement formulation.

% Initialize arrays for saving data for plotting.
xarray = x;
xhatarray = xhat;
Constrarray = D * xhat;
xhat1array = xhat1;
Constr1array = D * xhat1;
ConstrTildearray = D * xtilde;
ConstrTildeParray = D * xtildeP;
xtildearray = xtilde;
xtildeParray = xtildeP;

randn('state', sum(100*clock)); % initialize random number generator

% Begin the simulation.
for t = T : T : tf
    % Set the known input.
    if u == 1
        if (x(3) > 30) | (x(4) > 30)
            u = -1;
        end
    else
        if (x(3) < 5) | (x(4) < 5)
            u = 1;
        end
    end
    % Simulate the system.
    x = F*x + B*u + Qsqrt*randn(size(x));
    % Constrain the vehicle (i.e., the true state) to the straight road.
    x(1) = x(2) * tantheta;
    x(3) = x(4) * tantheta;
    % Get the measurement.
    y = H * x + Rsqrt * randn(size(Rsqrt,1),1);
    % Form the measurement vector for the perfect measurement formulation.
    y1 = [y ; 0 ; 0];
    P = F * P * F' + Q; % time update for standard Kalman filter
    P1 = F * P1 * F' + Q; % time update for perfect measurement filter
    K = P * H' * inv(H * P * H' + R); % Kalman gain for standard filter
    K1 = P1 * H1' * inv(H1 * P1 * H1' + R1); % Kalman gain for perfect measurement filter
    xhat = F*xhat + B*u; % time update for standard Kalman filter
    xhat1 = F*xhat1 + B*u; % time update for perfect measurement filter
    xhat = xhat + K * (y - H * xhat); % measurement update for standard Kalman filter
    xhat1 = xhat1 + K1 * (y1 - H1 * xhat1); % measurement update for perfect measurement filter
    P = P - K * H * P; % measurement update for standard Kalman filter
    P1 = P1 - K1 * H1 * P1; % measurement update for perfect measurement filter
    % Find the constrained Kalman filter estimates.
    xtilde = xhat - D' * inv(D*D') * D * xhat;
    xtildeP = xhat - P * D' * inv(D*P*D') * D * xhat;
    % Save data in arrays.
    xhatarray = [xhatarray xhat];
    xhat1array = [xhat1array xhat1];
    xtildearray = [xtildearray xtilde];
    xtildeParray = [xtildeParray xtildeP];
    ConstrErr = D * xhat;
    Constrarray = [Constrarray ConstrErr];
    Constr1Err = D * xhat1;
    Constr1array = [Constr1array Constr1Err];
    ConstrTilde = D * xtilde;
    ConstrTildearray = [ConstrTildearray ConstrTilde];
    ConstrTildeP = D * xtildeP;
    ConstrTildeParray = [ConstrTildeParray ConstrTildeP];
    xarray = [xarray x];
end

% Compute averages.
EstError = xarray - xhatarray;
EstError = sqrt(EstError(1,:).^2 + EstError(2,:).^2);
EstError = mean(EstError);
disp(['RMS Unconstrained Estimation Error = ', num2str(EstError)]);

EstError1 = xarray - xhat1array;
EstError1 = sqrt(EstError1(1,:).^2 + EstError1(2,:).^2);
EstError1 = mean(EstError1);
disp(['RMS Estimation Error (Perfect Meas) = ', num2str(EstError1)]);

EstErrorConstr = xarray - xtildearray;
EstErrorConstr = sqrt(EstErrorConstr(1,:).^2 + EstErrorConstr(2,:).^2);
EstErrorConstr = mean(EstErrorConstr);
disp(['RMS Constrained Estimation Error (W=I) = ', num2str(EstErrorConstr)]);

EstErrorConstrP = xarray - xtildeParray;
EstErrorConstrP = sqrt(EstErrorConstrP(1,:).^2 + EstErrorConstrP(2,:).^2);
EstErrorConstrP = mean(EstErrorConstrP);
disp(['RMS Constrained Estimation Error (W=inv(P)) = ', num2str(EstErrorConstrP)]);

disp(' ');

Constr = sqrt(Constrarray(1,:).^2 + Constrarray(2,:).^2);
Constr = mean(Constr);
disp(['Average Constraint Error (Unconstrained) = ', num2str(Constr)]);

Constr1 = sqrt(Constr1array(1,:).^2 + Constr1array(2,:).^2);
Constr1 = mean(Constr1);
disp(['Average Constraint Error (Perfect Meas) = ', num2str(Constr1)]);

ConstrTilde = sqrt(ConstrTildearray(1,:).^2 + ConstrTildearray(2,:).^2);
ConstrTilde = mean(ConstrTilde);
disp(['Average Constraint Error (W=I) = ', num2str(ConstrTilde)]);

ConstrTildeP = sqrt(ConstrTildeParray(1,:).^2 + ConstrTildeParray(2,:).^2);
ConstrTildeP = mean(ConstrTildeP);
disp(['Average Constraint Error (W=inv(P)) = ', num2str(ConstrTildeP)]);

% Plot data.
close all;
t = 0 : T : tf;
figure;
plot(t, xarray(1, :), 'r:', t, xarray(2, :), 'b-');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('seconds'); ylabel('position (meters)');
legend('north position', 'east position');

figure;
plot(t, xarray(1, :) - xhatarray(1, :), 'r:', t, xarray(2, :) - xhatarray(2, :), 'b-');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('Position Estimation Error (Unconstrained Filter)');
xlabel('seconds'); ylabel('estimation error');
legend('north position estimation error', 'east position estimation error');

figure;
plot(t, xarray(1, :) - xhat1array(1, :), 'r:', t, xarray(2, :) - xhat1array(2, :), 'b-');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('Position Estimation Error (Perfect Measurements)');
xlabel('seconds'); ylabel('estimation error');
legend('north position estimation error', 'east position estimation error');

figure;
plot(t, xarray(1, :) - xtildearray(1, :), 'r:', t, xarray(2, :) - xtildearray(2, :), 'b-');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('Position Estimation Error (Constrained, W=I)');
xlabel('seconds'); ylabel('estimation error');
legend('north position estimation error', 'east position estimation error');

figure;
plot(t, xarray(1, :) - xtildeParray(1, :), 'r:', t, xarray(2, :) - xtildeParray(2, :), 'b-');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('Position Estimation Error (Constrained, W=inv(P))');
xlabel('seconds'); ylabel('estimation error');
legend('north position estimation error', 'east position estimation error');

figure;
plot(t, abs(xarray(1, :) - xhatarray(1, :)), 'r:', t, abs(xarray(1, :) - xtildearray(1, :)), 'b-');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('seconds'); ylabel('north position estimation error');
legend('unconstrained', 'constrained');
