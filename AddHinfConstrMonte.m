function AddHinfConstrMonte

% Obtain Monte Carlo based statistics for the unconstrained and constrained
% Kalman and H-infinity state estimates

nMonte = 100; % number of Monte Carlo simulations
for i = 1 : nMonte
    [ErrKarray(i,:,:), ErrKCarray(i,:,:), ErrHinfarray(i,:,:), ErrHinfCarray(i,:,:)] = AddHinfConstr;
end
% At each time step compute the RMS position and velocity estimation errors
% at each time step.
nTime = size(ErrKarray,3); % number of time steps in each simulation
ErrKPos = zeros(1, nTime); % position error of unconstrained Kalman filter
ErrKVel = ErrKPos; % velocity error of unconstrained Kalman filter
ErrKCPos = ErrKPos; % position error of constrained Kalman filter
ErrKCVel = ErrKPos; % velocity error of constrained Kalman filter
ErrHinfPos = ErrKPos; % position error of unconstrained H-infinity filter
ErrHinfVel = ErrKPos; % velocity error of unconstrained H-infinity filter
ErrHinfCPos = ErrKPos; % position error of constrained H-infinity filter
ErrHinfCVel = ErrKPos; % velocity error constrained H-infinity filter
for t = 1 : nTime
    for i = 1 : nMonte
        ErrKPos(t) = ErrKPos(t) + ErrKarray(i,1,t)^2 + ErrKarray(i,2,t)^2;
        ErrKVel(t) = ErrKVel(t) + ErrKarray(i,3,t)^2 + ErrKarray(i,4,t)^2;
        ErrKCPos(t) = ErrKCPos(t) + ErrKCarray(i,1,t)^2 + ErrKCarray(i,2,t)^2;
        ErrKCVel(t) = ErrKCVel(t) + ErrKCarray(i,3,t)^2 + ErrKCarray(i,4,t)^2;
        ErrHinfPos(t) = ErrHinfPos(t) + ErrHinfarray(i,1,t)^2 + ErrHinfarray(i,2,t)^2;
        ErrHinfVel(t) = ErrHinfVel(t) + ErrHinfarray(i,3,t)^2 + ErrHinfarray(i,4,t)^2;
        ErrHinfCPos(t) = ErrHinfCPos(t) + ErrHinfCarray(i,1,t)^2 + ErrHinfCarray(i,2,t)^2;
        ErrHinfCVel(t) = ErrHinfCVel(t) + ErrHinfCarray(i,3,t)^2 + ErrHinfCarray(i,4,t)^2;
    end
    ErrKPos(t) = sqrt( ErrKPos(t) / nMonte ) ;
    ErrKVel(t) = sqrt( ErrKVel(t) / nMonte ) ;
    ErrKCPos(t) = sqrt( ErrKCPos(t) / nMonte ) ;
    ErrKCVel(t) = sqrt( ErrKCVel(t) / nMonte ) ;
    ErrHinfPos(t) = sqrt( ErrHinfPos(t) / nMonte ) ;
    ErrHinfVel(t) = sqrt( ErrHinfVel(t) / nMonte ) ;
    ErrHinfCPos(t) = sqrt( ErrHinfCPos(t) / nMonte ) ;
    ErrHinfCVel(t) = sqrt( ErrHinfCVel(t) / nMonte ) ;
end

% Plot the results
close all;
t = 1 : nTime;

figure;
plot(t, ErrKPos, 'r-', t, ErrKCPos, 'b:');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('Kalman filter position estimation error');
legend('unconstrained', 'constrained');
xlabel('seconds'); ylabel('meters');

figure;
plot(t, ErrKVel, 'r-', t, ErrKCVel, 'b:');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('Kalman filter velocity estimation error');
legend('unconstrained', 'constrained');
xlabel('seconds'); ylabel('meters/s');

figure;
plot(t, ErrHinfPos, 'r-', t, ErrHinfCPos, 'b:');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('H_{\infty} filter position estimation error');
legend('unconstrained', 'constrained');
xlabel('seconds'); ylabel('meters');

figure;
plot(t, ErrHinfVel, 'r-', t, ErrHinfCVel, 'b:');
set(gca,'FontSize',12); set(gcf,'Color','White');
title('H_{\infty} filter velocity estimation error');
legend('unconstrained', 'constrained');
xlabel('seconds'); ylabel('meters/s');

disp('average RMS position estimation errors:');
disp([num2str(mean(ErrKPos)), ', unconstrained Kalman filter']);
disp([num2str(mean(ErrKCPos)), ', constrained Kalman filter']);
disp([num2str(mean(ErrHinfPos)), ', unconstrained H_\infty filter']);
disp([num2str(mean(ErrHinfCPos)), ', constrained H_\infty filter']);
disp('average RMS velocity estimation errors:');
disp([num2str(mean(ErrKVel)), ', unconstrained Kalman filter']);
disp([num2str(mean(ErrKCVel)), ', constrained Kalman filter']);
disp([num2str(mean(ErrHinfVel)), ', unconstrained H_\infty filter']);
disp([num2str(mean(ErrHinfCVel)), ', constrained H_\infty filter']);
