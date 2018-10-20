function HinfEx1b(SSFlag, BiasFlag)

% Compare the Kalman and H-infinity estimators for a scalar system.
% SSFlag = 1 if the steady state estimator is used.
% BiasFlag = 1 if the process noise is biased

if ~exist('SSFlag', 'var')
    SSFlag = 1;
end

if SSFlag ~= 0
    KK = (1+sqrt(5))/(3+sqrt(5));
    theta = 1/3;
    %theta = 1/10;
    PH = (1-theta+sqrt((theta-1)*(theta-5)))/2/(1-theta);
    KH = PH * inv(1 - theta * PH + PH);
end
    
kf = 20;
PK = 1;
PH = PK;
Q = 10;
R = 10;
x = 0;
xhatK = x + sqrt(PK); % Kalman filter estimate (in error by one sigma)
xhatH = xhatK; % H-infinity estimate
xArr = [x];
xhatKArr = [xhatK];
xhatHArr = [xhatH];
for k = 1 : kf
    y = x + sqrt(R) * randn;
    if SSFlag == 0
        KK = PK * inv(1 + PK) * R;
        PK = PK * inv(1 + PK) + Q;
        KH = PH * inv(1 - PH/2 + PH) * R;
        PH = PH * inv(1 - PH/2 + PH) + Q;
    end
    x = x + sqrt(Q) * randn;
    if BiasFlag 
        x = x + 10;
    end
    xhatK = xhatK + KK * (y - xhatK);
    xhatH = xhatH + KH * (y - xhatH);
    xArr = [xArr x];
    xhatKArr = [xhatKArr xhatK];
    xhatHArr = [xhatHArr xhatH];
end

k = 0 : kf;
close all; 
figure; hold on;
plot(k, xArr, 'k', 'LineWidth', 2.5);
plot(k, xhatKArr, 'r--', k, xhatHArr, 'b:');
set(gca,'FontSize',12); set(gcf,'Color','White');
legend('true state', 'Kalman estimate', 'H_{\infty} estimate');
xlabel('time'); ylabel('state value');
set(gca,'box','on');

RMSK = sqrt(norm(xArr-xhatKArr)^2/(kf+1));
RMSH = sqrt(norm(xArr-xhatHArr)^2/(kf+1));
disp(['Kalman RMS estimation error = ', num2str(RMSK)]);
disp(['H-infinity RMS estimation error = ', num2str(RMSH)]);
