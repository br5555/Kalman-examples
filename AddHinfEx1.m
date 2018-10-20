function AddHinfEx1

Q = 1;
R = 1;
thetaMin = 0;
thetaMax = 1;
dtheta = 0.01;
KArr = [];
PArr = [];
for theta = thetaMin : dtheta : thetaMax
    c(1) = theta^2 - theta^4 * R;
    c(2) = Q * theta^4 * R - Q * theta^2 + R * theta^2 - 1;
    c(3) = Q * (1 - 2 * theta^2 * R);
    c(4) = Q * R;
    Pall = roots(c);
    % Find a real positive root of the ARE that results in a stable estimator.
    P = inf;
    for i = 1 : length(Pall)
        if abs(theta^2 * Pall(i) - 1) < 1e-12
            continue;
        end
        Pa = Pall(i) / (1 - theta^2 * Pall(i));
        V = R + Pa;
        Fhat = 1 - Pa / V;
        if isreal(Pall(i)) && (Pall(i) >= 0) && (Pall(i) < P) && (abs(Fhat) < 1)
            P = Pall(i);
            K = Pa / V;
        end
    end
    if P == inf
        thetaMax = theta - dtheta;
        break;
    end
    PArr = [PArr P];
    KArr = [KArr K];
end

close all;
theta = thetaMin : dtheta : thetaMax;

figure;
plot(theta, KArr); 
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('H_\infty performance bound \theta'); ylabel('Estimator gain K');

figure;
plot(theta, PArr); 
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('H_\infty performance bound \theta'); ylabel('Kalman performance bound P');
