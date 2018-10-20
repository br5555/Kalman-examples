function ContEx(a1, a2, q11, q12, q22, r1, r2)

% Continuous Kalman filter example for a two-state problem.

% Try different initial conditions.
p11 = 2; p12 = 1; p22 = 2;
p11 = 0; p12 = 0; p22 = 0;

PArr = [];
dt = 0.01;
tf = 3;
% Plot the time varying solution.
for t = 0 : dt : tf
    p11dot = 2 * a1 * p11 - p11^2 / r1 - p12^2 / r2 + q11;
    p12dot = (a1 + a2) * p12 - p11 * p12 / r1 - p12 * p22 / r2 + q12;
    p22dot = 2 * a2 * p22 - p12^2 / r1 - p22^2 / r2 + q22;
    p11 = p11 + p11dot * dt;
    p12 = p12 + p12dot * dt;
    p22 = p22 + p22dot * dt;
    PArr = [PArr ; p11 p12 p22];
end
close all
t = 0 : dt : tf;
plot(t, PArr(:, 1), t, PArr(:, 2), t, PArr(:, 3));
grid;
legend('p11', 'p12', 'p22');
% Compute the steady state solution.
Cond1 = (a1 ~= a2) && (q12 ~= 0);
Cond2 = (a1 == a2) && (a1 < 0) && (q12 ~= 0) && (q11*q22-q12*q12 == 0);
Cond3 = (a1 == a2) && (a1 > 0) && (q12 ~= 0) && (q11*q22-q12*q12 == 0);
if Cond1 || Cond3
    gamma1 = q11 / r1 + a1^2;
    gamma2 = q22 / r2 + a2^2;
    p12 = q12 / ( gamma1 + gamma2 + 2 * ( gamma1 * gamma2 - q12^2 / r1 / r2 )^(1/2) )^(1/2);
    p11 = r1 * ( a1 + ( gamma1 - p12^2 / r1 / r2 )^(1/2) );
    p22 = r2 * ( a2 + ( gamma2 - p12^2 / r1 / r2 )^(1/2) );
    disp(['p11 = ', num2str(p11), ', p12 = ', num2str(p12), ', p22 = ', num2str(p22)]);
    lambda = eig([p11 p12; p12 p22]);
    disp(['Eigenvalues of P = ', num2str(lambda(1)), ', ', num2str(lambda(2))]);
end
if Cond2 || Cond3
    gamma3 = -a1 + ( a1^2 + q11 / r1 + q22 / r2 )^(1/2);
    p11 = q11 / gamma3;
    p22 = q22 / gamma3;
    p12 = q12 / gamma3;
    disp(['p11 = ', num2str(p11), ', p12 = ', num2str(p12), ', p22 = ', num2str(p22)]);
    lambda = eig([p11 p12; p12 p22]);
    disp(['Eigenvalues of P = ', num2str(lambda(1)), ', ', num2str(lambda(2))]);
end
if ~Cond1 && ~Cond2 && ~Cond3 
    disp(['p11 = ', num2str(p11), ', p12 = ', num2str(p12), ', p22 = ', num2str(p22)]);
    lambda = eig([p11 p12; p12 p22]);
    disp(['Eigenvalues of P = ', num2str(lambda(1)), ', ', num2str(lambda(2))]);
end