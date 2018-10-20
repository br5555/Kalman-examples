function LinearSimEx1

% Compare rectangular, trapezoidal, and fourth order Runge Kutta integration.
% Integrate xdot = cos(t).

tf = 1;
TArray = [0.1 0.05 0.025];
for i = 1 : length(TArray)
    T = TArray(i);
    % Rectangular integration
    j = 1;
    xRect(j) = 0;
    for t = 0 : T : tf - T + T/10
        j = j + 1;
        xRect(j) = xRect(j-1) + cos(t) * T;
    end
    % Trapezoidal integration
    j = 1;
    xTrap(j) = 0;
    for t = 0 : T : tf - T + T/10
        j = j + 1;
        xTrap(j) = xTrap(j-1) + (cos(t) + cos(t+T)) * T / 2;
    end
    % Fourth order Runge Kutta integration
    j = 1;
    xRK(j) = 0;
    for t = 0 : T : tf - T + T/10
        t1 = t + T/2;
        d1 = cos(t);
        d2 = cos(t1);
        d3 = cos(t1);
        d4 = cos(t+T);
        j = j + 1;
        xRK(j) = xRK(j-1) + (d1 + 2 * d2 + 2 * d3 + d4) * T / 6;
    end
    errRect = 100 * abs(xRect(end) - sin(tf)) / sin(tf);
    errTrap = 100 * abs(xTrap(end) - sin(tf)) / sin(tf);
    errRK = 100 * abs(xRK(end) - sin(tf)) / sin(tf);
    disp(['T = ', num2str(T), ': ', num2str(errRect), ', ', num2str(errTrap), ', ', num2str(errRK)]);
end