function MotorSim(ControlChange)

% Two-phase step motor simulation.
% This file can be used to investigate the effects of linearization.
% If ControlChange = 0 then the nonlinear and linearized simulations should
% match exactly (in steady state). As ControlChange increases, the linearized simulation
% diverges more and more from the nonlinear simulation.

if ~exist('ControlChange', 'var')
    ControlChange = 0;
end

Ra = 1.9; % Winding resistance
L = 0.003; % Winding inductance
lambda = 0.1; % Motor constant
J = 0.00018; % Moment of inertia
B = 0.001; % Coefficient of viscous friction

dt = 0.0005; % Integration step size
tf = 1; % Simulation length 

x = [0; 0; 0; 0]; % Initial state
xlin = x; % Linearized approximation of state
w = 2 * pi; % Control input frequency

dtPlot = 0.005; % How often to plot results
tPlot = -inf;

% Initialize arrays for plotting at the end of the program
xArray = [];
xlinArray = [];
tArray = [];

dx = x - xlin; % Difference between true state and linearized state

tStart = cputime;

% Begin simulation loop
for t = 0 : dt : tf
    if t >= tPlot + dtPlot
        % Save data for plotting
        tPlot = t + dtPlot - eps;
        xArray = [xArray x];
        xlinArray = [xlinArray xlin];
        tArray = [tArray t];
    end
    % Nonlinear simulation
    ua0 = sin(w*t); % nominal winding A control input
    ub0 = cos(w*t); % nominal winding B control input
    ua = (1 + ControlChange) * sin(w*t); % true winding A control input
    ub = (1 + ControlChange) * cos(w*t); % true winding B control input
    xdot = [-Ra/L*x(1) + x(3)*lambda/L*sin(x(4)) + ua/L;
        -Ra/L*x(2) - x(3)*lambda/L*cos(x(4)) + ub/L;
        -3/2*lambda/J*x(1)*sin(x(4)) + 3/2*lambda/J*x(2)*cos(x(4)) - B/J*x(3);
        x(3)];
    x = x + xdot * dt;
    x(4) = mod(x(4), 2*pi);
    % Linear simulation
    w0 = -6.2832; % nominal rotor speed
    theta0 = -6.2835 * t + 2.3679; % nominal rotor position
    ia0 = 0.3708 * cos(2*pi*(t-1.36)); % nominal winding a current
    ib0 = -0.3708 * sin(2*pi*(t-1.36)); % nominal winding b current
    du = [ua - ua0; ub - ub0];
    F = [-Ra/L 0 lambda/L*sin(theta0) w0*lambda/L*cos(theta0);
        0 -Ra/L -lambda/L*cos(theta0) w0*lambda/L*sin(theta0);
        -3/2*lambda/J*sin(theta0) 3/2*lambda/J*cos(theta0) -B/J -3/2*lambda/J*(ia0*cos(theta0)+ib0*sin(theta0));
        0 0 1 0];
    G = [1/L 0; 0 1/L; 0 0; 0 0];
    dxdot = F * dx + G * du;
    dx = dx + dxdot * dt;
    xlin = [ia0; ib0; w0; theta0] + dx;
    xlin(4) = mod(xlin(4), 2*pi);
end

TCPU = cputime - tStart;
disp(['CPU time = ', num2str(TCPU)]);

% Compare linear and nonlinear simulation.
close all;
figure; hold on;
plot(tArray,xArray(1,:),'b-','LineWidth',1.0);
plot(tArray,xlinArray(1,:),'r--','LineWidth',1.0); 
set(gca,'FontSize',12); set(gcf,'Color','White'); set(gca,'Box','on');
xlabel('Seconds'); ylabel('Winding A Current (Amps)');
legend('Nonlinear', 'Linearized');

figure; hold on;
plot(tArray,xArray(2,:),'b-','LineWidth',1.0);
plot(tArray,xlinArray(2,:),'r--','LineWidth',1.0); 
set(gca,'FontSize',12); set(gcf,'Color','White'); set(gca,'Box','on');
xlabel('Seconds'); ylabel('Winding B Current (Amps)');
legend('Nonlinear', 'Linearized');

figure; hold on;
plot(tArray,xArray(3,:),'b-','LineWidth',1.0);
plot(tArray,xlinArray(3,:),'r--','LineWidth',1.0); 
set(gca,'FontSize',12); set(gcf,'Color','White'); set(gca,'Box','on');
xlabel('Seconds'); ylabel('Rotor Speed (Rad / Sec)');
legend('Nonlinear', 'Linearized');

figure; hold on;
plot(tArray,xArray(4,:),'b-','LineWidth',1.0);
plot(tArray,xlinArray(4,:),'r--','LineWidth',1.0); 
set(gca,'FontSize',12); set(gcf,'Color','White'); set(gca,'Box','on');
xlabel('Seconds'); ylabel('Rotor Position (Radians)');
legend('Nonlinear', 'Linearized');

% Put all four plots in a single figure
figure;
set(gcf,'Color','White'); 

subplot(2,2,1); hold on;
plot(tArray,xArray(3,:),'b-','LineWidth',1.5);
plot(tArray,xlinArray(3,:),'r--','LineWidth',1.5); 
set(gca,'FontSize',12); set(gca,'Box','on');
ylabel('Speed (Rad / Sec)');
legend('Nonlinear', 'Linearized');

subplot(2,2,2); hold on;
plot(tArray,xArray(4,:),'b-','LineWidth',1.5);
plot(tArray,xlinArray(4,:),'r--','LineWidth',1.5); 
set(gca,'FontSize',12); set(gca,'Box','on');
ylabel('Position (Radians)');

subplot(2,2,3); hold on;
plot(tArray,xArray(1,:),'b-','LineWidth',1.5);
plot(tArray,xlinArray(1,:),'r--','LineWidth',1.5); 
set(gca,'FontSize',12); set(gca,'Box','on');
xlabel('Seconds'); ylabel('Current A (Amps)');

subplot(2,2,4); hold on;
plot(tArray,xArray(2,:),'b-','LineWidth',1.5);
plot(tArray,xlinArray(2,:),'r--','LineWidth',1.5); 
set(gca,'FontSize',12); set(gca,'Box','on');
xlabel('Seconds'); ylabel('Current B (Amps)');
