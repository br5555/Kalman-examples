function UnscentedEx

x = [0; 1]; % x,y coordinates of target
r = sqrt(x(1)^2 + x(2)^2);
theta = atan2(x(2), x(1)); % radians
rand('state',sum(100*clock));

% uniform distribution
rtilde = 0.01;
thetatilde = 0.35; % radians
sigmar = rtilde / sqrt(3);
sigmatheta = thetatilde / sqrt(3);
zarr = [r; theta];

for k = 1 : 300
    z = [r + 2 * rtilde * (rand-0.5); theta + 2 * thetatilde * (rand-0.5)];
    zarr = [zarr z];
end

close all;
figure;
xarr = zarr(1,:) .* cos(zarr(2,:));
yarr = zarr(1,:) .* sin(zarr(2,:));
plot(xarr, yarr, '.');
xbar = 0;
ybar = sin(thetatilde) / thetatilde;
hold;
plot(xbar, ybar, 'ro');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('y_1'); ylabel('y_2');

figure;
plot(xarr, yarr, '.');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('y_1'); ylabel('y_2');
hold;
plot(xbar, ybar, 'ro');

Px = 0.5 * (1 + sigmar^2) * (1 - sin(2 * thetatilde) / 2 / thetatilde);
Px = sqrt(Px);
Py = 0.5 * (1 + sigmar^2) * (1 + sin(2 * thetatilde) / 2 / thetatilde) - (sin(thetatilde) / thetatilde)^2;
Py = sqrt(Py);

% plot the true uncertainty ellipse
xarr = [];
yarr1 = [];
yarr2 = [];
for x = xbar-Px : Px / 100 : xbar+Px
    xarr = [xarr x];
    y = ybar + Py * sqrt(1 - ((x - xbar) / Px)^2);
    yarr1 = [yarr1 y];
    y = ybar - Py * sqrt(1 - ((x - xbar) / Px)^2);
    yarr2 = [yarr2 y];
end
plot(xbar, ybar, 'ro');
plot(xarr, yarr1, 'r');
plot(xarr, yarr2, 'r');
    
% plot the linearized uncertainty ellipse
xbar = 0;
ybar = 1;
Px = sigmatheta;
Py = sigmar;
xarr = [];
yarr1 = [];
yarr2 = [];
for x = xbar-Px : Px / 100 : xbar+Px
    xarr = [xarr x];
    y = ybar + Py * sqrt(1 - ((x - xbar) / Px)^2);
    yarr1 = [yarr1 y];
    y = ybar - Py * sqrt(1 - ((x - xbar) / Px)^2);
    yarr2 = [yarr2 y];
end
plot(xbar, ybar, 'kx');
plot(xarr, yarr1, 'k');
plot(xarr, yarr2, 'k');

% unscented transformation
clear x;
x(:,1) = [1 + sigmar * sqrt(2) ; pi/2];
x(:,2) = [1 ; pi/2 + sigmatheta * sqrt(2)];
x(:,3) = [1 - sigmar * sqrt(2) ; pi/2];
x(:,4) = [1 ; pi/2 - sigmatheta * sqrt(2)];
z(:,1) = [0 ; 1 + sigmar * sqrt(2)];
z(:,2) = [cos(pi/2+sigmatheta*sqrt(2)) ; sin(pi/2+sigmatheta*sqrt(2))];
z(:,3) = [0 ; 1 - sigmar * sqrt(2)];
z(:,4) = [cos(pi/2-sigmatheta*sqrt(2)) ; sin(pi/2-sigmatheta*sqrt(2))];
zbar = 1/4 * (z(:,1) + z(:,2) + z(:,3) + z(:,4));
Sigmaz = 0;
for i = 1 : 4
    Sigmaz = Sigmaz + 1/4 * (z(:,i) - zbar) * (z(:,i) - zbar)';
end
Px = sqrt(Sigmaz(1,1));
Py = sqrt(Sigmaz(2,2));
% plot the unscented transformation
xarr = [];
yarr1 = [];
yarr2 = [];
for x = zbar(1)-Px : Px / 100 : zbar(1)+Px
    xarr = [xarr x];
    y = zbar(2) + Py * sqrt(1 - ((x - zbar(1)) / Px)^2);
    yarr1 = [yarr1 y];
    y = zbar(2) - Py * sqrt(1 - ((x - zbar(1)) / Px)^2);
    yarr2 = [yarr2 y];
end

plot(zbar(1), zbar(2), 'mx');
plot(xarr, yarr1, 'm');
plot(xarr, yarr2, 'm');
