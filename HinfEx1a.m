function HinfEx1a

th = 0 : .01 : .99;
%th = 5 : .01 : 10;

if 1 == 1
    p1 = (1-th+sqrt((th-1).*(th-5)))/2./(1-th);
    p2 = (1-th-sqrt((th-1).*(th-5)))/2./(1-th);
    close all;
    plot(th, p1, 'r'); hold;
    plot(th, p2, 'b');
    return;
end

Bound1 = 2 * (1 - th) ./ ( 1 - th + sqrt( (th - 1) .* (th - 5) ) );
Bound2 = 2 * (1 - th) ./ ( 1 - th - sqrt( (th - 1) .* (th - 5) ) );
Limit = th - 1;

close all;
plot(th, Limit, 'k'); hold;
plot(th, Bound1, 'r');
plot(th, Bound2, 'b');