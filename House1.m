function [W, T] = House1(A)

% Compute the Householder transformation T A = [ W ; 0 ] where A is a given
% 2n x n matrix, and T is an orthogonal 2n x 2n matrix, and W is an n x n
% matrix.

n = size(A,2);
if size(A,1) ~= 2*n then
    disp('A must be a 2n x n matrix');
    return
end

T = eye(2*n);

for k = 1 : n
    temp = 0;
    for i = k : 2*n
        temp = temp + A(i,k)^2;
    end
    signAkk = sign(A(k,k));
    if (signAkk == 0)
        signAkk = 1;
    end;
    sigma = signAkk * sqrt(temp);
    if abs(sigma) < 100 * eps
        disp('The input matrix is not full rank');
        return;
    end
    beta = 1 / sigma / (sigma + A(k,k));
    u = zeros(2*n, 1);
    u(k) = sigma + A(k,k);
    for i = k+1 : 2*n
        u(i) = A(i,k);
    end
    y = zeros(n, 1);
    y(k) = 1;
    for i = k+1 : n
        y(i) = beta * u' * A(:, i);
    end
    A = A - u * y';
    T = (eye(2*n) - beta * u * u') * T;
end
W = A(1:n, 1:n);