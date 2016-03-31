D = 1;
u = 0;
L = 10;
N = 100;

theta0 = 1;
theta1 = 0;
Q = -0.05;

dx = L/(N-1);

A = zeros(N-2);
b = ones(N-2,1)*Q;
% x = zeros(N-2,1);

for ii = 1:N-2
    A(ii,ii) = (2*D+u*dx)/dx^2;
    if ii~=1
        A(ii,ii-1) = -(D+u*dx)/dx^2;
    end
    if ii~=N-2
        A(ii,ii+1) = -(D+u*dx)/dx^2;
    end
    
end

b(1) = b(1) + (D+u*dx)/dx^2*theta0;
b(end) = b(end) + (D+u*dx)/dx^2*theta1;

x = A\b;

plot(x)

