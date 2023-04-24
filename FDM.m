
%% Set up System
a = 0;
b = 2;

alpha = 0;
beta = -4;

h = 0.0125;
N = (b-a)/h;

p = @(x) 2;
q = @(x) -1;
r = @(x) x*exp(x) -x;

%% Exact Solution

y = @(x)  (1/6) * (x^3) * exp(x) - (5/3) * x * exp(x) + 2 * exp(x)- x -2;

%% Get x

x = zeros(N+1,1);

for i = 0:N
    x(i+1) = a + i*h;
end

%% Set up system

A = zeros(N-1,N-1);
z = zeros(N-1,1);

A(1,1) = -(2+(h^2)*q(x(2)));
A(1,2) =  1 - (h/2)*p(x(2));

for i = 2:N-2
    A(i,i-1) = 1 + (h/2)*p(x(i+1));
    A(i,i)= -2 - (h^2)* q(x(i+1));
    A(i,i+1) = 1 - (h/2)*p(x(i+1));
end

A(N-1,N-2) = 1 + (h/2)*p(x(N));
A(N-1,N-1)= -2 - (h^2)* q(x(N));

for i=1:N-1
    z(i) = (h^2)*r(x(i+1));
end

z(1) = z(1) - (1 + (h/2)*p(x(2)))*alpha;
z(N-1) = z(N-1) - (1 - (h/2)*p(x(N)))*beta;

%% Solve the system

w = A \ z;

%% Add on the boundary values to w

w = [alpha; w; beta];

%% Calculate errors
eer = zeros(1,N+1);
for i=1:N+1
    eer(i) = abs(w(i)-y(x(i)));
end
% %% Print
% fprintf('x \t\t w \t\t y \t\t error \n')
% for i=1:N+1
%     fprintf('%f \t %f \t %f \t %f \n', x(i),w(i),y(x(i)), eer(i))
% end
% 
% %% Plot the solution
% plot(x,eer,'-*')
     
%% Compute l1 norm

l1=0;
for i =2:N
    l1 = l1 + h*(eer(i));
end

%% Print
fprintf('L1 Norm of Error %f \t \n', l1)

%% Compute l2 norm

l2=0;
for i =2:N
    l2 = l2 + h*((eer(i))^2);
end

l2 = sqrt(l2);

%% Print
fprintf('L2 Norm of Error %f \t  \n', l2)



