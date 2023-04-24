%% Input 

h=0.05;

bL = @ (y) sin(2*pi*y);
bR = @ (y) 0;
bT = @ (x) 0;
bB = @ (x) 0;

%% Initialization

N = 1/h;
xh = 0:h:1;
yh = xh;
M = N - 1;

w = zeros(M);
%v = zeros(M*M,1);
b = zeros(M*M,1);

D = zeros(M);
I = eye(M);

%% Construction of D

for i = 1:M
    if i == 1
        D(i,i)=-4;
        D(i,i+1)=1;
    elseif i == M
        D(i,i-1) = 1;
        D(i,i) = -4;
    else
        D(i,i-1) = 1;
        D(i,i) = -4;
        D(i,i+1) = 1;
    end
end

%% Construction of A

for i = 1:M
    if i == 1
        A((i-1)*M + 1: (i-1)*M + M, (i-1)*M + 1: (i-1)*M + M) = D;
        A((i-1)*M + 1: (i-1)*M + M, i*M + 1: i*M + M) = I;
    elseif i == M
        A((i-1)*M + 1: (i-1)*M + M, (i-1)*M + 1: (i-1)*M + M) = D;
        A((i-1)*M + 1: (i-1)*M + M, (i-2)*M + 1: (i-2)*M + M) = I;
    else
        A((i-1)*M + 1: (i-1)*M + M, (i-1)*M + 1: (i-1)*M + M) = D;
        A((i-1)*M + 1: (i-1)*M + M, i*M + 1: i*M + M) = I;
        A((i-1)*M + 1: (i-1)*M + M, (i-2)*M + 1: (i-2)*M + M) = I;
    end
end

%% Construction of b

for i = 1:M
    if i == 1
        b((i-1)*M + 1, 1) = -bL(yh(2)) - bB(xh(2));
        b((i-1)*M + 2, (i-1)*M + M - 1, 1) = - bB(xh(3:M));
        b((i-1)*M + M, 1) = -bR(yh(2)) -bB(xh(M+1));
    elseif i == M
        b((i-1)*M + 1, 1) = -bL(yh(M+1)) - bB(xh(2));
        b((i-1)*M + 2, (i-1)*M + M - 1, 1) = - bB(xh(3:M));
        b((i-1)*M + M, 1) = -bR(yh(M+1)) -bB(xh(M+1));
    else
        b((i-1)*M + 1, 1) = -bL(yh(i+1));
        b((i-1)*M + M, 1) = -bR(yh(i+1)); 
    end
end   

v=inv(A)*b;

%% Vector v to matrix w

for j = 1:M
    for i = 1:M
        w(i,j) = v((j-1)*M +i);
    end
end

surf(xh(2:N),xh(2:N)',w')
