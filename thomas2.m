function x=thomas2(a,b,c,d)
% This function solves a tridiagonal system of linear equations
% using Thomas Algorithm.

% Calculate the size of system.
N=length(b);

% Perform Forward Sweep.
cprime(1)=c(1)/b(1);
dprime(1)=d(1)/b(1);

for i=1:(N-1)
    cprime(i+1)=c(i+1)/(b(i+1)-a(i+1)*cprime(i)) ;
    dprime(i+1)=(d(i+1)-a(i+1)*dprime(i))/(b(i+1)-a(i+1)*cprime(i)) ;
end

% Perform back Substitution to calculate x.
x(N)=dprime(N) ;
for i =N-1:-1:1
    x(i)=dprime(i)-cprime(i)*x(i+1) ;
end