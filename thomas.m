function omeganew=thomas(l,D,up,rhs)
% thomas.m by Pallu Thapa

% This function solves a tridiagonal system of linear equations
% using Thomas Algorithm. Here, The inputs l is the lower diagonal,
% D is the main diagonal and up is the Upper diagonal of the 
% tridiagonal matrix and rhs is the known column vector of the 
% matrix system Ax=b.

% Calculates the size of the system.
N=length(D);

% Perform Forward Sweep.
up(1)=up(1)/D(1);
rhs(1)=rhs(1)/D(1);
for i=2:N-1;
    up(i)=up(i)/(D(i)-l(i)*up(i-1));
    rhs(i)=(rhs(i)-l(i)*rhs(i-1))/(D(i)-l(i)*up(i-1));
end
rhs(N)=(rhs(N)-l(N)*rhs(N-1))/(D(N)-l(N)*up(N-1));

% Perform back Substitution to calculate x.
omeganew(N)=rhs(N);
for i =N-1:-1:1
    omeganew(i)=rhs(i)-up(i)*omeganew(i+1);
end