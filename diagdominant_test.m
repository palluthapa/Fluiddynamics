function diagdominant_test(l,D,up)
% This function test whether the Matrix is Diagonally Dominant.

% Calculate Size of System.
N=length(D) ;

% Generate Tridiagonal Matrix.
A=zeros(N,N);

for i=2:N
    A(i,i-1)=l(i) ;
    A(i-1,i)=up(i-1) ;    
end

for i=1:N
    A(i,i)=D(i) ;
end

% Test for Diagonally Dominant.
if all((2*abs(diag(A))) >= sum(abs(A),2))
    disp('Matrix is Diagonally Dominant')
else
    disp('Matrix is Not Diagonally Dominant')
end
