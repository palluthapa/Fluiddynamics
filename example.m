clear,clc
A=[1 2 0 0;2 1 3 0;0 2 1 3;0 0 2 1]
d=[2 1 5 6]';
x0=A\d
rhs=d';
l(2:4)=2;
D(1:4)=1;
up(1)=2;
up(2:3)=3;
x(1:4)=thomas(l,D,up,rhs)

