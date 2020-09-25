clear;clc
omega0=0.1;
j=3;
Nmax=100;
Mmax=100;
xmin=0;
xmax=6;
ymin=0;ymax=6;
dx=(xmax-xmin)/(Nmax-1);
dy=(ymax-ymin)/(Mmax-1);
x=xmin:dx:xmax;
y=ymin:dy:ymax;
unew=zeros(1,Mmax);
unew(Mmax)=0.5;
umax=unew(Mmax);
x=x(j);
newton_raphson(omega0,x,umax)


