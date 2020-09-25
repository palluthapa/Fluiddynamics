function omega0=newton_raphson(omega0,x,umax)
difference=1;
tolerance=1e-5;

while difference>tolerance
    omega0_new=omega0-t(x,umax)/dt(x);
    difference=abs(omega0_new-omega0);
    omega0=omega0_new;
end

function test=t(x,umax) 
test=umax-6-exp(-x^2);

function dT=dt(x) 
dT=2*x*exp(-x^2);
