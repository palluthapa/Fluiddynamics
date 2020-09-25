% single_fluid.m by Pallu Thapa

% This Code solves the system of Navier-Stokes/Prandtl
% equations using finite difference scheme.

% Clear Workspaces.
close all ; clear ; clc

% Start Timer
tic

% Define Variables.
Nmax=500 ;                 % Number of nodes in X-direction.
Mmax=500 ;                 % Number of nodes in Y-direction.
xmin=-3 ;                  % Minimum value of X.
xmax=3 ;                   % Maximum value of X.
ymin=0 ;                   % Minimum value of Y.
ymax=6 ;                   % Maximum value of Y.
dx=(xmax-xmin)/(Nmax-1) ;  % Step size in X-direction. 
dy=(ymax-ymin)/(Mmax-1) ;  % Step size in Y-direction.
err=1 ;                    % Error between Iterations.
tol=1e-5 ;                 % Convergence tolerance.
alpha=0.2 ;                % Relaxation Parameter.
k=0 ;                      % Iteration Counter.

% Define Spatial Array and function for the constriction.
x=xmin:dx:xmax ;
y=ymin:dy:ymax ;
f=zeros(1,Nmax);
f0=-3;
f=f0*exp(-x.^2);

% Define Global Matrices.
u=zeros(Mmax,Nmax) ;
v=zeros(Mmax,Nmax) ;
omega=zeros(Mmax,Nmax) ;
psi=zeros(Mmax,Nmax) ;
pprime=zeros(1,Nmax) ;
uold=zeros(1,Mmax);
vold=zeros(1,Mmax);
omegaold=zeros(1,Mmax);
a=zeros(1,Mmax-2);
b=zeros(1,Mmax-2);
c=zeros(1,Mmax-2);
d=zeros(1,Mmax-2);

% At n=1 or j=1 (first column).
for j=1
    for i=1:Mmax
        u(i,j)=y(i) ;
        omega(i,j)=1 ;
        psi(i,j)=((y(i))^2)/2 ;
    end
end

% At n=2 or j=2 (second column).
for j=2
    for i=1:Mmax
        u(i,j)=y(i) ;
        omega(i,j)=1 ;
        psi(i,j)=((y(i))^2)/2 ;
    end
end

% Loop through n= 3 to Nmax to or j= 3 to Nmax
for j=3:Nmax
    for i= 1:Mmax
        uold(i)=2*u(i,j-1)-u(i,j-2) ;
        vold(i)=2*v(i,j-1)-v(i,j-2) ; 
        omegaold(i)=2*omega(i,j-1)-omega(i,j-2) ;
    end
   
    err=1.;
    k=0;
    while err>tol
       k=k+1;
        for i=2:Mmax-1
            a(i)=1+(dy/2)*vold(i) ;
            b(i)=-2-uold(i)*(3/2)*((dy^2)/dx) ;
            c(i)=2-a(i) ;
            d(i)=((dy^2)/(2*dx))*uold(i)*(-4*omega(i,j-1)+omega(i,j-2));
        end
        
        % Define Lower, Main and Upper Diagonals and RHS.
        omega0=omegaold(1);
        l=zeros(1,Mmax-2) ;
        D=zeros(1,Mmax-2) ;
        up=zeros(1,Mmax-3) ;
        rhs=zeros(1,Mmax-2) ;
        l(2:Mmax-2)=a(3:Mmax-1) ;
        D(1:Mmax-2)=b(2:Mmax-1) ;
        up(1:Mmax-3)=c(2:Mmax-2) ;
        rhs(1)=d(2)-a(2)*omega0 ;
        rhs(2:Mmax-3)=d(3:Mmax-2) ;
        rhs(Mmax-2)=d(Mmax-1)-omegaold(Mmax)*c(Mmax-1);
        
        % Calculate omeganew using Thomas Algorithm and compute unew.
        % Match the Boundary condition by finding omega0 such that the
        % test function is very close to zero using Newton-Raphson.
        
        omeganew=zeros(1,Mmax) ;
        omeganew(1)=omega0 ;
        omeganew(Mmax)=omegaold(Mmax) ;
        diagdominant_test(l,D,up) ;
        omeganew(2:Mmax-1)=thomas(l,D,up,rhs) ;
        unew=zeros(1,Mmax) ;
        unew(1)=0 ;
        for i=2:Mmax
            unew(i)=unew(i-1)+0.5*dy*(omeganew(i)+omeganew(i-1)) ;  
        end
        tn=unew(Mmax)-y(Mmax)-f(j) ;
        tnp1=0 ;
        
        if abs(tn)>tol
        omeganew(1)=omega0+(1e-4) ;
        rhs(1)=d(2)-a(2)*(omega0+(1e-4)) ;
        diagdominant_test(l,D,up) ;
        omeganew(2:Mmax-1)=thomas(l,D,up,rhs) ;
        unew(1)=0 ;
          for i=2:Mmax
              unew(i)=unew(i-1)+0.5*dy*(omeganew(i)+omeganew(i-1)) ;  
          end
        tnp1=unew(Mmax)-y(Mmax)-f(j);
        domega0=(tnp1-tn)/(1e-4) ;
        omega0p1=omega0-(tn/domega0) ;
        else
        end

        while abs(tnp1)>tol
        omeganew(1)=omega0p1 ;
        rhs(1)=d(2)-a(2)*omega0p1 ;
        diagdominant_test(l,D,up) ;
        omeganew(2:Mmax-1)=thomas(l,D,up,rhs) ;
        unew(1)=0 ;
          for i=2:Mmax
              unew(i)=unew(i-1)+0.5*dy*(omeganew(i)+omeganew(i-1)) ; 
          end
        tnp1=unew(Mmax)-y(Mmax)-f(j) ;
        domega0=(tnp1-tn)/(omega0p1-omega0) ;
        omega0=omega0p1 ;
        omega0p1=omega0p1-(tnp1/domega0) ;
        tn=tnp1 ;
        end
        
        %calculate psinew and vnew.
        psinew=zeros(1,Mmax) ;
        vnew=zeros(1,Mmax) ;
        psinew(1)=0 ;
        vnew(1)=0 ;
        for i=2:Mmax
        psinew(i)=psinew(i-1)+0.5*dy*(unew(i)+unew(i-1)) ;
        vnew(i)=-((3*psinew(i)-4*psi(i,j-1)+psi(i,j-2))/(2*dx)) ;
        end
        
        % Update error and u.
        err=max(abs(unew-uold)) ;
        uold=uold+alpha*(unew-uold) ;
        vold=vold+alpha*(vnew-vold) ;
        omegaold=omegaold+alpha*(omeganew-omegaold) ;
    end
    
    % Fill the Global Matrix of u,v,omega and psi.
    u(:,j)=unew ;
    v(:,j)=vnew ;
    omega(:,j)=omeganew ;
    psi(:,j)=psinew ;
    k   
end
    
% Update Pressure.
for j=1:Nmax
pprime(j)=(omega(2,j)-omega(1,j))/dy ;
end

% Plot Shear Stress and Pressure Gradient along the Channel
plot(x,omega(1,:),'r');
xlabel('x','fontsize',16)
ylabel('Omega','fontsize',16)
title('Shear Stress','fontsize',16)
figure;
plot(x,pprime,'Color',[0.2,0.7,0.6]);
xlabel('x','fontsize',16)
ylabel('Pressure Gradient','fontsize',16)
title('Pressure Gradient along the Channel','fontsize',16)
figure;

% Contour Plot of Streamlines.
X=zeros(Mmax,Nmax);
Y=zeros(Mmax,Nmax);
for i=1:Mmax
    for j=1:Nmax
        X(i,j)=x(j);
        Y(i,j)=y(i)+f(j);
    end
end
levels=-1.:0.01:1.3;
contour(X,Y,psi,levels),hold on
plot(x,f,'r')
xlabel('x','fontsize',16)
ylabel('psi','fontsize',16)
title('Streamlines','fontsize',16)
hold off

cpu=toc;
fprintf('\n CPU Time =%1.4f seconds\n\n',cpu)