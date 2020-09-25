# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 10:14:42 2017

@author: Sergei
"""
""" this is condensed bump problem with a slip on the wall
with speed uslip
"""
import numpy as np
import matplotlib.pyplot as plt

#choose slip on the wall
""" slip velocity must be constant"""
uslip=0
'mesh size'
xmin=-3.; xmax=3.; ymin=0.; ymax=6.
nmax=200; mmax=200
dx=(xmax-xmin)/(nmax-1.); dy=(ymax-ymin)/(mmax-1.)
x=np.linspace(xmin,xmax,nmax); y=np.linspace(ymin,ymax,mmax)
g=np.linspace(xmin,xmax,nmax); f=np.linspace(xmin,xmax,nmax)
''' initial conditions; s(mmax,nmax) is the streamfunction
w(,) is vorticity
'''
w=np.zeros((nmax,mmax));u=np.zeros((nmax,mmax));v=np.zeros((nmax,mmax))
s=np.zeros((nmax,mmax));
g[0]=0.; g[1]=0.
for m in range(mmax):
    w[0,m]=1.; w[1,m]=1.
    v[0,m]=0; v[1,m]=0.
    u[0,m]=y[m]+uslip; u[1,m]=y[m]+uslip
    s[0,m]=y[m]**2/2.+uslip*y[m];s[1,m]=y[m]**2/2.+uslip*y[m]

#'roughness, h0 is height'
#h0=-2.3
h0=-3.5
for n in range(nmax):
    f[n]=h0*np.exp(-x[n]**2)
    

a=np.zeros(mmax);b=np.zeros(mmax);c=np.zeros(mmax);d=np.zeros(mmax)
wold=np.zeros(mmax);uold=np.zeros(mmax);vold=np.zeros(mmax)
sold=np.zeros(mmax);
pp=np.zeros(mmax); qq=np.zeros(mmax)
wnew=np.zeros(mmax); unew=np.zeros(mmax)
vnew=np.zeros(mmax); snew=np.zeros(mmax)




def thomas(uold,vold,gxx):
    
#    'FUDGE ON U<0@'
#    for m in range (mmax):
#        if uold[m]<0.:
#            print m,uold[m]
#            uold[m]=0.001
#            vold[m]=0.
                
            
#        a[m]=1.-dy/2.*vold[m]
#        b[m]=-2.-1.5*dy*dy/dx*uold[m]
#        c[m]=-a[m]+2.
#        d[m]=dy*dy/2./dx*uold[m]*(-4.*w[n-1,m]+w[n-2,m])
    a=1.-dy/2.*vold
    b=-2.-1.5*dy*dy/dx*uold
    c=-a+2.
    d=dy*dy/2./dx*uold*(-4.*w[n-1,:]+w[n-2,:])
#for m in range(mmax):
#    a=1.-dy/2.*vold
#    b=-2.-1.5*dy*dy/dx*uold
#    c=-a+2.
#    d=dy*dy/2./dx*uold*(-4.*w[n-1,:]+w[n-2,:])
#    'thomas for new w'
    den=4.*c[1]+3*b[1]
    pp[2]=(c[1]-3.*a[1])/den
    qq[2]=(2.*dx*c[1]*gxx+3.*d[1])/den
    for m in range(2,mmax-1):
        den=b[m]+c[m]*pp[m]
        pp[m+1]=-a[m]/den
        qq[m+1]=(d[m]-c[m]*qq[m])/den
    ' w at outer edge --- specify here'       
    wnew[mmax-1]=1.
    for m in range(mmax-1,1,-1):
        wnew[m-1]=pp[m]*wnew[m]+qq[m]
    wnew[0]=(d[1]-a[1]*wnew[2]-b[1]*wnew[1])/c[1]

#    'compute new velocities and streamfunction'
    unew[0]=uslip; snew[0]=0.


    for m in range(mmax-1):
        unew[m+1]=unew[m]+dy/2.*(wnew[m]+wnew[m+1])
        snew[m+1]=snew[m]+dy/2.*(unew[m]+unew[m+1])
    vnew[0]=0.
    for m in range(1,mmax):
        vnew[m]=-(3*snew[m]-4*s[n-1,m]+s[n-2,m])/2./dx
# test is the boundary condition on u at infinity
    test=unew[mmax-1]-y[mmax-1]-f[n]-uslip
    return wnew, unew, vnew, snew, test

#'new x point'

for n in range(2,nmax):
#    'iterations for nonlinearity'    
#    'initial guesses'

    wold[:]=2*w[n-1,:]-w[n-2,:]; uold[:]=2*u[n-1,:]-u[n-2,:]
    vold[:]=2*v[n-1,:]-v[n-2,:]; sold[:]=2*s[n-1,:]-s[n-2,:]
    gnew=2*g[n-1]-g[n-2]
    itermax=1601
    for noniter in range(itermax):
#        gnew=g[n-1]
        
        if noniter==itermax-1:
                print 'max outer iter', 'x=', x[n], 'tau=', wnew[0]
#        'iterations for g'
        gold=gnew
        wnew, unew, vnew, snew, testold=thomas(uold,vold,gold)
        delg=0.001
        gnew=gold+delg
        wnew,unew,vnew,snew,testnew=thomas(uold,vold,gnew)
        
        for iter in range(1,itermax):
            if iter==itermax-1:
                print 'max inn iter', x[n]

            dtdg=(testnew-testold)/delg
            gold=gnew
            testold=testnew
            delg=-testold/dtdg
            gnew=gold+delg
            wnew,unew,vnew,snew,testnew=thomas(uold,vold,gnew)

#            print 'end inner iter=', iter, 'gnew=', gnew

            if np.absolute(delg) <0.00001:
                break
        
        eps=np.linalg.norm(wnew-wold)

        if eps < 0.000001:
            break
        'relaxation'
        rel=-0.8
        wold[:]=wnew+rel*(wnew-wold); uold[:]=unew+rel*(unew-uold)
        vold[:]=vnew+rel*(vnew-vold); sold[:]=snew+rel*(snew-sold)
    print 'end outer', noniter,'x=',x[n], 'tau=', wnew[0]
    w[n,:]=wnew;u[n,:]=unew;v[n,:]=vnew;s[n,:]=snew;g[n]=gnew

    
#    print x[n], wnew[0]


#print wnew[mmax-1],pp[mmax-1]
plt.plot (x, w[:,0])
plt.show() 
#plt.plot(y,u[101,:]-y)
#plt.show()

#for plotting only
splot=np.zeros((nmax,mmax));yplot=np.zeros((nmax,mmax))
xplot=np.zeros((nmax,mmax))
xplot=np.tile(x,(mmax,1))
xplot=xplot.transpose()
for n in range(nmax):
    for kk in range(mmax):
        yplot[n,kk]=f[n]+y[kk]
#define streamfunction
for n in range(nmax):
    
    for kk in range(mmax):
        splot[n,kk]=s[n,kk]
#start plotting
#xout, yout = np.meshgrid(x, yplot[:,0])
levels = np.arange(-1.0,2.,0.01)
plt.contour(xplot, yplot, splot,levels=levels)

plt.plot(x,f)
plt.show()

