# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 10:14:42 2017

@author:
"""
""" Condensed 2 fluids, lower in Poiseuille flow
no pressure jump across interface
"""

import numpy as np
import matplotlib.pyplot as plt
import time

'Start timer'
tic=time.time()

'mesh size'
xmin=-3.; xmax=3.; ymin=0.; ymax=6.
nmax=100; mmax=100
""" a0 is incoming film thickness, Ka is the pressure parameter
in bottom fluid, mu and ro are for bottom fluid
kmax is number of nodes in bottom fluid
"""

kmax=11; ro=1.5**3; mu=1.5**2; a0=.2; nu=mu/ro
Ka=0.5
dx=(xmax-xmin)/(nmax-1.); dy=(ymax-ymin)/(mmax-1.)
x=np.linspace(xmin,xmax,nmax); y=np.linspace(ymin,ymax,mmax)
z=np.linspace(0,1.,kmax); dz=1./(kmax-1)




h=np.zeros(nmax); f=np.zeros(nmax)
''' 
h is wall, f is interface
initial conditions; s(mmax,nmax) is the streamfunction
w(,) is vorticity
'''
wp=np.zeros((nmax,mmax));up=np.zeros((nmax,mmax));vp=np.zeros((nmax,mmax))
sp=np.zeros((nmax,mmax))
wm=np.zeros((nmax,kmax));um=np.zeros((nmax,kmax));vm=np.zeros((nmax,kmax))
sm=np.zeros((nmax,kmax))
us=np.zeros(nmax)
dumdx=np.zeros(kmax); dupdx=np.zeros(mmax)
""" srtarting profiles """
lamb=1./mu+Ka*a0/mu
us0=a0/mu+Ka*a0*a0/2./mu
for m in range(mmax):
    wp[0,m]=1.; wp[1,m]=1.
    vp[0,m]=0; vp[1,m]=0.
    up[0,m]=y[m]+us0; up[1,m]=y[m]+us0
#    sp[0,m]=y[m]**2/2.;sp[1,m]=y[m]**2/2.

for k in range(kmax):
    wm[0,k]=lamb-Ka/mu*a0*z[k]; wm[1,k]=wm[0,k]
    vm[1,k]=0.; vm[1,k]=0.
    um[0,k]=lamb*a0*z[k]-Ka/2./mu*a0*a0*z[k]*z[k]
    um[1,k]=um[0,k]
#'roughness, h0 is height'
h0=-3.0


for n in range(nmax):
    h[n]=h0*np.exp(-x[n]**2)
    
f[0]=a0; f[1]=a0; us[0]=us0; us[1]=us0
    
ap=np.zeros(mmax);bp=np.zeros(mmax);cp=np.zeros(mmax);dp=np.zeros(mmax)
am=np.zeros(kmax);bm=np.zeros(kmax);cm=np.zeros(kmax);dm=np.zeros(kmax)
wpold=np.zeros(mmax);upold=np.zeros(mmax);vpold=np.zeros(mmax)
wmold=np.zeros(kmax);umold=np.zeros(kmax);vmold=np.zeros(kmax)
#sold=np.zeros(mmax);
pp=np.zeros(mmax); qp=np.zeros(mmax)
pm=np.zeros(kmax); qm=np.zeros(kmax)

wpnew=np.zeros(mmax); upnew=np.zeros(mmax); vpnew=np.zeros(mmax)
wmnew=np.zeros(kmax); umnew=np.zeros(kmax); vmnew=np.zeros(kmax)
#vnew=np.zeros(mmax); snew=np.zeros(mmax)

#auxiliary
wplus=np.zeros(mmax); wminus=np.zeros(kmax)

def thomas2fl(ap,bp,cp,dp,am,bm,cm,dm,r1,r2,wwall):
# computes w in both layers with w =1 at the top and  
# wm=wwall at the bottom; interfacial conditions require r1 and r2
#qp[mmax-2] is taken from boundary condition
    wplus[mmax-1]=1.
    pp[mmax-2]=0.; qp[mmax-2]=wplus[mmax-1] 
    for m in range(mmax-2,0,-1):
        den=ap[m]*pp[m]+bp[m]
        pp[m-1]=-cp[m]/den
        qp[m-1]=(dp[m]-ap[m]*qp[m])/den
    pm[1]=0.;qm[1]=wwall
    for k in range(1,kmax-1):
        den=cm[k]*pm[k]+bm[k]
        pm[k+1]=-am[k]/den
        qm[k+1]=(dm[k]-cm[k]*qm[k])/den
    """ interfacial conditions used here
    keep an eye on mu and r1, r2 """
    den=pp[0]-1.+r2*(pm[kmax-1]-1.)/mu
    wplus[0]=(r1-qp[0]-r2*qm[kmax-1])/den
#    print 'rs', r1, r2
#    wplus[0]=1.
    wminus[kmax-1]=wplus[0]/mu
#    wminus[kmax-1]=1./mu
    for m in range(mmax-2):
        wplus[m+1]=pp[m]*wplus[m]+qp[m]
    for k in range(kmax-1,1,-1):
        wminus[k-1]=pm[k]*wminus[k]+qm[k]
    wminus[0]=wwall
    
    return wplus,wminus
    
def coeff(HH,HH1,up,vp,um,vm):
#compute coefficients in momentum
    ap=1.-dy/2.*vp
    bp=-2.-1.5*dy*dy/dx*up
    cp=-ap+2.
    dp=dy*dy/2./dx*up*(-4.*wp[n-1,:]+wp[n-2,:])
    am=1.-HH1*dz/2.*vm
    bm=-2.-1.5*HH*dz*dz/dx*um
    cm=-am+2.
    dm=HH/2.*dz*dz/dx*um*(-4.*wm[n-1,:]+wm[n-2,:])  
    return am,bm,cm,dm,ap,bp,cp,dp

def solver(fold,wwallold):
#    global umnew,upnew
    H=fold-h[n]    
    HH=H*H/nu; HH1=HH/H
    am,bm,cm,dm,ap,bp,cp,dp=coeff(HH,HH1,upold,vpold,umold,vmold)

    """ PRESSURE JUMP SPECIFIED HERE
    The zero in r1 means no pressure jump"""
    r1=dy*(0+(1-ro)*usold*(3.*usold-4*us[n-1]+us[n-2])/2./dx+Ka)
    r2=mu*dy/dz/H
    wpnew,wmnew=thomas2fl(ap,bp,cp,dp,am,bm,cm,dm,r1,r2,wwallold)
    
  
    #    'compute new velocities'
    umnew[0]=0.
    for k in range(kmax-1):
        umnew[k+1]=umnew[k]+H*dz/2.*(wmnew[k]+wmnew[k+1])
    upnew[0]=umnew[kmax-1]
    for m in range(mmax-1):
        upnew[m+1]=upnew[m]+dy/2.*(wpnew[m]+wpnew[m+1])
    vmnew[0]=0.
    dumdx=(3.*H*umnew-4.*(f[n-1]-h[n-1])*um[n-1,:]+(f[n-2]-h[n-2])*um[n-2,:])/2./dx
    for k in range(kmax-1):
        vmnew[k+1]=vmnew[k]-dz/2*(dumdx[k]+dumdx[k+1])
    vpnew[0]=0.
    dupdx=(3.*upnew-4.*up[n-1,:]+up[n-2,:])/2./dx
    for m in range(mmax-1):
        vpnew[m+1]=vpnew[m]-dy/2.*(dupdx[m]+dupdx[m+1])

#tests 
    """ v at inerface"""
    T1=vmnew[kmax-1]
    """ boundary condition in upper"""
    T2=upnew[mmax-1]-(y[mmax-1]+us0)-(fold-a0)
    return T1,T2,wpnew,upnew,vpnew,wmnew,umnew,vmnew    

def Newton(fold,wwallold):
#Newton's iterations to satisfy interfacial conditions
    for itnewt in range(itermax):
        if itnewt==itermax-1:
            print 'Newton iter = MAX'
            break
        eps=0.0001
        T11,T21,wpnew,upnew,vpnew,wmnew,umnew,vmnew =solver(fold+eps,wwallold)
        T12,T22,wpnew,upnew,vpnew,wmnew,umnew,vmnew =solver(fold,wwallold+eps)
        T1,T2,wpnew,upnew,vpnew,wmnew,umnew,vmnew =solver(fold,wwallold)
        A1=(T11-T1)/eps;A3=(T21-T2)/eps;A2=(T12-T1)/eps;A4=(T22-T2)/eps
        det=A1*A4-A2*A3
        deltf=(T2*A2-T1*A4)/det;deltw=(T1*A3-T2*A1)/det
        T1,T2,wpnew,upnew,vpnew,wmnew,umnew,vmnew =solver(fold+deltf,wwallold+deltw)
        if(abs(deltf)+abs(deltw)) < 0.00001:
            break
        fold=fold+deltf; wwallold=wwallold+deltw
#    print 'iter=, T1,T2' ,itnewt,T1,T2
    return fold,wwallold
#new x
itermax=301
for n in range(2,nmax):
    wpold[:]=2*wp[n-1,:]-wp[n-2,:]; upold[:]=2*up[n-1,:]-up[n-2,:]
    vpold[:]=2*vp[n-1,:]-vp[n-2,:]
    wmold[:]=2*wm[n-1,:]-wm[n-2,:]; umold[:]=2*um[n-1,:]-um[n-2,:]
    vmold[:]=2*vm[n-1,:]-vm[n-2,:]
    usold=2*us[n-1]-us[n-2]
    fold=2*f[n-1]-f[n-2]
    wwallold=2*wm[n-1,0]-wm[n-2,0]
#start nonlinear iterations here
    for noniter in range (itermax):
        if noniter==itermax-1:
            print 'Nonlinear iter = MAX'
            break
#run Newton to find fold and wwallold
        fold,wwallold=Newton(fold,wwallold)
        T1,T2,wpnew,upnew,vpnew,wmnew,umnew,vmnew =solver(fold,wwallold)
        epsil=np.linalg.norm(wpnew-wpold)+np.linalg.norm(wmnew-wmold)
#        print wpnew
        if epsil < 0.0001:
            break
#        'relaxation'
        rel=-0.8
        wpold[:]=wpnew+rel*(wpnew-wpold); upold[:]=upnew+rel*(upnew-upold)
        vpold[:]=vpnew+rel*(vpnew-vpold)
        wmold[:]=wmnew+rel*(wmnew-wmold); umold[:]=umnew+rel*(umnew-umold)
        vmold[:]=vmnew+rel*(vmnew-vmold)
        usold=upold[0]
        
#    print 'end outer', noniter,'x=',x[n], 'tau=', wnew[0]
#    w[n,:]=wnew;u[n,:]=unew;v[n,:]=vnew;s[n,:]=snew;g[n]=gnew
    f[n]=fold; us[n]=usold
    wp[n,:]=wpnew; up[n,:]=upnew; vp[n,:]=vpnew
    wm[n,:]=wmnew; um[n,:]=umnew; vm[n,:]=vmnew
    print 'x=',x[n], 'noniter=', noniter,' wall shear=' ,wwallold

#plt.plot(x,wm[:,0])
#plt.plot(x,f-h)
fig=plt.figure()
fig.set_size_inches(10,6)
plt.plot(x,us)
plt.title("Shear Stress")
plt.xlabel("x")
plt.ylabel("omega")
plt.show()

#for plotting only
s=np.zeros((nmax,kmax+mmax-1));yplot=np.zeros((nmax,kmax+mmax-1))
xplot=np.zeros((nmax,kmax+mmax-1))
xplot=np.tile(x,(kmax+mmax-1,1))
xplot=xplot.transpose()
for n in range(nmax):
    for kk in range(kmax-1):
        yplot[n,kk]=h[n]+z[kk]*(f[n]-h[n])
    for kk in range(kmax-1,kmax+mmax-1):
        yplot[n,kk]=f[n]+y[kk-kmax+1]
#define streamfunction
for n in range(nmax):
    s[n,0]=0.
    for kk in range(0,kmax-2):
        s[n,kk+1]=s[n,kk]+(f[n]-h[n])*dz/2.*(um[n,kk]+um[n,kk+1])
    s[n,kmax-1]=a0*a0/2./mu
    for kk in range(kmax-1,kmax+mmax-2):
        s[n,kk+1]=s[n,kk]+dy/2.*(up[n,kk-kmax+1]+up[n,kk+1-kmax+1])
#start plotting
#xout, yout = np.meshgrid(x, yplot[:,0])
#levels = np.arange(-1.0,.1,0.001)
levels = np.arange(-1.0,2.,0.001)
fig=plt.figure()
plt.contour(xplot, yplot, s,levels=levels)
plt.title("Stream Lines")
plt.xlabel("x")
plt.ylabel("psi")
fig.set_size_inches(10,6)

plt.plot(x,h)
plt.plot(x,f)
plt.show()
  
toc=time.time()
cpu_time=toc-tic
print 'CPU time=',cpu_time, 'sec Elapsed'          
        
    