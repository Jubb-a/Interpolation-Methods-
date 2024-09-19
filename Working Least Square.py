#Math 501A Least Squares
import math
import sympy as sym
import numpy as np
import pandas as pd
import matplotlib.pyplot as mpl
from datetime import date
from matplotlib.pyplot import *

x=sym.symbols('x')
# pre data
x_pre = np.array([1., 2., 3., 4., 5., 6., 7., 8., 9.])
y_pre = np.array([0.22026432, 0.3960396 , 0.67264574, 0.46403712, 0.43196544,
                   0.        , 0.22123894, 0.21367521, 0.25      ])

x_post = np.array([10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20.,
                   21., 22., 23., 24., 25., 26., 27.])
y_post = np.array([0.40733198, 0.        , 0.21881838, 0.        , 0.41666667,
                0.21978022, 0.        , 0.21691974, 0.23255814, 0.22779043,
                0.22371365, 0.22988506, 0.41067762, 0.        , 0.        ,
                0.23148148, 0.23923445, 0.        ])
poly1=-0.00233849924069444*x**7 + 0.0729136605613889*x**6 - 0.929862571099444*x**5 + 6.24744184443472*x**4 - 23.6730195779315*x**3 + 49.9020966177539*x**2 - 53.0018366951783*x + 21.6048695379
poly2=-3.52788755858676e-11*x**16 + 1.03083749621601e-8*x**15 - 1.40397723571569e-6*x**14 + 0.000118299381386001*x**13 - 0.00690148648674111*x**12 + 0.295560503663906*x**11 - 9.61061384029875*x**10 + 242.015383808481*x**9 - 4769.45468646023*x**8 + 73795.6403184011*x**7 - 893395.669982893*x**6 + 8372961.21460506*x**5 - 59548243.1150479*x**4 + 310653432.274052*x**3 - 1121028760.90263*x**2 + 2499940162.40462*x - 2595273634.18736


def leastsquare(f,a,b,n):#f= function, a=lower, b= upper, n=number of steps
    x=sym.symbols('x')
    A=np.zeros((n+1,n+1))
    c=np.zeros((n+1,1))
    pn=0
    for i in range(n+1):
        for j in range(n+1):
            A[i,j]=sym.integrate(x**(i+j),(x,a,b))
            c[i,:]=sym.integrate(f*(x**i),(x,a,b))
            
    y=np.linalg.solve(A,c)
    vector= np.arange(n+1)[:,np.newaxis]
    exp=x**(vector[0:len(vector)])
    p= sum(y*exp)
    return p  

polynomial= leastsquare(poly1,1,9,3)
t=np.linspace(1,9,10)
y= polynomial[0].subs(x,t)

polynomial2= leastsquare(poly2,10,27,3)
t2=np.linspace(10,27,10)
y2= polynomial2[0].subs(x,t2)

mpl.figure(1)
mpl.plot(x_pre,y_pre,'o')
#mpl.vlines(10, 0, max(data[:,5]),colors='r')
#mpl.text(11, max(data[:,5]), 'post CNL')
#mpl.text(2, max(data[:,5]), 'pre CNL')
mpl.ylabel('Patient Fall Rate (%)')
mpl.xlabel('Month')
mpl.plot(t,y)
#mpl.ylim(-1,1)

mpl.figure(2)
mpl.plot(x_post,y_post,'o')
#mpl.vlines(10, 0, max(data[:,5]),colors='r')
#mpl.text(11, max(data[:,5]), 'post CNL')
#mpl.text(2, max(data[:,5]), 'pre CNL')
mpl.ylabel('Patient Fall Rate (%)')
mpl.xlabel('Month')
mpl.plot(t2,y2)
#mpl.ylim(-1,1)

