#Math 501A Divided Differences for Project
import math
import sympy as sym
import numpy as np
import pandas as pd
import matplotlib.pyplot as mpl
from datetime import date

t=sym.symbols('t')
# pre data
x_pre = np.array([1., 2., 3., 4., 5., 6., 7., 8., 9.])
y_pre = np.array([0.22026432, 0.3960396 , 0.67264574, 0.46403712, 0.43196544,
                   0.        , 0.22123894, 0.21367521, 0.25      ])

## two options for post-data
## 10 - 27 -- the graph on desmos looked real wierd at the higher numbers
x_post = np.array([10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20.,
                   21., 22., 23., 24., 25., 26., 27.])
## 1-18 -- if we graph seperately from pre
x_post2 = np.array([ 1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12., 13.,
                    14., 15., 16., 17., 18.])
## post-rates the same either way
y_post = np.array([0.40733198, 0.        , 0.21881838, 0.        , 0.41666667,
                0.21978022, 0.        , 0.21691974, 0.23255814, 0.22779043,
                0.22371365, 0.22988506, 0.41067762, 0.        , 0.        ,
                0.23148148, 0.23923445, 0.        ])

def Newton(x_month,y_rate):
    #create interval with evenly spaced subdivisions
    x=x_month
    y= y_rate
    #build array of zeros for size defined by interval
    n=len(x)
    a=np.zeros((n,n+1))
    # define variables
    t=sym.symbols('t')
    count=0
    polyprod=[]
    c=[]
    k=1
   
    #filling array with x and y values
    for i in range(n):
        a[i,0]=x[i]
        a[i,1]=y[i]
   
    #calculate 1st order DD and insert into array
    for i in range(n-1):
        if(x[i+1]-x[i])!=0:
            a[i,2]=((y[i+1]-y[i])/(x[i+1]-x[i]))
   
    #calculating higher order DD and inserting result into array
    for j in range(3,n+1):
        count+=1
        for i in range(0,n-1-count):
                a[i,j]=((a[i+1,j-1]-a[i,j-1])/(x[i+1+count]-x[i]))
   
    #build polynomial
    for i in range(n-1):
        prod=(t-x[i])
        k=k*prod
        polyprod.append(k)
    for i in range(1,n):
        coeff=a[0,i]
        c.append(coeff)
    #I appended two zeros so that the lists for c and polyprod would be of equal size during multiplication loop
    c.append(0)
    c.append(0)
    pn=c[0]
    for i in range(n-1):
        poly=c[i+1]*polyprod[i]
        pn+=poly
    redpn=sym.simplify(pn)
    return redpn

poly_pre= Newton(x_pre,y_pre)
poly_post= Newton(x_post,y_post)
poly_post2= Newton(x_post2,y_post)
xint=np.linspace(0,10,100)
yint=[]
xint2=np.linspace(9,27,100)
yint2=[]
for x in xint:
    y=sym.N(poly_pre.subs(t,x))
    yint.append(y)
for x in xint2:
    y=sym.N(poly_post2.subs(t,x))
    yint2.append(y)

mpl.figure(1)
mpl.plot(x_pre,y_pre,'o')
mpl.vlines(10, 0, max(data[:,5]),colors='r')
mpl.text(11, max(data[:,5]), 'post CNL')
mpl.text(2, max(data[:,5]), 'pre CNL')
mpl.ylabel('Patient Fall Rate (%)')
mpl.xlabel('Month')
mpl.plot(xint,yint)
mpl.ylim(-1,1)

mpl.figure(2)
mpl.plot(x_post,y_post,'o')
#mpl.vlines(10, 0, max(data[:,5]),colors='r')
#mpl.text(11, max(data[:,5]), 'post CNL')
#mpl.text(2, max(data[:,5]), 'pre CNL')
mpl.ylabel('Patient Fall Rate (%)')
mpl.xlabel('Month')
mpl.plot(xint2,yint2)
mpl.ylim(-1,1)
