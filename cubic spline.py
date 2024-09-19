# Math501A − Project − Cubic Spline − Ren Padgett and Ashley Jubb
import numpy as np
from sympy import symbols , lambdify
from numpy import zeros , arange , diagflat
import matplotlib.pyplot as plt
import pandas as pd
# df=pd.readexcel(r"C:\Users\ajubb\OneDrive - UCI Health\Documents\personal\Math 501A\De-identified Pt Fall Data for Analysis.xlsx" )
# dfarray=pd.DataFrame.tonumpy(df)
# X = dfarray[: , 0]
# y = dfarray[: , 1]
X=np.array([1,2,5,8,9])
y=np.array([0.22026432, 0.3960396, 0.43196544, 0.21367521, 0.25])
t = X
n = len(t)
def zvals (t) :
    h=zeros(n-1)
    for i in arange(n-1):
        h[i] = t[i+1] - t[i]
    # u −−> diagonal of the matrix
    u = zeros (n-1)
    for i in arange(1 , n-1):
        u[i] = 2 * (h[i]+ h[i-1]) # left u[0] blank so indeces match h and b
    # calculate b values
    b = zeros(n-1)
    for i in arange(n-1):
        b[i] = 6*(y[i +1]-y[i])/h[i]
    # calc v values −−> right hand side of matrix Az=v
    v = zeros (n-1)
    for i in arange(1 ,n -1):
        v[i] = b[i]-b[i -1]
    # drop zero for u & v −− they should be size n−2
    u = u[ 1 : ]
    v = v[ 1 : ]
    # create matrix
    A = diagflat(u)
    for i in arange(1, n-2):
        A[i, i -1] = h[i]
        A[i -1, i] = h[i]
    z = zeros(n)
    for i in arange(n -2 ,0 , -1):
        z[i] = (v[i-1] - h[i] *z[i +1])/u[i-1]
    return z

z = zvals(t)
x = symbols ('x')

def Si (z,y,t,i) : # create spline polynomials
    hi = t[i +1] - t[i] #only need one h−val
    S= z[i]*(t[i+1]-x)**3/(6*hi) + z[i +1]*(x-t[i])**3/(6*hi)+ ((y[i+1])/hi-(z[i +1]*hi)/6)*(x-t[i]) + ((y[i])/hi-(z[i]*hi)/6)*(t[i+1]-x)
    return S

def Seval (xi,i) : #evaluate Si at a given x
    S_i = Si(z,y,t,i)
    s = lambdify(x, S_i)
    return (xi)
# all data points
xval= np.linspace(1,27, num=200)
nums = len(xval)
yval= zeros(nums)
for j in arange(n -1):
    for i in arange(nums):
        if xval[i] >= t[j] and xval[i]< t[j+1]:
            yval[i]= Seval(xval[i],j)
        if xval[i]>=t[-2] :
            yval[i] = Seval(xval[i], n-2)
## pre−CNL data points
xpre= np.linspace (1,10, num=200)
ypre = zeros (len(xpre))
for j in arange ( 9 ) :
    for i in range (len(xpre)) :
        if xpre[i]>= t[j] and xpre[i] <t[j+1]:
            ypre[i] = Seval(xpre[i], j)
        if xpre[i]>= t[-2]:
            ypre[i] = Seval(xpre[i], 7)