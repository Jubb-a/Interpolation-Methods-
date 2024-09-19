#Math 501A Divided Differences for Project
import math
import sympy as sym
import numpy as np
import pandas as pd
import matplotlib.pyplot as mpl
from datetime import date

#import data from csv as numpy array
data = np.genfromtxt('C:\\Users\\ajubb\\OneDrive - UCI Health\\Documents\\personal\\Math 501A\\Project\\De-identified Pt Fall Data for Analysis.csv',skip_header=True,delimiter=',')


y_rate=data[0:28,5]
x_month= data[0:28,0]
t=sym.symbols('t')

def Newton(x_month,y_rate):
    #create interval with evenly spaced subdivisions
    x=x_month
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
        a[i,1]=y_rate[i]
        
    #calculate 1st order DD and insert into array
    for i in range(n-1):
        if(x[i+1]-x[i])!=0:
            a[i,2]=((y_rate[i+1]-y_rate[i])/(x[i+1]-x[i]))
            
    #calculating higher order DD and inserting result into array
    for j in range(2,n+1):
        count+=1
        for i in range(n-1-count):
            if (x[i+1-count]-x[i])!=0:
                a[i,j]=((a[i+1,j-1]-a[i,j-1])/(x[i+1+count]-x[i]))
    
    #build polynomial
    for i in range(n-1):
        prod=(t-x[i])
        k=k*prod
        polyprod.append(k)
    for i in range(1,n+1):
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
    return(redpn)
#End of function code

polynomial= Newton(x_month[0:10],y_rate[0:10])
polynomial2= Newton(x_month[9:28],y_rate[9:28])
xint=np.linspace(0,10,100)
yint=[]
xint2=np.linspace(9,27,100)
yint2=[]
for x in xint:
    y=sym.N(polynomial.subs(t,x))
    yint.append(y)
for x in xint2:
    y=sym.N(polynomial.subs(t,x))
    yint2.append(y)

print(polynomial)
print(polynomial2)


#raw data full
mpl.figure(1)
mpl.plot(data[:,0],data[:,5],'o')
mpl.vlines(10, 0, max(data[:,5]),colors='r')
mpl.text(11, max(data[:,5]), 'post CNL')
mpl.text(2, max(data[:,5]), 'pre CNL')
mpl.ylabel('Patient Fall Rate (%)')
mpl.xlabel('Month')

#raw data months 0-27 for analysis
mpl.figure(2)
mpl.plot(data[0:10,0],data[0:10,5],'o')
mpl.vlines(10, 0, max(data[:,5]),colors='r')
mpl.text(11, max(data[:,5]), 'post CNL')
mpl.text(2, max(data[:,5]), 'pre CNL')
mpl.ylabel('Patient Fall Rate (%)')
mpl.xlabel('Month')
mpl.plot(xint,yint)
#mpl.ylim(-10,10)

mpl.figure(3)
mpl.plot(data[9:28,0],data[9:28,5],'o')
mpl.vlines(10, 0, max(data[:,5]),colors='r')
mpl.text(11, max(data[:,5]), 'post CNL')
mpl.text(2, max(data[:,5]), 'pre CNL')
mpl.ylabel('Patient Fall Rate (%)')
mpl.xlabel('Month')
mpl.plot(xint2,yint2)

#this can be commented out to show full polynomial range
#mpl.ylim(-10,10)
