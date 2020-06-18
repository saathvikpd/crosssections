#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:00:46 2020

@author: saathvikdirisala
"""

from warnings import filterwarnings
filterwarnings("ignore")
import matplotlib.pyplot as plt
import pandas as pd
from numpy import sqrt
from math import isnan
from math import atan
from math import sin
from math import cos
import random
import numpy as np

'''
Cross-sectional area distribution generator for three-dimensional models formed by
drawing two polynomials of choice on the x-y plane and projecting two parabolas
from one curve to the other along the y-z plane, on both sides of the x-y plane
'''

data = [1]

Px = np.array([-1,0,-0.05]) #NOTE: Start with the intercept of the polynomial
'''Array that represents a polynomial that bounds the shape of the model
from below, on the x-y plane'''
Qx = np.array([1,0,-0.08]) #NOTE: Start with the intercept of the polynomial
'''Array that represents a polynomial that bounds the shape of the model
from above, on the x-y plane'''

'''This function finds the product of two polynomials'''
def prpoly(P,Q):
    a = [1]
    c =[1]
    beta = 0
    for i in range(len(P)):
       for j in range(len(Q)):
           c[beta] = i+j
           c.append(1)
           if(beta != 0):
              d = c.index(c[beta])
              if(d != beta):
                   a[d]=a[d]+(P[i]*Q[j])
              else:
                   a[beta] = P[i]*Q[j]
                   a.append(1)
                   beta=beta+1
           else:
              a[beta] = P[i]*Q[j]
              a.append(1) 
              beta=beta+1
    a = a[:(len(a)-1)]
    a = np.array(a)
    return a

prPP = prpoly(Px,Px)
prQQ = prpoly(Qx,Qx)

'''This function finds the sum of two polynomials'''
def smpoly(pP,pQ):
    gy = len(pQ)-len(pP)
    if(gy>0):
        pP = list(pP)
        pP.extend(abs(gy)*[0])
        pP = np.array(pP)
    else:
        pQ = list(pQ)
        pQ.extend(abs(gy)*[0])
        pQ = np.array(pQ)
    return(pQ+pP)

for treq in range(1):
    '''
    Equation for the plane cutting the cross-section:
        c1*x + c2*y + c3*z = 0
    '''
    '''Random generation of coefficients from a normal distribution'''
    c1 = round(random.normalvariate(0,5),ndigits = 5) 
    c2 = round(random.normalvariate(0,5),ndigits = 5)
    '''c3 is held constant because the randomization of three coefficients
    would be redundant'''
    c3 = -1


    '''Expanding equations to be presented as polynomial arrays:
        -Check poly0, poly1 and poly2 below'''
#=======================================================================================
    
    art = prpoly(Px,Qx)
    prPQ = (c2**2)*(art)
    pPx = (c1*c2)*Px
    pQx = (c1*c2)*Qx
    '''pPx = list(pPx)
    pQx = list(pQx)
    pPx.append(0)
    pQx.append(0)
    pPx = np.roll(pPx,1)
    pQx = np.roll(pQx,1)'''
    pPx = prpoly(np.array([0,1]),pPx)
    pQx = prpoly(np.array([0,1]),pQx) 
    fx = smpoly(pPx,pQx)
    hx = smpoly(fx,prPQ)
    hx = smpoly(hx,np.array([0,0,c1**2]))
    hx = hx[::-1]
    
#=======================================================================================
    
    ix = smpoly(prQQ,prPP)
    jx = smpoly(-2*art,ix)
    lx = smpoly(2*c2*Px,2*c2*Qx)
    mx = smpoly(lx,jx)
    nx = smpoly(mx,np.array([c2**2,4*c1]))
    nx = nx[::-1]
    
#=======================================================================================
    
    ly = smpoly(-2*c2*Px,-2*c2*Qx)
    my = smpoly(ly,jx)
    ny = smpoly(my,np.array([c2**2,-4*c1]))
    ny = ny[::-1]
    
    '''Initializing the expanding lists used in these iterations'''
    crit0 = [1]
    crit1 = [1]
    crit2 = [1]
    redy = [1]
    bluey = [1]
    orangey = [1]
    purpley = [1]
    redx2 = [1]
    redy2 = [1]
    bluex2 = [1]
    bluey2 = [1]
    orangex2 = [1]
    orangey2 = [1]
    purplex2 = [1]
    purpley2 = [1]
    redz = [1]
    bluez = [1]
    orangez = [1]
    purplez = [1]
    
    
    '''The two base polynomials specified above (Px,Qx) in functions form'''
    def base1(x):
        d = 0
        for u in range(len(Qx)):
            d = d + Qx[u]*(x**u)
        return d
    def base2(x):
        d = 0
        for u in range(len(Px)):
            d = d + Px[u]*(x**u)
        return d
    
    '''The functions that define the three-dimensional model'''
    def redorange(x,y):
        d = -(y-(base1(x)))*(y-(base2(x)))
        return d
    def bluepurple(x,y):
        d = (y-(base1(x)))*(y-(base2(x)))
        return d
    
    '''The functions obtained by setting the equation for the three-dimensional model and
    the equation for the sectioning plane in terms of x and y equal to one another 
    (They were color-coded on the graphing calculator according to their namesake)'''
    def orange(x):
        d = smpoly(Px,Qx)
        d[0]=d[0]-c2
        res = 0
        bes = 0
        psny = ny[::-1]
        for u in range(len(d)):
            res = res + d[u]*(x**u)
        for t in range(len(psny)):
            bes = bes + psny[t]*(x**t)
        end = (res + sqrt(bes))/2
        return end
    def red(x):
        d = smpoly(Px,Qx)
        d[0]=d[0]-c2
        res = 0
        bes = 0
        psny = ny[::-1]
        for u in range(len(d)):
            res = res + d[u]*(x**u)
        for t in range(len(psny)):
            bes = bes + psny[t]*(x**t)
        end = (res - sqrt(bes))/2
        return end
    def purple(x): 
        d = smpoly(Px,Qx)
        d[0]=d[0]+c2
        res = 0
        bes = 0
        psny = nx[::-1]
        for u in range(len(d)):
            res = res + d[u]*(x**u)
        for t in range(len(psny)):
            bes = bes + psny[t]*(x**t)
        end = (res + sqrt(bes))/2
        return end
    def blue(x): 
        d = smpoly(Px,Qx)
        d[0]=d[0]+c2
        res = 0
        bes = 0
        psny = nx[::-1]
        for u in range(len(d)):
            res = res + d[u]*(x**u)
        for t in range(len(psny)):
            bes = bes + psny[t]*(x**t)
        end = (res - sqrt(bes))/2
        return end
    
    '''basecrit stores the intersection points of the two polynomials on the x-y plane'''
    polybase = np.poly1d(smpoly(Px,-1*Qx)[::-1]).r
    pb = np.array(polybase, dtype=float)
    basecrit = [1]
    for w in range(len(polybase)):
        if(polybase[w]==pb[w]):
            basecrit[len(basecrit)-1]=pb[w]
            basecrit.append(1)
    basecrit = basecrit[:(len(basecrit)-1)]
    basecrit = sorted(polybase)
    basecrit = [min(basecrit),max(basecrit)]
    '''We need only the smallest and the largest values to bound the other critical points'''
    
    '''hx is the polynomial that has its zeroes located at the points of intersection
    of the red or orange polynomial and the blue or purple polynomial'''
    poly0 = np.poly1d(hx).r
    p = np.array(poly0, dtype=float)
    a1=0
    for m in range(len(p)):
        if(isnan(round(red(p[m]), ndigits = 5))==False and isnan(round(orange(p[m]), ndigits = 5))==False and isnan(round(blue(p[m]), ndigits = 5))==False and isnan(round(purple(p[m]), ndigits = 5))==False):
            if (round(poly0[m]-p[m])==0) and (p[m]>basecrit[0]) and (p[m]<basecrit[1]): 
                crit0[a1] = p[m]
                a1+=1
                crit0.append(1)
    crit0 = crit0[:(len(crit0)-1)]
    crit0 = sorted(crit0)
    
    '''ny is the polynomial that has its zeroes located at the points of intersection
    of the red and orange polynomials'''
    acc = 3
    poly1 = np.poly1d(ny).r
    q = np.array(poly1, dtype=float)
    b1=0
    for n in range(len(q)):
        if (round(poly1[n]-q[n])==0):
            while ((isnan(red(round(q[n], ndigits = acc)))==True or isnan(orange(round(q[n], ndigits = acc)))==True) and acc<50):
                acc = acc + 1
             #print(acc)
            if(acc>50):
                break
            if(isnan(red(round(q[n], ndigits = acc)))==False and isnan(orange(round(q[n], ndigits = acc)))==False):
                if q[n]>basecrit[0] and q[n]<basecrit[1]: 
                    crit1[b1] = q[n]
                    b1+=1
                    crit1.append(1)
    crit1 = crit1[:(len(crit1)-1)]
    crit1 = sorted(crit1)
    
    '''nx is the polynomial that has its zeroes located at the points of intersection
    of the blue and purple polynomials'''
    acc1 = 3
    poly2 = np.poly1d(nx).r
    r = np.array(poly2, dtype = float)
    ce=0
    for o in range(len(r)):
        if (round(poly2[o]-r[o])==0):
            while ((isnan(blue(round(r[o], ndigits = acc1)))==True or isnan(purple(round(r[o], ndigits = acc1)))==True) and acc1<50):
                acc1 = acc1 + 1
             
            if (acc1>50):
                break
          
            if (isnan(blue(round(r[o], ndigits = acc1)))==False and isnan(purple(round(r[o], ndigits = acc1)))==False):
                if r[o]>basecrit[0] and r[o]<basecrit[1]: 
                    crit2[ce] = r[o]
                    ce+=1
                    crit2.append(1)
    
    crit2 = crit2[:(len(crit2)-1)]
    crit2 = sorted(crit2)
    
    '''The critical values found above are filtered out based on whether they lie within the bounds of the basecrit points'''
    
    '''The critical values are now compiled and stored in crits after they are sorted'''
    crits = [1]
    v = 0
    w = crit0
    w.extend(crit1)
    w.extend(crit2)
    for ab in range(len(w)):
         crits[v] = w[ab]
         v+=1
         crits.append(1)
       
    crits = crits[:(len(crits)-1)]
    crits = sorted(crits)
    
    '''======================================================================================
    This part builds the transformation matrix that will bring the cross section to be parallel 
    with the x-y plane for easier integration
    '''
    
    cs = (c3/(sqrt(c1**2 + c2**2 + c3**2)))
    sn = (sqrt((c1**2 + c2**2)/(c1**2 + c2**2 + c3**2)))
    u1 = (c2/(sqrt(c1**2 + c2**2 + c3**2)))
    u2 = (-c1/(sqrt(c1**2 + c2**2 + c3**2)))
    
    a = cs + (u1**2)*(1-cs)
    b = u1*u2*(1-cs)
    c = u2*sn
    d = u1*u2*(1-cs)
    e = cs + (u2**2)*(1-cs)
    f = -u1*sn
    g = -u2*sn
    h = u1*sn
    i = cs
    
    '''
    The transformation elements would be arranged in a matrix as such:
        [a,  b,  c]
        [d,  e,  f]
        [g,  h,  i]
    '''
    
    '''Transformation function that returns only the x and y values because z will be 0'''
    def trans(x,y,z):
        x1 = (round(x*a + y*b + z*c,ndigits = 20))
        y1 = (round(x*d + y*e + z*f,ndigits = 20))
        x2 = x1*cos(-atan(-c1/c2)) - y1*sin(-atan(-c1/c2))
        y2 = x1*sin(-atan(-c1/c2)) + y1*cos(-atan(-c1/c2))
        return x2, y2
    '''====================================================================================='''
    
    '''xspan defines the domain within which we want to integrate the 4 functions, while
    xinc increments through the domain'''
    xspan = crits[len(crits)-1] - crits[0]
    ks = 0
    w2 = 0
    mm = 0
    vk = 0
    areablue = 0
    areaorange = 0
    areared = 0
    areapurple = 0
    xinc = [1]
    
    '''This integrator loop runs through each x-value in the xinc list, transforms 
    the corresponding x, y and z values, and finds the area enclosed by the transformed 
    x and y values for each of the 4 functions'''
    prt = 50 #No. of partitions of xspan | No. of increments
    for to in range(prt):
       xinc[to] = crits[0]+((xspan/prt)*(to+1))
       
       bluey[to] = blue(xinc[to])
       bluez[to] = bluepurple(xinc[to],bluey[to])

       redy[to] = red(xinc[to])
       redz[to] = redorange(xinc[to],redy[to])
       
       orangey[to] = orange(xinc[to])
       orangez[to] = redorange(xinc[to],orangey[to])
           
       purpley[to] = purple(xinc[to])
       purplez[to] = bluepurple(xinc[to],purpley[to])
       
       xinc.append(1)
       bluey.append(1)
       redy.append(1)
       orangey.append(1)
       purpley.append(1)
       bluez.append(1)
       redz.append(1)
       orangez.append(1)
       purplez.append(1)
       
       if(isnan(bluey[to])==False):
              if(bluey[to]<=base1(xinc[to]) and bluey[to]>=base2(xinc[to])):
                 bluex2[ks] = trans(xinc[to], bluey[to], bluez[to])[0]
                 bluey2[ks] = trans(xinc[to], bluey[to], bluez[to])[1]
                 if(ks!=0):
                  areablue = abs(bluex2[ks]-bluex2[ks-1])*abs(bluey2[ks]) + areablue
                 ks = ks + 1
                 bluex2.append(1)
                 bluey2.append(1)
         
       if(isnan(redy[to])==False):
              if(redy[to]<=base1(xinc[to]) and redy[to]>=base2(xinc[to])):
                 redx2[w2] = trans(xinc[to], redy[to], redz[to])[0]
                 redy2[w2] = trans(xinc[to], redy[to], redz[to])[1]
                 if(w2!=0):
                  areared = abs(redx2[w2]-redx2[w2-1])*abs(redy2[w2]) + areared
                 w2 = w2 + 1
                 redx2.append(1)
                 redy2.append(1)
               
       if(isnan(orangey[to])==False):
              if(orangey[to]<=base1(xinc[to]) and orangey[to]>=base2(xinc[to])):
                 orangex2[mm] = trans(xinc[to], orangey[to], orangez[to])[0]
                 orangey2[mm] = trans(xinc[to], orangey[to], orangez[to])[1]
                 if(mm!=0):
                     areaorange = abs(orangex2[mm]-orangex2[mm-1])*abs(orangey2[mm]) + areaorange
                 mm = mm + 1
                 orangex2.append(1)
                 orangey2.append(1)
          
       if(isnan(purpley[to])==False):
              if(purpley[to]<=base1(xinc[to]) and purpley[to]>=base2(xinc[to])):
                 purplex2[vk] = trans(xinc[to], purpley[to], purplez[to])[0]
                 purpley2[vk] = trans(xinc[to], purpley[to], purplez[to])[1]
                 if(vk!=0):
                    areapurple = abs(purplex2[vk]-purplex2[vk-1])*abs(purpley2[vk]) + areapurple 
                 vk = vk + 1
                 purplex2.append(1)
                 purpley2.append(1)
       
       totalarea = areared + areapurple + areablue + areaorange

    redx2 = redx2[:len(redx2)-1]
    orangex2 = orangex2[:len(orangex2)-1]
    bluex2 = bluex2[:len(bluex2)-1]
    purplex2 = purplex2[:len(purplex2)-1]
    redy2 = redy2[:len(redy2)-1]
    orangey2 = orangey2[:len(orangey2)-1]
    bluey2 = bluey2[:len(bluey2)-1]
    purpley2 = purpley2[:len(purpley2)-1]
    
    '''Plots the shape of the cross section:
        -Activate only when working with a single iteration'''
    shapex = list(redx2)
    shapex.extend(list(orangex2[::-1]))
    shapex.extend(list(purplex2[::-1]))
    shapex.extend(list(bluex2))
    shapey = list(redy2)
    shapey.extend(list(orangey2[::-1]))
    shapey.extend(list(purpley2[::-1]))
    shapey.extend(list(bluey2))
    shape = pd.DataFrame({'x':shapex,'y':shapey})
    shape.sort_values(by=['x','y'])
    plt.plot(shape['x'],shape['y'])
    
    totalarea=round(totalarea, ndigits = 5)
    #print(c1,c2,totalarea)
    data[treq] = totalarea
    data.append(1)


'''This piece of code builds the histogram for the sample of cross-sectional areas 
found for the three-dimensional model'''    
data = data[:(len(data)-1)]
'''In order to view the cross-sectional area distribution for the model, remove 
the hashtag on the line below, add a hashtag to line 449, and change the iteration 
number on line 77 from 1 to whatever sample size you prefer'''
#plt.hist(data, bins = round(len(data)/10), color = 'red')




'''
IGNORE:
===========================================================================================
Remnants of older pieces of code that have been replaced by newer, more efficient ones:


def rnd(x):
    for i in range(len(x)):
        x[i] = int(round(x[i], ndigits=-2))
    return(x)
#data = rnd(data)
data = tuple(data)
def freq(x):
    hs = {}
    for j in data:
        hs[j] = hs.get(j,0) + 1
    return(hs)
count = freq(data)
unique = tuple(set(data))
ct = [1]
for k in range(len(unique)-1):
    ct[k] = count.get(unique[k])
    ct.append(1)
ct = ct[:len(ct)]'''
'''bluex1[0]= (round(xinc[to]*a + bluey[to]*b + bluepurple(xinc[to],bluey[to])*c, ndigits = 20))
bluey1[0]= (round(xinc[to]*d + bluey[to]*e + bluepurple(xinc[to],bluey[to])*f, ndigits = 20))
bluex2[0] = bluex1[0]*cos(-atan(-c1/c2)) - bluey1[0]*sin(-atan(-c1/c2))
bluey2[0] = bluex1[0]*sin(-atan(-c1/c2)) + bluey1[0]*cos(-atan(-c1/c2))
ks = ks + 1
bluex1.append(1)
bluey1.append(1)
bluex2.append(1)
bluey2.append(1)
elif(bluey[to]<=base1(xinc[to]) and bluey[to]>=base2(xinc[to]) and ks!=0):'''
'''bluex1[ks]= (round(xinc[to]*a + bluey[to]*b + bluepurple(xinc[to],bluey[to])*c,ndigits = 20))
bluey1[ks]= (round(xinc[to]*d + bluey[to]*e + bluepurple(xinc[to],bluey[to])*f,ndigits = 20))
bluex2[ks] = bluex1[ks]*cos(-atan(-c1/c2)) - bluey1[ks]*sin(-atan(-c1/c2))
bluey2[ks] = bluex1[ks]*sin(-atan(-c1/c2)) + bluey1[ks]*cos(-atan(-c1/c2))
                 '''

'''redx1[0]= (round(xinc[to]*a + redy[to]*b + redorange(xinc[to],redy[to])*c,ndigits = 20))
redy1[0]= (round(xinc[to]*d + redy[to]*e + redorange(xinc[to],redy[to])*f,ndigits = 20))
redx2[0] = redx1[0]*cos(-atan(-c1/c2)) - redy1[0]*sin(-atan(-c1/c2))
redy2[0] = redx1[0]*sin(-atan(-c1/c2)) + redy1[0]*cos(-atan(-c1/c2))
w2 = w2 + 1
redx1.append(1)
redy1.append(1)
redx2.append(1)
redy2.append(1)
elif(redy[to]<=base1(xinc[to]) and redy[to]>=base2(xinc[to]) and w2!=0):
             redx1[w2]= (round(xinc[to]*a + redy[to]*b + redorange(xinc[to],redy[to])*c,ndigits = 20))
             redy1[w2]= (round(xinc[to]*d + redy[to]*e + redorange(xinc[to],redy[to])*f,ndigits = 20))
             redx2[w2] = redx1[w2]*cos(-atan(-c1/c2)) - redy1[w2]*sin(-atan(-c1/c2))
             redy2[w2] = redx1[w2]*sin(-atan(-c1/c2)) + redy1[w2]*cos(-atan(-c1/c2))'''


'''orangex1[0]= (round(xinc[to]*a + orangey[to]*b + redorange(xinc[to],orangey[to])*c,ndigits = 20))
orangey1[0]= (round(xinc[to]*d + orangey[to]*e + redorange(xinc[to],orangey[to])*f,ndigits = 20))
             orangex2[0] = orangex1[0]*cos(-atan(-c1/c2)) - orangey1[0]*sin(-atan(-c1/c2))
             orangey2[0] = orangex1[0]*sin(-atan(-c1/c2)) + orangey1[0]*cos(-atan(-c1/c2))
             mm = mm + 1
             orangex1.append(1)
             orangey1.append(1)
             orangex2.append(1)
             orangey2.append(1)
          elif(orangey[to]<=base1(xinc[to]) and orangey[to]>=base2(xinc[to]) and mm!=0):
             
             orangex1[mm]= (round(xinc[to]*a + orangey[to]*b + redorange(xinc[to],orangey[to])*c,ndigits = 20))
             orangey1[mm]= (round(xinc[to]*d + orangey[to]*e + redorange(xinc[to],orangey[to])*f,ndigits = 20))
             orangex2[mm] = orangex1[mm]*cos(-atan(-c1/c2)) - orangey1[mm]*sin(-atan(-c1/c2))
             orangey2[mm] = orangex1[mm]*sin(-atan(-c1/c2)) + orangey1[mm]*cos(-atan(-c1/c2))
'''
''' purplex1[0]= (round(xinc[to]*a + purpley[to]*b + bluepurple(xinc[to],purpley[to])*c,ndigits = 20))
             purpley1[0]= (round(xinc[to]*d + purpley[to]*e + bluepurple(xinc[to],purpley[to])*f,ndigits = 20))
             purplex2[0] = purplex1[0]*cos(-atan(-c1/c2)) - purpley1[0]*sin(-atan(-c1/c2))
             purpley2[0] = purplex1[0]*sin(-atan(-c1/c2)) + purpley1[0]*cos(-atan(-c1/c2))
             vk = vk + 1
             purplex1.append(1)
             purpley1.append(1)
             purplex2.append(1)
             purpley2.append(1)
          elif(purpley[to]<=base1(xinc[to]) and purpley[to]>=base2(xinc[to]) and vk!=0):
             
             purplex1[vk]= (round(xinc[to]*a + purpley[to]*b + bluepurple(xinc[to],purpley[to])*c,ndigits = 20))
             purpley1[vk]= (round(xinc[to]*d + purpley[to]*e + bluepurple(xinc[to],purpley[to])*f,ndigits = 20))
             purplex2[vk] = purplex1[vk]*cos(-atan(-c1/c2)) - purpley1[vk]*sin(-atan(-c1/c2))
             purpley2[vk] = purplex1[vk]*sin(-atan(-c1/c2)) + purpley1[vk]*cos(-atan(-c1/c2))'''
