#!/usr/bin/env python
#https://gist.github.com/devries/11405101
import math
import sys

def main(files=None,array1=None,array2=None):
    if argv is None:
        argv = sys.argv

    if len(argv) != 3:
        print("Usage:",sys.argv[0],"<datafile1> <datafile2>", file=sys.stderr)
        sys.exit(1)

    print("The KS test measures the maximum deviation (D) in the cumulative")
    print("probabilities. Then the probability that the data sets are drawn")
    print("from the same distribution can be calculated. Small values")
    print("of the probability indicate that the distrubution of sets one")
    print("and two are significantly different.")
    f = open(sys.argv[1],'r')
    d1 = []
    for ln in f:
        try:
            d1.append(float(ln))
        except:
            pass
    f.close()

    f = open(sys.argv[2],'r')
    d2 = []
    for ln in f:
        try:
            d2.append(float(ln))
        except:
            pass
    f.close()

    (d,prob,ne) = kstest(d1,d2)
    print("D =",d)
    print("Prob =",prob)
    print("Ne =",ne)
    return 0

def kstest(datalist1, datalist2):
    n1 = len(datalist1)
    n2 = len(datalist2)
    datalist1.sort()
    datalist2.sort()

    j1 = 0
    j2 = 0
    d = 0.0
    fn1=0.0
    fn2=0.0
    while j1<n1 and j2<n2:
        d1 = datalist1[j1]
        d2 = datalist2[j2]
        if d1 <= d2:
            fn1 = (float(j1)+1.0)/float(n1)
            j1+=1
        if d2 <= d1:
            fn2 = (float(j2)+1.0)/float(n2)
            j2+=1
        dtemp = math.fabs(fn2-fn1)
        if dtemp>d:
            d=dtemp

    ne = float(n1*n2)/float(n1+n2)
    nesq = math.sqrt(ne)
    prob = ksprob((nesq+0.12+0.11/nesq)*d)
    return d,prob,ne

def ksprob(alam):
    fac = 2.0
    sum = 0.0
    termbf = 0.0

    a2 = -2.0*alam*alam
    for j in range(1,101):
        term = fac*math.exp(a2*j*j)
        sum += term
        if math.fabs(term) <= 0.001*termbf or math.fabs(term) <= 1.0e-8*sum:
            return sum
        fac = -fac
        termbf = math.fabs(term)

    return 1.0

            
