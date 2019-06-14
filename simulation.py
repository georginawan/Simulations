from scipy.stats import binom
from scipy.stats import poisson
from scipy.stats import beta
import matplotlib.pyplot as plt
import numpy as np

def Input():
    custom = input("Define parameters by yourself?(Y/N):")
    if custom == "Y" or custom == "y":
        #true methylation rate p
        CoV_p = input("Constant rate p or variable rate p?(C/V):")
        if CoV_p == "C":
            p = float(input("True methylation rate p:"))
            variable_p = False
        elif CoV_p == "V":
            p = float(input("Average true methylation rate p:"))
            variable_p = True
        else:
            exit("Wrong input!")
        # coverage C 
        CoV_c = input("Constant coverage c or variable coverage c?(C/V):")
        if CoV_c == "C":
            c = int(input("coverage c:"))
            variable_c = False
        elif CoV_c == "V":
            c = int(input("average coverage c:"))
            variable_c = True
        else:
            exit("Wrong input!")
        #number of CpGs
        m = int(input("number of CpGs:"))
        #number of sampele per class(n1,n2)
        n1 = int(input("number of samples of class 1:"))
        n2 = int(input("number of samples of class 2:"))
    elif custom == "N" or custom == "n":
        p = 0.5
        c = 10
        m = 4
        n1 = n2 = 3
        variable_c = variable_p = False
    else:
        exit("Wrong input!")
    return p,c,m,n1,n2,variable_c,variable_p

#if coverage C is variable with given average, e.g.Poisson(n)
def variableCoverage(n,s,s1,s2,p):
    c = poisson.rvs(n,size = s)
    class1 = class2 = 0
    c1_sum = c2_sum = 0
    for i in range(s):
        r = binom.rvs(c[i],p)
        if i < s1:
            class1 += r
            c1_sum += c[i]
        else:
            class2 += r
            c2_sum += c[i]
    #average methylation level for class 1 and class 2
    aml1 = class1/c1_sum
    aml2 = class2/c2_sum
    return abs(aml1-aml2) # =delta

def constantCoverage(c,p,s,s1,s2):
    r = binom.rvs(c,p,size = s)
    class1 = class2 = 0
    for i in range(s):
        if i < s1:
            class1 += r[i]
        else:
            class2 += r[i]
    aml1 = class1/s1/c
    aml2 = class2/s2/c
    return abs(aml1-aml2) # =delta

#if rate p is variable, using beta distribution, e.g. a=b=2, average p =0.5
def variableP(c,p,s,s1,s2):
    a = 2
    b = a/p -a
    p_array = beta.rvs(a,b,size = s)
    class1 = class2 = 0
    for i in range(s):
        r = binom.rvs(c,p_array[i])
        if i < s1:
            class1 += r
        else:
            class2 += r
    aml1 = class1/s1/c
    aml2 = class2/s2/c
    return abs(aml1-aml2)  # = delta

p,c,m,n1,n2,vc,vp = Input()
s1 = n1*m
s2 = n2*m
s = s1+s2
times = 0
deltas = list()
while times < 1000000:
    if vc == True:
        delta = variableCoverage(c,s,s1,s2,p)
    elif vp == True:
        delta = variableP(c,p,s,s1,s2)
    else:
        delta = constantCoverage(c,p,s,s1,s2)
    deltas.append(delta)
    times += 1
    #if times%10000 == 0:
        #print(times)
deltas = sorted(deltas, reverse = True)
plt.plot(range(1000000),deltas)
plt.xlabel("index")
plt.ylabel("delta")
plt.show()
