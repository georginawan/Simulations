from scipy.stats import binom
from scipy.stats import poisson
from scipy.stats import beta
import matplotlib.pyplot as plt
import numpy as np

#input mean and variance, output alpha and beta.
def mv_to_ab(mean,var):
    alpha = mean*mean*(1-mean)/var - mean
    beta = (1-mean)/mean*alpha
    return alpha,beta

def variableCP(a1,b1,a2,b2,s1,s2,c1,c2):
    p1_array = beta.rvs(a1,b1,size = s1)
    p2_array = beta.rvs(a2,b2,size = s2)
    c1_array = poisson.rvs(c1,size = s1)
    c2_array = poisson.rvs(c2,size = s2)
    class1 = class2 = 0
    c1_sum = c2_sum = 0
    for i in range(s1+s2):
        if i < s1:
            r = binom.rvs(c1_array[i],p1_array[i])
            class1 += r
            c1_sum += c1_array[i]
        else:
            r = binom.rvs(c2_array[i-s1],p2_array[i-s1])
            class2 += r
            c2_sum += c2_array[i-s1]
    aml1 = class1/c1_sum
    aml2 = class2/c2_sum
    return abs(aml1-aml2)

def deltas():
    deltas1 = list()
    deltas2 = list()
    times = 0
    a1,b1 = mv_to_ab(0.3,0.01)
    a2,b2 = mv_to_ab(0.7,0.01)
    while times < 1000000:
        delta1 = variableCP(a1,b1,a2,b2,12,12,10,10) #real DMR
        delta2 = variableCP(2,2,2,2,12,12,10,10)     #random
        deltas1.append(delta1)
        deltas2.append(delta2)
        times += 1
        if times%10000 == 0:
            print(times)
    deltas1 = sorted(deltas1, reverse = True)
    deltas2 = sorted(deltas2, reverse = True)
    return deltas1,deltas2

def picture(deltas1,deltas2):
    plt.plot(np.array(range(1000000))/1000000,deltas1,color="deepskyblue")
    plt.plot(np.array(range(1000000))/1000000,deltas2,color="pink")
    plt.legend(["real DMR","random DMR"])
    plt.xlabel("index")
    plt.ylabel("delta")
    plt.show()

def t_to_i(deltalist,t):
    count = 0
    for i in range(len(deltalist)):
        if deltalist[i] >= float(t):
            count += 1
        else:
            break
    return "index = "+str(count/len(deltalist))

def i_to_t(deltalist,i):
    t = deltalist[int(i*len(deltalist))]
    return "t = "+str(t)

def test(deltas1,deltas2):
    input_ = input("What do you want to test?(t=0.4/i=0.8):") # t for delta threshold, i for index, enter space to get out of loop.
    if input_ == " ":
        exit("See you!")
    else:
        input_ = input_.replace(" ","")   
        if input_[0] == "t":
            t = input_[3:]      
            print("Test result for",input_,"is:",)
            print("real DMR:",t_to_i(deltas1,t))
            print("random DMR:",t_to_i(deltas2,t))        
        elif input_[0] == "i":
            i = float(input_[2:])
            print("Test result for",input_,"is:",)
            print("real DMR:",i_to_t(deltas1,i))
            print("random DMR:",i_to_t(deltas2,i))
        else:
            exit("Wrong input!")

def main():
    deltas1,deltas2 = deltas()
    picture(deltas1,deltas2)
    while True:
        test(deltas1,deltas2)

main()
