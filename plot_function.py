import numpy as np
import matplotlib.pyplot as plt

def functionLJ(x,eps,sigma):
    return eps*(5*((sigma/x)**12) - 6*((sigma/x)**10))

def functionRep(x,eps,sigma):
    return eps*((sigma/x)**12)

if __name__ == "__main__":
    x = np.arange(0.01, 2, 0.01)
    yreg = np.zeros(len(x))
    ysmall = np.zeros(len(x))
    ysmalls = np.zeros(len(x))
    yrep = np.zeros(len(x))

    for i in range(len(x)):
        yreg[i] = functionLJ(x[i], 1.0, 0.4)
        ysmall[i] = functionLJ(x[i], 0.1, 0.4)
        ysmalls[i] = functionLJ(x[i],0.01,0.4)
        yrep[i] = functionRep(x[i], 1.0, 0.3)
 
    plt.plot(x,yreg,color="b")
    plt.plot(x,ysmall,color="r")
    plt.plot(x,ysmalls,color="g")
    plt.plot(x,yrep,color="k")
    plt.axis([0, 2, -10, 100])
    plt.savefig("test.png")
    plt.show()
