""" Script for reading an outputted Jacobian Matrix and printing a color map

Justin Chen
November 16, 2014
"""
 
import numpy as np
import model_builder as mdb
import os as os_n
import matplotlib.pyplot as plt
import matplotlib.mpl as mpl

def read_jacobian(temp):
    cwd = os.getcwd()
    fwd = "/home/jchen/projects/2014/10_31_1PB7/1PB7"

    os.chdir(fwd)
    modelA = mdb.check_inputs.load_model(fwd, False)
    
    os.chdir(cwd)

    contacts =  modelA.contacts
    Jacobian = np.loadtxt(cwd+"/Jac_1PB7_"+temp+"_0.txt")
    print "Jacobian Shape for model T = "+temp+" is = " + str(np.shape(Jacobian))
    x = np.arange(np.shape(Jacobian)[1])

    #cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', ['blue','red','black'], 256)

    '''
    img = plt.imshow(Jacobian, interpolation="nearest", cmap=cmap, origin=0, aspect="auto")
    plt.colorbar(img, ticks = [0, 5, 10,15])
    plt.show()


    for i in np.arange(np.shape(Jacobian)[0]):
        plt.figure()
        plt.plot(x, Jacobian[i,:])
        plt.title("Figure " + str(i))
    '''
    print Jacobian

    umat, smat, vmat = np.linalg.svd(Jacobian)

    print np.shape(smat)
    print smat

    #print smat
    
    return smat[:]

    os.chdir(cwd)

if __name__ == "__main__":
    temparray = ["135", "140", "145", "150", "155", "160"]
    print "Starting the loop array sequence"
    smat_range = []
    for i in temparray:
        smat_range.append(read_jacobian(i))
        print "Finished Temperature: " + i +" Kelvin"
    
    colors = ["k","b","g","r","c","m","y"]
    count = 0
    for i in smat_range:
        plt.plot(np.arange(np.shape(i)[0]), i/i[0], color=colors[count], linewidth=2, label=temparray[count])
        count += 1
    plt.axis([0, 18, 10**-16, 100000])
    plt.yscale('log')
    plt.legend()
    plt.show()
    print "Don't panic, trust Data and disregard Riker"



