""" Script for calculating the next set of epsilon values for a FRET
 distribution with Arbitrary number of pairs

assumes you are in the subdirectory containing the model.info file
"""
 
import numpy as np
import model_builder as mdb
import os as os
import project_tools.parameter_fitting.FRET.compute_Jacobian as compJ

GAS_CONSTANT_KJ_MOL = 0.0083144621

def Calc_Jacobian(temp):

    ##Define the known pairs
    pairs = np.array([[115,193]])
    pairs -= np.ones(pairs.shape)
    cwd = os.getcwd()
    spacing = 0.1

    modelA = mdb.check_inputs.load_model(cwd, False)

    #print modelA.residues[192]

    beta = float(temp) * GAS_CONSTANT_KJ_MOL * -1.0

    os.chdir("Tf_0/"+temp+"_0" )

    Jacobian, simparams = compJ.compute_Jacobian_for_directory(modelA, beta, pairs, spacing)

    np.savetxt("/home/jchen/analysis/2014/11-21-Jacobian_Attempts/1_oldJac_1PB7_"+temp+"_0.txt", Jacobian)
    np.savetxt("/home/jchen/analysis/2014/11-21-Jacobian_Attempts/1_oldsimparams_1PB7_"+temp+"_0.txt", simparams)
    
    os.chdir(cwd)
    
if __name__ == "__main__":
    temparray = ["135", "140", "145", "150", "155", "160"]
    print "Starting the loop array sequence"
    for i in temparray:
        Calc_Jacobian(i)
        print "Finished Temperature: " + i +" Kelvin"
    
    print "Open the door Hal"

