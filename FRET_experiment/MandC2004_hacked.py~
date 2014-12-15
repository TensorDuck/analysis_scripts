"""
Method for running MandC2004, while loading contact epsilons from a separate file

"""

from project_tools.recipes import MatysiakClementi2004 as mc2004
from project_tools.recipes.MatysiakClementi2004 import MatysiakClementi2004 as mc_class

class MandC2004hack(mc_class):
    def __init__(self,paths):
        self.path = paths

if __name__ == "__main__":
    args, options = mc2004.get_args()
    options["Contact_Epsilons"] = np.loadtxt("epsilons.dat")
    mc_class(args, options)
    
    
