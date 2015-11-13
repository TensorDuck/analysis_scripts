"""
Method for plotting the distirbution of Epsilon values from a simulation
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

def direc_run(args):
    cwd = os.getcwd()
    maxz = 0
    minz = 1000
    if not args.reference == None:
        ref_data = np.loadtxt(args.reference)
    else:
        ref_data = None
    for temp in args.temps:
        for i in np.arange(args.iterations[0], args.iterations[1]+1, 1):
            fwd = "%s/%d/%s/iteration_%d/%d_0"% (args.filedir, temp, args.subdir, i, temp) 
            pairs = np.loadtxt("%s/%s"%(fwd, args.pairs_name), skiprows=args.pairs_skip)[:,0:2]
            params = np.loadtxt("%s/%s"%(fwd, args.mparams_name))
            y = pairs[:,0]
            x = pairs[:,1]
            title = "T%d-I%d" % (temp, i)
            os.chdir(args.savedir)
            plot_it(x,y, params, title, args, reference=ref_data)
            plot_spread(params, title, args)
            os.chdir(cwd)
            if np.max(params) > maxz:
                maxz = np.max(params)
            if np.min(params) < minz:
                minz = np.min(params)
    
    os.chdir(args.savedir)
    fmax = open("max_value.dat", "a")
    fmax.write("For iterations %d-%d and temps %d-%d is: max = %f min = %f\n" % (args.iterations[0], args.iterations[1], np.min(args.temps), np.max(args.temps), maxz, minz))
    
    
def single_run(args):
    cwd = os.getcwd()
    
    pairs = np.loadtxt(args.pairs_name, skiprows=args.pairs_skip)[:,0:2]
    params = np.loadtxt(args.mparams_name)
    if not args.reference == None:
        ref_data = np.loadtxt(args.reference)
    else:
        ref_data = None
        
    y = pairs[:,0]
    x = pairs[:,1]
    os.chdir(args.savedir)
    plot_spread(params, args.title, args)
    plot_it(x,y, params, args.title, args, reference=ref_data)
    os.chdir(cwd)

def plot_spread(epsilons, title, args):
    plt.figure()
    plt.hist(epsilons, 50, alpha=0.75)
    plt.xlabel("epsilon",fontsize=20)
    plt.ylabel("number",fontsize=20)
    plt.title("spread of epsilons", fontsize=20)
    plt.savefig("%s_spread.png"%title)

def plot_it(x, y, z, title, args, reference=None):
    ctype = args.ctype
    
    if args.log:
        z = np.log(z)
        center = 0.0
    else:
        center = 1.0
    
    if args.zmin == None:
        zmin = np.min(z)
    else:
        zmin = args.zmin
    if args.zmax == None:
        zmax = np.max(z)
    else:
        zmax = args.zmax
        
    if args.center:
        zmin = 2*center-(zmax)
    
    z[z>zmax] = zmax #for setting all above a certain value mapped to the maximum value, for ease of plotting
    
    #make the array
    #set edges
    edges = np.array([0])
    edges = np.append(edges, np.arange(0.5, args.max_residue, 1.0))
    edges = np.append(edges, args.max_residue)
    residues = np.zeros((args.max_residue+1, args.max_residue+1))
    
    if reference == None:
        for i in range(np.shape(z)[0]):
            residues[x.astype(int)[i], y.astype(int)[i]] = z[i]
            residues[y.astype(int)[i], x.astype(int)[i]] = z[i]
    else:
        for i in range(np.shape(z)[0]):
            residues[x.astype(int)[i], y.astype(int)[i]] = z[i]
        for i in range(np.shape(reference)[0]):
            residues[reference.astype(int)[i,1],reference.astype(int)[i,0]] = 1    
    
    residues_masked = np.ma.masked_where(residues==0, residues)
    
    plt.figure()
    qmesh = plt.pcolormesh(edges, edges, residues_masked, vmin=zmin, vmax=zmax, cmap=ctype)
    #cp = plt.scatter(x, y, s=5, c=z, cmap=ctype, marker='o', linewidth=0., vmin=zmin, vmax=zmax)
    #cp = plt.scatter(y, x, s=5, c=z, cmap=ctype, marker='o', linewidth=0., vmin=zmin, vmax=zmax)

    plt.axis([0, args.max_residue, 0, args.max_residue])
    cb = plt.colorbar(qmesh)
    cb.set_label("epsilons", size=16)
    plt.xlabel("i", size=18)
    plt.ylabel("j", size=18)
    plt.title(title, size=20)
    
    plt.savefig("%s.png"%title)
    plt.show()



def sanitize_args(args):
    if not os.path.isdir(args.savedir):
        os.mkdir(args.savedir)
    
    
    return args

def get_args():
    par = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    par.add_argument("--filedir", default=os.getcwd(), type=str, help="Starting location")
    par.add_argument("--mparams_name", default="model_params", type=str, help="name of the model_params file(s)")
    par.add_argument("--pairs_name", default="pairwise_params", type=str, help="name of the file containing the contact pairs")
    par.add_argument("--pairs_skip", default=1, type=int, help="specify number of rows to skip when reading the pairs file")
    par.add_argument("--zmin", default=None, type=float, help="minimum value of z for color mapping")
    par.add_argument("--zmax", default=None, type=float, help="maximum value of z for color mapping")
    par.add_argument("--ctype", default="jet", type=str, help="type of color-map to use")
    par.add_argument("--center", default=False, action="store_true")
    par.add_argument("--log", default=False, action="store_true")
    par.add_argument("--max_residue", default=292, type=int)
    par.add_argument("--only", default=None, args="+", help="specify which types of atom contacts you want")
    par.add_argument("--reference", default=None, help="Specify a reference set of pairs")
    
    parser = argparse.ArgumentParser(description="For Deciding how to plot the results")
    sub = parser.add_subparsers(dest="method")
    
    
    
    #for analyzing directories
    direc_sub = sub.add_parser("direc", parents=[par], help="for plotting a whole directory")
    direc_sub.add_argument("--iterations", nargs=2, type=int, help="iteration range for plotting")
    direc_sub.add_argument("--temps", nargs="+", type=int, help="temperatures for plotting")
    direc_sub.add_argument("--savedir", default="%s/eps_plots"%os.getcwd(), type=str, help="Saving location")
    direc_sub.add_argument("--subdir", default="1PBQ", type=str, help="subdir of the simulations")
    
    single_sub = sub.add_parser("single", parents=[par], help="for plotting only a single plot")
    single_sub.add_argument("--savedir", default=os.getcwd(), type=str, help="Saving location")
    single_sub.add_argument("--title", default="test", type=str, help="Plot title and file-save name")
    
    args = parser.parse_args()
    
    args = sanitize_args(args)
    
    return args
    
if __name__ == "__main__":
    
    
    args = get_args()
    
    if args.method == "direc":
        direc_run(args)
    elif args.method == "single":
        single_run(args)
        
        
