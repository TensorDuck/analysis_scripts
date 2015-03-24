""" script for making a .ndx file """

import argparse
import numpy as np

def gen_index(name, labels, ranges):
    f = open("%s.ndx"%name,"w")
    for idx, label in enumerate(labels):
        f.write("[ %s ]\n" % label)
        for ran in ranges[idx]:
            
            values = np.arange(ran[0], ran[1]+1, 1)
            for i in values:
                f.write("%d\n" % i)
    f.close()    

def get_args():
    parser = argparse.ArgumentParser(description="specify ranges and labels")
    
    parser.add_argument("--labels", nargs="+", help="Specify names of labels, default is to do A, B, C, etc.")
    parser.add_argument("--ranges", nargs="+", help="Specify the residue ranges")
    
    args = parser.parse_args()
    
    return args

if __name__ == "__main__":
    args = get_args()
    
    name="lobes"
    labels=["A", "B"]
    
    ranges=[  [ [1, 143],[249, 285] ], [ [144,248],[286,292] ]  ]
    
    gen_index(name, labels, ranges)
