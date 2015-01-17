"""
Simple method for merging files together for the results from a dmdmd simulation
"""

def merge(start, stop):
    #fout = open("iter%d-%d.gro"%(start,stop), "w")
    fwout = open("iter%d-%d.w"%(start,stop), "w")
    for i in range(start, stop+1, 2):
        #fin = open("iter%d.gro"%i,"r")
        fwin = open("iter%d.w"%i,"r")
        #for lines in fin:
           # fout.write(lines)
        for lines in fwin:
            fwout.write(lines)
        #fin.close()
        fwin.close()
    #fout.close()
    fwout.close()
    
if __name__ == "__main__":
    
    merge(6,8)
        
