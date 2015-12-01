""" for adding nonnative contacts to the mix"""

fold = open("pairwise_params", "r")
fnew = open("pairwise_params_nonnative", "w")
fparams = open("model_params_nonnative", "w")
pdb = open("Native.pdb", "r")

residue_list = []
for line in pdb:
    residue_list.append(pdb.strip.split()[3])

fnew.write(fold.readline())
fparams.write("# model params\n")
maxresidue = 58

residue_volumes = {'ALA': 67.0, 'ARG': 148.0, 'ASN': 96.0,
                'ASP': 91.0, 'CYS': 86.0, 'GLN': 114.0,
                'GLU': 109.0, 'GLY': 48.0, 'HIS': 118.0,
                'ILE': 124.0, 'LEU': 124.0, 'LYS': 135.0,
                'MET': 124.0, 'PHE': 135.0, 'PRO': 90.0,
                'SER': 73.0, 'THR': 93.0, 'TRP': 163.0,
                'TYR': 141.0, 'VAL': 105.0}

residue_volume_mod = {key:(residue_volumes[key]/ ((4.0*3.14159)/3.0))**(1.0/3.0) for key in residue_volumes} 

residue_radii_new = {key:residue_volume_mod*1.4 for key in residue_volume_mod}

count = 0
parcount = 0
cur_list = fold.readline().strip().split()

for i in range(maxresidue):
    for j in range(maxresidue):
        radii_a = residue_radii_new[residue_list[i+1]]
        radii_b = residue_radii_new[residue_list[j+1]]
        excvol = (radii_a*radii_b)**0.5
        if i+1 == native_contact[0] and  j+1 == native_contact[1]:
            native = True
            minima = float(cur_list[5])
            fold.readline()
            cur_list = fold.readline().strip().split()
            try:
                native_contact = [int(cur_list[0]), int(cur_list[1])]
            except:
                native_contact=[max_residue+1, max_residue+1]
            
        else:
            native=False
            minima = excvol+0.2
        if j-i > 3:
            fnew.write("%5d%5d%7d%5d%11.5f%12.5f%12.5f\n" %(i+1, j+1, parcount, 8, excvol, minima, 0.05)
            parcount += 1 
            fnew.write("%5d%5d%7d%5d%11.5f%12.5f\n" %(i+1, j+1, parcount, 4, minima, 0.05)
            parcount += 1
            if native:
                fparams.write("   1.00000\n   1.00000\n")
            else:
                fparams.write("   1.00000\n   0.00000\n")



fold.close()
fnew.close()
fparams.close()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        