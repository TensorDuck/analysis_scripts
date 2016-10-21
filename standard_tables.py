import numpy as np

amino_acid_3code = {"Cys":0, "Met":1, "Phe":2, "Ile":3, "Leu":4, "Val":5,
    "Trp":6,"Tyr":7, "Ala":8, "Gly":9, "Thr":10, "Ser":11, "Asn":12, "Gln":13,
    "Asp":14, "Glu":15, "His":16, "Arg":17, "Lys":18, "Pro":19}

aminio_acid_1code = {"G":"Gly", "C":"Cys", "A":"Ala", "M":"Met", "V":"Val",
    "K":"Lys", "L":"Leu", "R":"Arg", "I": "Ile", "H":"His", "F":"Phe",
    "W":"Trp", "P":"Pro", "D", "Asp", "S":"Ser", "E":"Glu", "T":"Thr",
    "N":"Asn", "Y":"Tyr", "Q":"Gln"}

mj_contact_energies = np.loadtxt("Miyazawa-Jernigan_Potentials.dat")

assert np.shape(mj_contact_energies)[0] == 20
assert np.shape(mj_contact_energies)[1] == 20
