import numpy as np
import os

analysis_file = os.path.realpath(__file__)
stuff = analysis_file.split("/")
analysis_dir = ""
for thing in stuff[1:-1]:
    analysis_dir += "/%s" % thing

amino_acid_3code_noncaps = {"Cys":0, "Met":1, "Phe":2, "Ile":3, "Leu":4, "Val":5,
    "Trp":6,"Tyr":7, "Ala":8, "Gly":9, "Thr":10, "Ser":11, "Asn":12, "Gln":13,
    "Asp":14, "Glu":15, "His":16, "Arg":17, "Lys":18, "Pro":19}

amino_acid_3code = {key.upper():amino_acid_3code_noncaps[key] for key in amino_acid_3code_noncaps}

amino_acid_1code_to_3code = {"G":"Gly", "C":"Cys", "A":"Ala", "M":"Met", "V":"Val",
    "K":"Lys", "L":"Leu", "R":"Arg", "I": "Ile", "H":"His", "F":"Phe",
    "W":"Trp", "P":"Pro", "D":"Asp", "S":"Ser", "E":"Glu", "T":"Thr",
    "N":"Asn", "Y":"Tyr", "Q":"Gln"}

amino_acid_1code_to_3code_uppercase = {key:amino_acid_1code_to_3code[key].upper() for key in amino_acid_1code_to_3code}
amino_acid_3code_to_1code_uppercase = {amino_acid_1code_to_3code_uppercase[key]:key for key in amino_acid_1code_to_3code_uppercase}

mj_contact_energies = np.loadtxt("%s/Miyazawa-Jernigan_Potentials.dat" % analysis_dir, comments="#")

assert np.shape(mj_contact_energies)[0] == 20
assert np.shape(mj_contact_energies)[1] == 20
