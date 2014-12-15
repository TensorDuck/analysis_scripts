import model_builder.models.beadbead_to_params as bbp

ints, others = bbp.convert_beadbead("BeadBead.dat")
open("pairwise_params","w").write(ints)
open("model_params","w").write(others)
