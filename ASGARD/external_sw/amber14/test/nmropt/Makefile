# This is a Makefile for the NMRopt tests.  This lower-level Makefile is used to
# simplify the division between those tests that can be run with pmemd and those
# that can only be run with sander.  This just tests nmropt=1 so far with most
# of the available options.

all: pmemd_compat pmemd_incompat

pmemd_compat: pmemd_compat_gb pmemd_compat_pme

pmemd_incompat: pmemd_incompat_gb pmemd_incompat_pme

pmemd_compat_gb:
	-cd gb/angle && ./Run.nmropt_1angle_gb
	-cd gb/distance && ./Run.dist_gb
	-cd gb/distance_COM && ./Run.distCOM_gb
	-cd gb/jar_distance && ./Run.jar_gb
	-cd gb/jar_distance_COM && ./Run.jar_gb
	-cd gb/jar_torsion && ./Run.jar_torsion
	-cd gb/nmropt_1_torsion && ./Run.nmropt_1_torsion
	-cd gb/tautp && ./Run.nmropt_1tautp_gb
	-cd gb/temp && ./Run.nmropt_1temp_gb

pmemd_incompat_gb:
	-cd gb/2xangle && ./Run.nmropt_1_2xangle_gb
	-cd gb/jar_comtorsion && ./Run.jar_comtorsion
	-cd gb/jar_lcod && ./Run.lcod_jar
	-cd gb/jar_lcod_dist && ./Run.jar_lcod
	-cd gb/lcod && ./Run.lcod
	-cd gb/lcod_com && ./Run.lcod_com
	-cd gb/nmropt_1_2xtorsion && ./Run.nmropt_1_2xtorsion

pmemd_compat_pme:
	-cd pme/angle && ./Run.nmropt_1angle_pbc
	-cd pme/distance && ./Run.dist_pbc
	-cd pme/distance_COM && ./Run.distCOM_pbc
	-cd pme/jar_torsion && ./Run.jar_torsion
	-cd pme/jar_distance && ./Run.jar_pbc
	-cd pme/jar_distance_COM && ./Run.jar_pbc
	-cd pme/nmropt_1_torsion && ./Run.nmropt_1_torsion
	-cd pme/tautp && ./Run.nmropt_1tautp_pbc
	-cd pme/temp && ./Run.nmropt_1temp_pbc

pmemd_incompat_pme:
	-cd pme/2xangle && ./Run.nmropt_1_2xangle_pbc
	-cd pme/jar_comtorsion && ./Run.jar_comtorsion
	-cd pme/jar_lcod && ./Run.jar_lcod
	-cd pme/jar_lcod_dist && ./Run.lcod_dist
	-cd pme/lcod && ./Run.lcod_pme
	-cd pme/lcod_com && ./Run.lcod_com
	-cd pme/nmropt_1_2xtorsion && ./Run.nmropt_1_2xtorsion
