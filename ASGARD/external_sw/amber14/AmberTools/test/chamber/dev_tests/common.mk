#Defines specific to the AMBER install

##To implement in future
#SANDER=$$AMBERHOME/bin/sander
#CHAMBER=$$AMBERHOME/bin/chamber
CHAMBER=../../../../bin/chamber

DACDIF=../../../dacdif
COMPENE=../comp_ene.awk
GDB=gdb

#needed for frcdump
#CHARMM=/server-home/mjw/code/CHARMM/CVS/c36a2/exec/em64t/charmm

#DEFAULT
#CHARMM=/server-home/netbin/c35b1_xxlarge_em64t_serial/exec/em64t/charmm
CHARMM="charmm"

#Default AMBER md engine
MD_ENGINE=../../../../bin/sander

##MPI related
MD_ENGINE_MPI=../../../../bin/sander.MPI

#MPIRUN=/server-home/netbin/mpi/mpich2-1.0.7-gfortran-4.1.2/bin/mpirun
MPIRUN=mpirun
