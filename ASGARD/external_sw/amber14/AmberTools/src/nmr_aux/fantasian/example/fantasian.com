rm obs.out res.out
../bin/fantasian << EOF          
pcshifts.in                      # pseudocontact shifts file (input)
parm.pdb                         # pdb file (AMBER) (input)
obs.out                          # observed out file (output)
res.out                          # main out file (output)
1782            #  number of atoms for each structure
1               #  number of superimposed structures
FE              #  name of the paramagnetic center(s)
y               #  Would you like reference system(s)? (y/n)
1               #  number of paramagnetic centers
23              #  residue number which contains the origin of the ref. sys. 
358             #  atom number which fixes the x direction of the ref. sys.
357             #  atom number which fixes the origin of the ref. sys.
379             #  atom number of the other point in the plane x-y
3               #  number of independent simplex calculations 
EOF             
