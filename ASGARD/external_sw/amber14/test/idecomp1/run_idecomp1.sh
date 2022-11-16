#!/bin/csh -f
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented
 
set sander = "../../bin/sander"
if( $?TESTsander ) then
   set sander = $TESTsander
endif
 
if( ! $?DO_PARALLEL ) then
        setenv DO_PARALLEL " "
endif
 
cat > decomp.in <<EOF
#decompose energy on a per-residue basis on the trajectory
 &cntrl
      imin = 5, igb = 5, gbsa = 2,
      ntx = 1, maxcyc = 1,
      ntc = 1, ntf = 3,
      ntb = 0, ntp = 0,
      ntwe = 0, ntpr = 500, ntwx=1,
      cut = 1000.0, idecomp = 1,
 /
Residues considered as REC
RRES  1  13
END
Residues considered as LIG
LRES  1  13
END
Residues to print
RES 1 13
END
END
EOF
 
touch dummy
$DO_PARALLEL $sander \
-O \
-i ./decomp.in \
-p ./test.top \
-c ./test.crd \
-o ./decomp.out \
-y ./test.traj \
-ref ./test.crd \
-inf ./decomp.info \
-r ./decomp.rst \
< dummy || goto error
 
../dacdif decomp.out.save decomp.out
/bin/rm -f decomp.in decomp.info  decomp.rst mdcrd  dummy
exit(0)
 
error:
echo "  ${0}:  Program error"
exit(1)
