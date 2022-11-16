#!/bin/csh -f
#TEST-PROGRAM nmode
#TEST-DESCRIP Newton Raphson minimization
#TEST-PURPOSE regression, basic
#TEST-STATE   ok

echo " "
echo NMODE: check Newton Raphson minimization
../../bin/nmode -O \
		-i nr.in \
		-c nmode_1.rst \
		-p parm.top \
		-o nr.out		|| goto error

../dacdif nr.out.save nr.out
/bin/rm -f restrt

exit(0)

error:
echo "  ${0}:  Program error"
echo " "
exit(1)

