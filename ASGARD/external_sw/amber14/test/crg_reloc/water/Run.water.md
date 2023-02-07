#!/bin/csh -f
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

if( ! $?TESTsander ) set TESTsander = "../../../bin/sander"

if( ! $?DO_PARALLEL ) then
   setenv DO_PARALLEL " "
endif

$DO_PARALLEL $TESTsander -O -i md.in -o md.out -c md.inpcrd || goto error

../../dacdif md.out.save md.out
/bin/rm -f mdinfo restrt mdcrd
exit(0)

error:
echo "  ${0}:  Program error"
exit(1)
