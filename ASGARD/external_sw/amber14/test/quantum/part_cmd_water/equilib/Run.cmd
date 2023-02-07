#!/bin/csh -f
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

set sander = "../../../../bin/sander.LES"

if( $?TESTsander ) then
   set sander = $TESTsander
endif

if( ! $?DO_PARALLEL ) then
        setenv DO_PARALLEL " "
endif


$DO_PARALLEL $sander -O -i cmd.in -p h2o_les.top -c h2o_les.crd -o cmd.out < /dev/null || goto error
../../../dacdif cmd.out.save cmd.out
/bin/rm -f mdcrd mdinfo restrt fort.277

exit(0)

error:
echo "  ${0}:  Program error"
exit(1)


