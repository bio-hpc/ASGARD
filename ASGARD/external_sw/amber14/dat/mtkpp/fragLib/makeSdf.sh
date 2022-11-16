#!/bin/sh

files="terminal.xml 2PtLinkers.xml 3PtLinkers.xml 4PtLinkers.xml 5PtLinkers.xml 3MemRings.xml 4MemRings.xml 5MemRings.xml 6MemRings.xml gt6MemRings.xml fusedRings.xml cores.xml hcaII.xml"

rm -rf fragment.sdf
touch fragLib.sdf
for f in ${files}
 do
   f_id=`echo $f | awk '{split($0,a,"."); print a[1]}'`
   stdLib2Sdf -i $f -o $f_id.sdf -a $f_id.log
   cat $f_id >> fragLib.sdf
 done
