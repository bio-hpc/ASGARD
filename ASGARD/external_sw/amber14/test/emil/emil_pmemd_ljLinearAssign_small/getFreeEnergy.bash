#!/bin/bash

nlambdas=8
lambdas="0.0 0.1 0.3 0.5 0.7 0.900 0.999 1.0"

##get the reference free energy from one of the outfiles
aRef=$(grep A-Solid lj2_1.0/emil.log | awk '{print $3}' | tail -1)

chunk=1000
size=$chunk


rm -f gF*.dat
for l in $lambdas
do
    rm -f gF$l.dat

    grep nstep lj2_$l/emil.log     |\
          awk -v l=$l '{print l, $6,$8}' > gF$l.dat
done




files=$(ls gF*.dat)
#echo $files

echo "#t  Aref  dA_mol dA_ref dA_tot" > A.dat
paste $files | awk -v Aref=$aRef -v n=$nlambdas \
 'NF==3*n{lf1=$2;lf2=$3;ll=0.0;t1=0.0;t2=0.0;\
          for(i=3;i<=NF;i+=3){l=$(i-2);f1=$(i-1);f2=$i;\
                              t1+=(0.5*(l-ll)*(f1+lf1));lf1=f1;\
                              t2+=(0.5*(l-ll)*(f2+lf2));lf2=f2;ll=l} print NR,Aref,t1,t2,t1+t2}' >> A.dat


##get the mean free energy
meanA=$(grep -v "#" A.dat | awk '{t+=($2-$5);c++}END{print t/c}')

echo "Estimated free energy (should be -1 + small d): $meanA"


