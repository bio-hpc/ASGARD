#!/bin/bash
##values taken from sander constants.F90
amber_BOLTZ=1.380658d-23
amber_AVOGADRO=6.0221367d+23
amber_JPKC=4.184d+3

##derived value for kB in amber units:
amber_KB=$(echo 1.380658 6.0221367 4184.0 |\
                awk '{print ($1 * $2) / $3}')

##build a cubic box full of Lennard-Jonesium using leap!

##coexistence parameters from Ferreira et al:
Tstar=0.7793
Pstar=1.0
n_particles=2

rhoStar_fluid=0.8699
rhoStar_fcc=0.9732
rCut_star=6.0
epsilon_star=1.0
sigma_star=1.0

Vstar=$(echo $n_particles $rhoStar_fcc | awk '{print $1/$2}')
Lstar=$(echo $Vstar | awk '{print $1**(1.0/3.0)}')

Astar=$(echo $epsilon_star $sigma_star | awk '{print 4.0*$1*($2**12)}')
Bstar=$(echo $epsilon_star $sigma_star | awk '{print 4.0*$1*($2**6)}')

rMinStar=$(echo $sigma_star | awk '{print $1*(2.0**(1./6.))}')

# Amber "Lennard Jones Radius" is half the distance between 
# the two particles at minimum energy epsilon.
halfRminStar=$(echo $rMinStar | awk '{print $1*0.5}')

echo "source leaprc.ff99SB"       >  makeLJ.tls
echo "m=sequence { WAT }"         >> makeLJ.tls
echo "solvatebox m TIP3PBOX 5"    >> makeLJ.tls
echo "savepdb m watBox.pdb"       >> makeLJ.tls
echo "quit"                       >> makeLJ.tls
 

export AMBERHOME=/usr/local/amber11
tleap -f  makeLJ.tls
export AMBERHOME=/home/josh/amberCheckout/amber


n_whole_cells=$[ $n_particles / 4 ]
n_extra_parts=$[ $n_particles % 4 ]

if [ "$n_extra_parts" != "0" ]
then
    echo "WARNING: non-integer number of fcc cells"
    exit
fi

cells_nx=$(echo $n_whole_cells | awk '{print $1**(1.0/3.0)}')
echo "CELLS_NX == $cells_nx"

cells_x=$(echo $n_whole_cells $Lstar | awk '{print $2/$1**(1.0/3.0)}')
echo "CELLS_X == $cells_x"


grep "O   WAT" watBox.pdb | head -n $n_particles |\
  awk -v L=$Lstar -v dx=$cells_x  \
                   'BEGIN{x1=-0.5*L;y1=-0.5*L;z1=-0.5*L;c=1}\
                    {x2=x1+0.5*dx;y2=y1+0.5*dx;z2=z1;\
                     x3=x1+0.5*dx;y3=y1;       z3=z1+0.5*dx;\
                     x4=x1;       y4=y1+0.5*dx;z4=z1+0.5*dx;\
                     tail=substr($0,56);\
                     print "ATOM  "    sprintf("%5d",c)\
                         "  O   WAT " sprintf("%5d",c++) "     "\
                         sprintf("%+7.3f %+7.3f %+7.3f ",x1,y1,z1) tail;\
                     print "ATOM  "    sprintf("%5d",c)\
                         "  O   WAT " sprintf("%5d",c++) "     "\
                         sprintf("%+7.3f %+7.3f %+7.3f ",x2,y2,z2) tail;\
                     print "ATOM  "    sprintf("%5d",c)\
                         "  O   WAT " sprintf("%5d",c++) "     "\
                         sprintf("%+7.3f %+7.3f %+7.3f ",x3,y3,z3) tail;\
                     print "ATOM  "    sprintf("%5d",c)\
                         "  O   WAT " sprintf("%5d",c++) "     "\
                         sprintf("%+7.3f %+7.3f %+7.3f ",x4,y4,z4) tail;\
                     x1+=dx;\
                     if(x1>=0.5*L){x1=-0.5*L;y1+=dx;\
                                       if(y1>=0.5*L){y1=-0.5*L;z1+=dx}}}'\
  | head -n $n_particles > hexBox.pdb


echo "source leaprc.ff99SB"       >  makeLJ.tls
echo "m=loadpdb hexBox.pdb"       >> makeLJ.tls
echo "solvatebox m TIP3PBOX 0"    >> makeLJ.tls
echo "saveamberparm m hexBox.top hexBox.crd" >> makeLJ.tls
echo "quit"                       >> makeLJ.tls
 

export AMBERHOME=/usr/local/amber11
tleap -f  makeLJ.tls
export AMBERHOME=/home/josh/amberCheckout/amber

###now strip the hydrogens 
###and set all charges to zero: hey presto! Lennard-Jones in natural units
rm -f ljBox.top
python $AMBERHOME/AmberTools/bin/parmed.py -p hexBox.top <<EOF
strip @H*
change charge * 0.0
change mass   * 1.0
changeLJSingleType * $halfRminStar $epsilon_star
printDetails @1
printLJMatrix *
parmout ljBox.top
EOF


##grep -A 2 "A coefficient   B coefficient" > lj_coeffs.dat

###now strip the hydrogens from the crd file
ptraj hexBox.top <<EOF
trajin hexBox.crd restart box
strip @H*
trajout hexBox.crd restart 
EOF
mv  hexBox.crd.1 ljBox.crd

##now characterise the system in reduced units:
mass=1.0 ##this was set in the parm.

# Amber temperature T_amber is Tstar*epsilon/k_B
# i.e. T* = k_B T_amber/epsilon
T_amber=$(echo $Tstar $epsilon_star $amber_KB | awk '{print $1*$2*$3}')

# P_amber works out equal to P star.
P_amber=$(echo $Pstar $epsilon_star $sigma_star | awk '{print $1*$2/($3**3)}')


echo "Set amber temperature to $T_amber to get Tstar=$Tstar"
echo "Set amber pressure to $P_amber to get Pstar=$Pstar"
echo "Should converge to give volume: $Vstar"


##fix up the box size
c=$(wc -l  ljBox.crd | awk '{print $1-1}')
head -n $c ljBox.crd > c
echo $Lstar | awk \
 '{printf " %11.7f %11.7f %11.7f  90.0000000  90.0000000  90.0000000\n",$1,$1,$1}'\
           >> c
mv c ljBox.crd
