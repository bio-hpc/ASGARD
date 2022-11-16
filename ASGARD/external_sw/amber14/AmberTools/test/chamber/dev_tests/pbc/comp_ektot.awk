BEGIN{impr=0.0000}


#CHARMM Vars
/DYNA>     6000/{
  cEktot=$5;
  #print "charmm", $3,$4;
}

#AMBER vars
/ Etot   /{aEktot=$6; }



#BOND          ANGLE          DIHED          IMPRP               ELECT            VDW
#/ENERGY: / {NAMD_bnd=$3; NAMD_ang=$4; NAMD_dih=$5; NAMD_impr=$6; NAMD_elec=$7; NAMD_vdw=$8; }

END{
    printf( "            amb         chm    difference \n");
    printf( "                               (amb - chm)\n");
    printf (" Ektot   %11.4f %11.4f %11.4f\n", (aEktot), (cEktot), (aEktot-cEktot));
  print "   " 
}
