BEGIN{impr=0.0000}


#CHARMM Vars
/ENER EXTERN>/{
  cVdwExt=$3;
  cEleExt=$4;
#  print "charmm", $3,$4;
}
/ENER IMAGES>/{
  cVdwImg=$3;
  cEleImg=$4;
#  print "charmm", $3,$4;
}
/ENER EWALD>/{
  cEwaldRec=$3;
  cEwaldSelf=$4;
  cEwaldAdj=$5;
#  print "charmm", $3,$4,$5;
}



#AMBER vars
/\(self\)/{aEwaldSelf=$5; }
/\(rec\)/{aEwaldRec=$3; }
/\(dir\)/{aEwaldDir=$3; }
/\(adj\)/{aEwaldAdj=$3; }

/ EELEC /{elec=$3; }
/ 1-4 NB /{nb14=$3;ee14=$8;vdw=$11; }



#BOND          ANGLE          DIHED          IMPRP               ELECT            VDW
#/ENERGY: / {NAMD_bnd=$3; NAMD_ang=$4; NAMD_dih=$5; NAMD_impr=$6; NAMD_elec=$7; NAMD_vdw=$8; }

END{
    printf( "            amb         chm    difference \n");
    printf( "                               (amb - chm)\n");
    printf ("ELEC  %11.4f %11.4f %11.4f \tamber:\n \t\t\t\t\t elec\t%11.4f\n \t\t\t\t\t ee14\t%11.4f\n",
	    (elec+eel14), 
            (cEwaldRec + cEwaldSelf + cEwaldAdj + cEleExt + cEleImg), 
            (elec+ee14)-(cEwaldRec + cEwaldSelf + cEwaldAdj + cEleExt + cEleImg),
            elec,ee14 );
    printf ("VDW   %11.4f %11.4f %11.4f \tamber:\n \t\t\t\t\t vdw\t%11.4f\n \t\t\t\t\t vdw14\t%11.4f\n",
	    (vdw+nb14),(cVdwExt + cVdwImg),(vdw+nb14)-(cVdwExt + cVdwImg), vdw,nb14);
    printf ("EWALD\n")
    printf (" self  %11.4f %11.4f %11.4f\n", (aEwaldSelf),(cEwaldSelf), (aEwaldSelf-cEwaldSelf) );
    printf (" rec   %11.4f %11.4f %11.4f\n", (aEwaldRec), (cEwaldRec), (aEwaldRec-cEwaldRec) )
    printf (" dir   %11.4f %11.4f %11.4f\n", (aEwaldDir), (cEleExt + cEleImg), (aEwaldDir-(cEleExt + cEleImg)) )
    printf (" adj   %11.4f %11.4f %11.4f\n", (aEwaldAdj), (cEwaldAdj), (aEwaldAdj-cEwaldAdj));
  print "   " 
}
