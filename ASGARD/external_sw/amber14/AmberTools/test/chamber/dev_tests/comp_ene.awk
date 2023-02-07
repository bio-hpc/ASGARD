BEGIN{impr=0.0000}
/ENER INTERN>/{
  cbnd=$3;
  cang=$4;
  cub=$5;
  cdih=$6;
  cimpr=$7;
#  print "charmm", $3,$4,$5,$6,$7;
}
/ENER EXTERN/{cvdw=$3;celec=$4;}
/ENER CROSS/{cmap=$3;}

/ BOND /{bnd=$3; ang=$6; dih=$9}
/ VDWAALS /{vdw=$3;elec=$6;}
/ 1-4 / {nb14=$4;ee14=$8;}

#BOND          ANGLE          DIHED          IMPRP               ELECT            VDW
#/ENERGY: / {NAMD_bnd=$3; NAMD_ang=$4; NAMD_dih=$5; NAMD_impr=$6; NAMD_elec=$7; NAMD_vdw=$8; }

END{
    printf( "            amb         chm    difference \n");
    printf( "                               (amb - chm)\n");
    printf( "BOND  %10.4f %10.4f %10.4f \n",bnd,cbnd,bnd-cbnd);
    printf ("ANGLE %10.4f %10.4f %10.4f \tcharmm:\n \t\t\t\t\t ang\t%10.4f\n \t\t\t\t\t ub\t%10.4f \n",ang, cang+cub, ang-(cang+cub), cang,cub);
    printf ("DIHED %10.4f %10.4f %10.4f \tcharmm:\n \t\t\t\t\t dihe\t%10.4f\n \t\t\t\t\t impr\t%10.4f\n \t\t\t\t\t cmap\t%10.4f\n",
            dih,(cdih+cmap+cimpr),dih-(cdih+cmap+cimpr), cdih,cimpr,cmap);
    printf ("ELEC  %10.4f %10.4f %10.4f \tamber:\n \t\t\t\t\t elec\t%10.4f\n \t\t\t\t\t ee14\t%10.4f\n",
	    (elec+ee14), celec, (elec+ee14)-celec ,elec,ee14 );
    printf ("VDW   %10.4f %10.4f %10.4f \tamber:\n \t\t\t\t\t vdw\t%10.4f\n \t\t\t\t\t vdw14\t%10.4f\n",
	    (vdw+nb14),cvdw,(vdw+nb14)-cvdw, vdw,nb14);
  print "   " 
}
