This is an initial DHFR in TIP3P water periodic test case for md engines supporting
Chamber prmtop and inpcrd files. The Charmm settings and energy are:

NOSHAKE
-------

  UPDATE      bycb           vswi           eps      1.0                -
              INBFRQ    -1   -
              cutnb   11.0   ctofnb    9.0   ctonnb   9.0 -
              Ewald          kappa   0.340  pmEwald        order      4 -
              fftx      96   ffty      80   fftz      64

ENER ENR:  Eval#     ENERgy      Delta-E         GRMS
ENER INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals IMPRopers
ENER CROSS:           CMAPs
ENER EXTERN:        VDWaals         ELEC       HBONds          ASP USER
ENER IMAGES:        IMNBvdw       IMELec       IMHBnd       RXNField EXTElec
ENER EWALD:          EWKSum       EWSElf       EWEXcl       EWQCor EWUTil
 ----------       ---------    ---------    ---------    --------- ---------
ENER>        0-226996.62307      0.00000      2.00904
ENER INTERN>     8578.98727   5018.32063     29.64895    740.94861 14.24181
ENER CROSS>      -216.23917
ENER EXTERN>    30384.03673-241055.92531      0.00000      0.00000 0.00000
ENER IMAGES>     -559.38057  -2893.64396      0.00000      0.00000 0.00000
ENER EWALD>       1122.3929-1202974.7245 1174814.7136       0.0000 0.0000
----------       ---------    ---------    ---------    --------- ---------

PMEMD 2009/10/30 output from step 1 of minimization is (note vdwmeth=0):

   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      1      -2.2699E+05     2.0095E+00     2.5326E+01     O       34449

 BOND    =     8578.9873  ANGLE   =     5018.3206  DIHED      =      740.9486
 UB      =       29.6490  IMP     =       14.2418  CMAP       =     -216.2392
 VDWAALS =    29478.9185  EEL     =  -277456.8839  HBOND      =        0.0000
 1-4 VDW =      345.7376  1-4 EEL =     6475.6373  RESTRAINT  =        0.0000

Summary
Energy Type              CHARMM               AMBER
BOND                 8578.98727           8578.9873   Pass
ANGLE                5018.32063           5018.3206   Pass
DIHED                 740.94861            740.9486   Pass
UB                     29.64895             29.6490   Pass
IMP                    14.24181             14.2418   Pass
CMAP                 -216.23917           -216.2392   Pass
VDW
        30384.03673 (EXT)     29478.9185 (VDWAALS)
         -559.38057 (IMG)       345.7376 (1-4 VDW)
                    29824.65616          29824.6561   Pass

EEL
      -241055.92531 (EXT)   -277456.8839 (EEL)
        -2893.64396 (IMG)      6475.6373 (1-4 EEL)
      -1202974.7245 (EWself)
       1174814.7136 (EWExc)
          1122.3929 (EWKsum)
                   -270987.1874        -270981.2466   Fail at 6 s.f.

SHAKE
-----

  UPDATE      bycb           vswi           eps      1.0                -
              INBFRQ    -1   -
              cutnb   11.0   ctofnb    9.0   ctonnb   9.0 -
              Ewald          kappa   0.340  pmEwald        order      4 -
              fftx      96   ffty      80   fftz      64

  shake fast bonh tol 1.0e-7 para

ENER ENR:  Eval#     ENERgy      Delta-E         GRMS
ENER INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals IMPRopers
ENER CROSS:           CMAPs
ENER EXTERN:        VDWaals         ELEC       HBONds          ASP USER
ENER IMAGES:        IMNBvdw       IMELec       IMHBnd       RXNField EXTElec
ENER EWALD:          EWKSum       EWSElf       EWEXcl       EWQCor EWUTil
 ----------       ---------    ---------    ---------    --------- ---------
ENER>        0-215052.98497      0.00000      0.66823
ENER INTERN>      142.53289    392.90653     29.14707    740.60758 14.18471
ENER CROSS>      -216.41219
ENER EXTERN>    30160.90586-218794.22189      0.00000      0.00000 0.00000
ENER IMAGES>     -560.11010  -2607.54559      0.00000      0.00000 0.00000
ENER EWALD>       1131.9834-1202974.7245 1177487.7612       0.0000 0.0000
 ----------       ---------    ---------    ---------    --------- ---------

PMEMD 2009/10/30 output from step 1 of minimization is (note vdwmeth=0):

   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      1      -2.1495E+05     1.8165E+01     5.5597E+01     O       25755

 BOND    =      139.6651  ANGLE   =      396.1407  DIHED      =      740.9034
 UB      =       30.2509  IMP     =       14.2200  CMAP       =     -216.2637
 VDWAALS =    28948.0282  EEL     =  -251824.9437  HBOND      =        0.0000
 1-4 VDW =      346.1847  1-4 EEL =     6473.3759  RESTRAINT  =        0.0000

Summary
Energy Type              CHARMM               AMBER
BOND                  142.53289            139.6651  Fail
ANGLE                 392.90653            396.1407  Fail
DIHED                 740.60758            740.9034  Fail
UB                     29.14707             30.2509  Fail
IMP                    14.18471             14.2200  Fail
CMAP                 -216.41219           -216.2637  Fail
VDW
        30160.90586 (EXT)     28948.0282 (VDWAALS)
         -560.11010 (IMG)       346.1847 (1-4 VDW)
                    29600.79576          29294.2129   Fail

EEL
      -218794.22189 (EXT)   -251824.9437 (EEL)
        -2607.54559 (IMG)      6473.3759 (1-4 EEL)
      -1202974.7245 (EWself)
       1177487.7612 (EWExc)
          1131.9834 (EWKsum)
                   -245756.747        -245351.5678   Fail

It is clear that CHARMM and PMEMD do shake in different ways so they are, I
believe, operating on different structures here. For the moment we will copy
the pmemd output as the saved test case outputs but this needs investigating. 

