      BLOCK DATA
      IMPLICIT none
#include "sebomd.dim"

#include "se_param1.h"
#include "se_param2.h"
#include "se_param3.h"
#include "se_param4.h"
#include "se_param5.h"
#include "se_param6.h"
#include "se_param7.h"

#include "se_block_mndo.h"
#include "se_block_am1.h"
#include "se_block_pm3.h"
#include "se_block_rm1.h"
#include "se_block_am1dphot.h"
#include "se_block_pm3pddg.h"

#include "se_logic.h"
#include "se_subpar.h"
#include "se_pbcbox.h"
#include "se_elemnt.h"
#include "se_pmeee.h"
#include "se_parproc.h"
#include "se_timings.h"
#include "se_crit.h"
#include "se_scfstuff.h"
#include "se_endecomp.h"
#include "se_bpairs.h"
#include "se_fopen.h"
#include "se_funit.h"
C
C
C     ASSIGN NONPARAMETRIC ELEMENTAL QUANTITIES.
C
C     ELEMENTAL SYMBOL:
C
      DATA SYMBOL
     .       /'XX','H ','HE','LI','BE','B ','C ','N ','O ','F ','NE',
     .             'NA','MG','AL','SI','P ','S ','CL','AR','K ','CA',
     .             'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     .             'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR',
     .             'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN',
     .             'SB','TE','I ','XE','CS','BA','LA','CE','PR','ND',
     .             'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB',
     .             'LU','HF','TA','W ','RE','OS','IR','PT','AU','HG',
     .             'TL','PB','BI'/

C     atomic masses
      data emass
     &     /0.0, 1.0079, 4.0026, 6.941, 9.01218, 10.81, 12.011,
     &     14.0067, 15.99994, 18.99840, 20.179, 22.9898, 24.305,
     &     26.98154, 28.0855, 30.97376, 32.06, 35.453, 39.948,
     &     39.0983, 40.08, 44.9559, 47.90, 50.9415, 51.996, 54.938,
     &     55.847, 58.9332, 58.71, 63.546, 65.38, 69.737, 72.59,
     &     74.9216, 78.96, 79.904, 83.80, 85.4678, 87.62, 88.9059,
     &     91.22, 92.9064, 95.94, 98.9062, 101.07, 102.9055, 106.4,
     &     107.868, 112.41, 114.82, 118.69, 121.75, 127.60, 126.9045,
     &     131.30, 132.9054, 137.33, 15*0.000, 178.49, 180.9479,
     &     183.850, 186.207, 190.20, 192.220, 195.090,
     &     196.9665, 200.590, 204.370, 207.200, 208.9804/

C
C     CORE NUCLEAR CHARGES:
C
      DATA ZCHG
     .      /0.0, 1.0, 2.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 0.0,
     .            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 0.0, 1.0, 2.0,
     .            3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,10.0,11.0, 2.0,
     .            3.0, 4.0, 5.0, 6.0, 7.0, 0.0, 1.0, 2.0, 3.0, 4.0,
     .            5.0, 6.0, 7.0, 8.0, 9.0,10.0,11.0, 2.0, 3.0, 4.0,
     .            5.0, 6.0, 7.0, 0.0,
     .            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,10.0,
     .            11.0,12.0,13.0,14.0,15.0,16.0,3.0, 4.0, 5.0, 6.0,
     .            7.0, 8.0, 9.0,10.0,11.0, 2.0, 3.0, 4.0, 5.0/

C
C     NUMBER OF VALENCE ATOMIC ORBITALS:
C     (K IS SPARKLE -- NO ORBITALS)
C
! 2007-06-08: put orbitals and 1 electron on Na (Ed's parameters) (GM-JT)
      DATA NATORB
     .       / 0,   1,   0,   4,   4,   4,   4,   4,   4,   4,   0,
     .              4,   4,   4,   4,   4,   4,   4,   0,   0,   4,
     .              4,   4,   4,   4,   4,   4,   4,   4,   4,   4,
     .              4,   4,   4,   4,   4,   0,   4,   4,   4,   4,
     .              4,   4,   4,   4,   4,   4,   4,   4,   4,   4,
     .              4,   4,   4,   0,
     .              4,   4,   4,   4,   4,   4,   4,   4,   4,   4,
     .              4,   4,   4,   4,   4,   4,   4,   4,   4,   4,
     .              4,   4,   4,   4,   4,   4,   4,   4,   4/

C
C     NUMBER OF VALENCE ELECTRONS:
C
      DATA NVELEC
     .       / 0,   1,   0,   1,   2,   3,   4,   5,   6,   7,   0,
     .              1,   2,   3,   4,   5,   6,   7,   0,   0,   2,
     .              3,   4,   5,   6,   7,   8,   9,  10,  11,   2,
     .              3,   4,   5,   6,   7,   0,   1,   2,   3,   4,
     .              5,   6,   7,   8,   9,  10,  11,   2,   3,   4,
     .              5,   6,   7,   0,
     .              1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
     .             11,  12,  13,  14,  15,  16,   3,   4,   5,   6,
     .              7,   8,   9,  10,  11,   2,   3,   4,   5/
C
C     PRINCIPAL QUANTUM NUMBER:
C
      DATA NQUANT
     .       / 0,   1,   1,   2,   2,   2,   2,   2,   2,   2,   2,
     .              3,   3,   3,   3,   3,   3,   3,   3,   4,   4,
     .              4,   4,   4,   4,   4,   4,   4,   4,   4,   4,
     .              4,   4,   4,   4,   4,   4,   5,   5,   5,   5,
     .              5,   5,   5,   5,   5,   5,   5,   5,   5,   5,
     .              5,   5,   5,   5,
     .              6,   6,   6,   6,   6,   6,   6,   6,   6,   6,
     .              6,   6,   6,   6,   6,   6,   6,   6,   6,   6,
     .              6,   6,   6,   6,   6,   6,   6,   6,   6/

C
C
C******************************************************************************
C
C     ASSIGN UNIT NUMBERS, FILENAMES, AND STATUS FOR FILES THAT MAY BE OPENED.
C     note that fname is the basename for the filename for parallel I/O
C     ext is the extension used for parallel I/O

C
C     A. FILES USED FOR PARALLEL I/O:
C     ------------------------------
C
C
C     A.1. FILES FOR THE MAIN NODE (MYID=0)
C     .    (so no EXT is needed):
C
C     MAIN INPUT FILE CONTAINING KEYWORDS, TITLE, COORDINATES, ETC.:
C
C     density matrix:
      data iunit(1) /27/
      data fname(1) /'sebomd.dmx          '/
      data fstat(1) /'unknown'/
C
C
C     B. FILES FOR SERIAL I/O:
C     -----------------------
C
C     file with guess density matrix (matrices)
      data iunit(2) /39/
      data fname(2) /'guess.dmx           '/
      data fstat(2) /'old'/

C     data for /subpar/
      data dcbuff1 /0.0/
      data dcbuff2 /0.0/
      data ncores /0/

C     energy decomposition stuff
      data noovrlp /.false./

C     ENERGY (eV) AND DENSITY MATRIX CONVERGENCE CRITERIA:
      data eecrt /1.0D-12/
      data dencrt /5.0D-4/
      data iprec/3/

C     box-dimensions: initialize them to be big, so you won't
C     get an error message if you don't do a pbc calculation, but
C     still want to do grid-subsetting.
      data dbox(1) /9.0D9/
      data dbox(2) /9.0D9/
      data dbox(3) /9.0D9/

C     data for parallel processing with MPI
      data myid /0/
      data nproc /1/

C     logicals /logic/

      data pbc /.false./
      data prtsub /.false./
      data screen /.false./
      data stand /.false./
      data fullscf /.false./
      data clust /.false./
      data autosub /.false./
      data am1 /.false./
      data pm3 /.false./
      data mndo /.false./
      data rm1 /.false./
      data pm3pddg /.false./
      data cart /.true./
      data gradient /.false./
      data pme /.false./
      data gridsub /.false./
      data inter /.false./
      data recipintr /.false./
      data pole /.false./
      data ionpot /.false./
      data wrtscr /.false./
      data guess /.false./
!     data ipolyn /.false./
      data wrtedmx /.true./
      data pdb/.false./

C     data for bpairs
      data bcut2 /0.0/
      data bcut /0.0/

C     data for timings
      data ndiag /0/

      data maxitscf /1000/

      END
