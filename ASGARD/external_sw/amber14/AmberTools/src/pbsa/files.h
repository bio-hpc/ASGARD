!+ Specification and control of Amber's Input/Output
!+ R. Luo added pqr support

! File names
character(len=512) groupbuffer
character(len=80) mdin, mdout, inpcrd, parm, restrt, &
      refc, mdvel, mden, mdcrd, mdinfo, nmr, mincor, &
      vecs, radii, freqe,redir(8), &
      rstdip,mddip,inpdip,groups,gpes, &
      cpin, cpout, cprestrt, &
      pqr

character owrite
common /files/ groupbuffer, mdin, mdout, inpcrd, parm, restrt, &
      refc, mdvel, mden, mdcrd, mdinfo, nmr, mincor, &
      vecs, radii, freqe, owrite, &
      rstdip,mddip,inpdip,groups,gpes, &
      cpin, cpout, cprestrt, &
      pqr

! put this in a separate common block to stop the compiler from
! complaining about misalignment
integer numgroup
common/nmgrp/ numgroup


! File units
! An I/O Unit resource manager does not exist.
integer     MDCRD_UNIT
integer     MDEN_UNIT
integer     MDINFO_UNIT
integer     MDVEL_UNIT
parameter ( MDINFO_UNIT =  7 )
parameter ( MDCRD_UNIT  = 12 )
parameter ( MDEN_UNIT   = 15 )
parameter ( MDVEL_UNIT  = 13 )
integer, parameter :: CNSTPH_UNIT = 18, CPOUT_UNIT = 19
integer, parameter :: PQR_UNIT = 20


! File related controls and options
character(len=80) title,title1
common/runhed/ title, title1

logical mdin_pb
common/mdin_flags/mdin_pb

integer BC_HULP  ! size in integers of common HULP
parameter ( BC_HULP = 9 )

integer     ntpr,ntwr,ntwx,ntwv,ntwe,ntpp,ioutfm,ntwprt,ntave
common/hulp/ntpr,ntwr,ntwx,ntwv,ntwe,ntpp,ioutfm,ntwprt,ntave

!      NMRRDR : Contains information about input/output file redirection
!               REDIR and IREDIR contain information regarding
!               LISTIN, LISTOUT, READNMR, NOESY, SHIFTS, DUMPAVE,
!               PCSHIFT and DIPOLE respectively. If IREDIR(I) > 0,
!               then that input/output has been redirected.

integer iredir(8)
common/nmrrdr/redir,iredir
