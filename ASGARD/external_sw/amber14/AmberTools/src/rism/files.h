!+ Specification and control of RISM Input/Output

character(len=256) xvvfil, guvfil, huvfil, cuvfil, &
      rismcrdfil, rismfrcfil, rismcrdrstfil, rismfrcrstfil
common /filesr/ xvvfil, guvfil, huvfil, cuvfil, &
      rismcrdfil, rismfrcfil, rismcrdrstfil, rismfrcrstfil

integer     RISMCRD_UNIT
integer     RISMFRC_UNIT
integer     RISMCRDRST_UNIT
integer     RISMFRCRST_UNIT
parameter ( RISMCRD_UNIT = 30 )
parameter ( RISMFRC_UNIT = 31 )
parameter ( RISMCRDRST_UNIT = 32 )
parameter ( RISMFRCRST_UNIT = 33 )
