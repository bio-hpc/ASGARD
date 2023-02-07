#include "copyright.h"
#include "../include/dprec.fh"
#include "ncsu-config.h"
#include "../include/assert.fh"
#ifndef PBSA
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Open input files and read cntrl namelist.
subroutine mdread1()
#undef API
#include "mdread1.F90"
end subroutine mdread1 
#endif /*ifndef PBSA*/

#ifndef PBSA
!======================================================================
!          MDREAD2
!======================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initialize to defaults and print the inputable variables.
subroutine mdread2(x,ix,ih)
#undef API
#include "mdread2.F90"
end subroutine mdread2 
#endif /* ifndef PBSA */


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit defined preprocessor names, ie, flags.
#include "mdreadutils.h"
