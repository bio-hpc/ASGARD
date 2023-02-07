      subroutine alert(mesg,meslen)
      character mesg*80
      integer meslen
#include "SIZE.h"
#include "TOP.h"
#include "MISC.h"

      write (6,*) mesg(1:meslen)
      stop
      end
