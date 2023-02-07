      subroutine se_debugmsg0(message)
      implicit none
#include "se_logic.h"
      character(*) message

      if (master.and.debug) then
        write(0,'(a)') message
      endif
      return
      end

      subroutine se_debugmsg1i(message, a)
      implicit none
#include "se_logic.h"
      character(*) message
      integer a 

      if (master.and.debug) then
        write(0,'(a,i6)') message, a
      endif
      return
      end

      subroutine se_debugmsg3i(message, a1,a2,a3)
      implicit none
#include "se_logic.h"
      character(*) message
      integer a1,a2,a3

      if (master.and.debug) then
        write(0,'(a,3i6)') message, a1,a2,a3
      endif
      return
      end

      subroutine se_debugmsg1f(message, a)
      implicit none
#include "se_logic.h"
      character(*) message
      double precision a 

      if (master.and.debug) then
        write(0,'(a,f25.16)') message, a
      endif
      return
      end

      subroutine se_debugmsg2f(message, a1,a2)
      implicit none
#include "se_logic.h"
      character(*) message
      double precision a1,a2

      if (master.and.debug) then
        write(0,'(a,2f25.16)') message, a1,a2
      endif
      return
      end

      subroutine se_debugmsg3f(message, a1,a2,a3)
      implicit none
#include "se_logic.h"
      character(*) message
      double precision a1,a2,a3

      if (master.and.debug) then
        write(0,'(a,3f25.16)') message, a1,a2,a3
      endif
      return
      end

