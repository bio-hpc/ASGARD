! PUPIL program
! 
! Original Author:       J. Torras
!                        joan.torras@upc.edu
! Last Modified By:      B. P. Roberts
!                        roberts@qtp.ufl.edu
! Version:               2.0
! Date:                  2011-10-10
!=======================================================================

      subroutine fixport
      
      character(len=30) :: argv, &
                           host, &
                           jxms, jxmx, jxss, &
                           optprint, &
                           port 
      integer :: i, iargc, ijxms, ijxmx, ijxss, ioptprint, &
                 m
      
      m = iargc()
      
      do i=1,m
         call getarg(i,argv)
         if (argv .eq. "-ORBInitialPort") then
            call getarg(i+1,port)
            write(*,*) 'Port: ', trim(port)
         else if (argv .eq. '-ORBInitialHost') then
            call getarg(i+1,host)
            write(*,*) 'Host: ', trim(host)
         else if (argv .eq. '-OptPrint') then
            call getarg(i+1,optprint)
            read(optprint,'(I30)') ioptprint
            write(*,*) 'Log level: ', ioptprint
         else if (argv .eq. '-jxms') then
            call getarg(i+1,jxms)
            read(jxms,'(I30)') ijxms
            write(*,*) 'Minimum Java memory: ', ijxms
         else if (argv .eq. '-jxmx') then
            call getarg(i+1,jxmx)
            read(jxmx,'(I30)') ijxmx
            write(*,*) 'Maximum Java memory: ', ijxmx
         else if (argv .eq. '-jxss') then
            call getarg(i+1,jxss)
            read(jxss,'(I30)') ijxss
            write(*,*) 'Stack size of Java: ', ijxss
         end if
      end do
      
      ! Assign host, port and log level to PUPIL system
      call setcorbanameserver(host,port,ioptprint)
      call setjavamemory(ijxms,ijxmx,ijxss)
      
      return
      
      end subroutine fixport
     
!=======================================================================
