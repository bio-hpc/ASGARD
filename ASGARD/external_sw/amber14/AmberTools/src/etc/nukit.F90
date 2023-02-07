program nukit
   
   ! PROGRAM TO INTERROGATE USER & GENERATE INPUT FILES
   ! NUC.DAT, LIN.DAT FOR NUCGEN & LINK ROUTINES 
   !
   ! BILL ROSS, TINOCO LAB, BERKELEY,1987
   
   implicit none
   
   character(len=80) :: inline, outline
   character(len=3) :: res
   character :: ch, flag
   character, parameter :: blank = ' '
   
   integer :: i, j, k, inew
   
   write(6,*) 'Residue naming convention? (O = pre-94, N = 94) O/N:'
   read(*,200) inline
   
   if (inline .eq. 'N') then
      inew = 1
   else if (inline .eq. 'O') then
      inew = 0
   else
      write(6,*) 'unknown option: ', inline
      stop
   end if
   
   write(6,420)
   read(*,200) inline
   
   open(unit=9,file='lin.in',status='NEW')
   open(unit=8,file='nuc.in',status='NEW')
   write(9,200) inline
   write(9,*) ' '
   write(9,500)' '
   outline(1:80) = 'DU' 
   write(9,200) outline
   outline(1:80) = '     0    0    0    0    0  '
   write(9,200) outline
   write(6,700)
   write(6,800)
   !write(6,810)
   
   ! -- do 2 strands
   do i=1,2
      write(9,400) i
      write(8,400) i
      
   10 inline(1:80) = ' '
      write(6,900)
      read(*,200) inline
      outline(1:80) = '   '
      
      if (inew .eq. 0) then
         outline(1:80) = 'HB '
         k = 1
      else
         k = 0
      end if
      
      write(6,920)
      read(*,500) flag
      
      if (i .eq. 1) then
         write(9,500) flag
      else
         write(9,500) flag, 1, 3
      end if
      
      write(8,500) flag
      j = 1
      
      ! primitive while() loop over residues in strand:
      
   20 if (inline(j:j) .eq. blank) goto 30
      
      if (k .gt. 10) then
         write(8,200) outline
         outline(5:5) = '2'
         write(9,200) outline
         outline(1:80) = ' '
         k = 0
      end if
      
      ch = inline(j:j)
      
      if (inew .eq. 0) then
         if (ch .eq. 'C') then
            res = 'CYT'
         else if (ch .eq. 'G') then
            res = 'GUA'
         else if (ch .eq. 'T') then
            res = 'THY'
         else if (ch .eq. 'U') then
            res = 'URA'
         else if (ch .eq. 'A') then        
            res = 'ADE'
         else 
            write(6,600) j
            goto 10
         end if
      else
         if (ch .eq. 'C') then
            res = 'C  '
         else if (ch .eq. 'G') then
            res = 'G  '
         else if (ch .eq. 'T') then
            res = 'T  '
         else if (ch .eq. 'U') then
            res = 'U  '
         else if (ch .eq. 'A') then
            res = 'A  '
         else
            write(6,600) j
            goto 10
         end if
         if (j .eq. 1) res(2:2) = '5'
      end if
      j = j + 1
      outline((5*k+1):(5*k+3)) = res
      k = k + 1
      if (inew .eq. 0 .and. inline(j:j) .ne. blank) then
         outline((5*k+1):(5*k+3)) = 'POM'
         k = k + 1
      end if
      goto 20
      
   30 continue
      if (inew .eq. 0) then
         outline((5*k+1):(5*k+2)) = 'HE'
      else
         outline((5*(k-1)+2):(5*(k-1)+2)) = '3'
      end if
      
      write(8,200) outline
      outline(5:5) = '2'
      write(9,200) outline
      outline(1:80) = ' '
      write(8,200) outline
      write(9,200) outline
   end do
   
   outline(1:4) = 'QUIT'
   write(9,200) outline
   outline(1:4) = 'END '
   write(8,200) outline
   write(6,901)
   write(6,902)
   write(6,903)
   write(6,904)
   write(6,905)
   write(6,906)
   write(6,907)
   write(6,908)
   write(6,909)
   write(6,930)
   write(6,932)
   write(6,910)
   inline(1:80) = ' '
   read(5,200) inline
   write(8,300) inline
   close(9)
   close(8)
   
   200 format(a80)
   300 format(t1, '$',a8)
   400 format('  NUC ',i3,i8,i4)
   420 format(1x,'JOB NAME?  ')
   500 format(a,i14,i5)
   600 format('***COULDN''T MATCH BASE ',i2,' TRY AGAIN:'/ /)
   700 format('--------(from here on, USE CAPITALS)       ---------')
   800 format('---------e.g. CGCGATAT                     ---------')
   810 format('--------(Note: "8" =brom8gua, "5" =brom5cyt)--------')
   900 format(/ /1x,'ENTER SEQUENCE {5-prime to 3-prime}:   ')
   901 format('             CONFORMATIONS:')
   902 format('     ARNA   right handed A RNA (arnott)')
   903 format('     APRNA  right handed A-prime RNA (arnott)')
   904 format('     LBDNA  right handed B DNA (langridge)')
   905 format('     ABDNA  right handed B DNA (arnott)')
   906 format('     SBDNA  left handed  B DNA (sasisekharan)')
   907 format('     ADNA   right handed A DNA (arnott)')
   908 format('     NIXON  none of above - nucgen.pdb user-supplied')
   909 format(' ')
   930 format('  A-forms may need work to place H1-primes properly. ')
   932 format(' ')
   910 format(1x,'CONFORMATION?  ')
   920 format(1x,'DNA OR RNA? (D/R): ')
   
end program nukit
