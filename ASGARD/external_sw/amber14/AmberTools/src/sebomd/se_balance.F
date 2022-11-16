      subroutine se_balance_bynorbs
      implicit none
#include "sebomd.dim"
#include "se_parproc.h"
#include "se_bal_mpi.h"
#include "se_subsys.h"

#ifdef MPI 
#include "mpif.h"
      integer target,total(maxproc),min_sum
      integer min_pe,sub_nums(maxSUB),sub_norbs(maxSUB),sum
      integer i, j, norbsk, isub
      logical used

      ! JV sort times into decreasing order to get better distribution
      sum = 0
      my_numsubs = 0
      do i = 1,nproc
        total(i) = 0
      enddo
      do i = 1,nsub
        sub_nums(i) = i
        norbsk = iorbpt(i+1)-iorbpt(i)
        sub_norbs(i) = norbsk
        sum = sum + time_subs(i)
      enddo
      target = sum/nproc
      call se_bsort(nsub,sub_norbs,sub_nums)
      
      ! Distribute 1 subsystem to each PE, biggest times first
      do i = 1,nproc
        isub = nsub-(i-1)
        total(i) = total(i) + sub_norbs(isub)
        ! make sure only pe getting the subsys updates it's list
        if (myid+1 .eq. i ) then
          my_numsubs = my_numsubs + 1
          my_subs(my_numsubs) = sub_nums(isub)
        endif
      enddo

      ! Place remaining subsystems in PE with least increase in total time
      do i = nsub-nproc,1,-1
        ! Place next subsystem in bin (PE) with smallest total time so far
        min_sum = total(1)
        min_pe = 1
        do j = 2,nproc
          if( total(j) .lt. min_sum) then
            min_sum = total(j)
            min_pe = j
          endif
        enddo
        total(min_pe) = total(min_pe) + sub_norbs(i)
        if (myid+1 .eq. min_pe ) then
          my_numsubs = my_numsubs + 1
          my_subs(my_numsubs) = sub_nums(i)
        endif
      enddo
      return
#endif
      end

      subroutine se_balance_bytime
      implicit none
#include "sebomd.dim"
#include "se_parproc.h"
#include "se_bal_mpi.h"
#include "se_subsys.h"

#ifdef MPI 
      double precision  target,total(maxproc),min_sum,sum
      integer min_pe,sub_nums(maxSUB)
      integer i, j, isub

      ! JV sort times into decreasing order to get better distribution
      sum = 0.0
      my_numsubs = 0
      do i = 1,nproc
        total(i) = 0.0
      enddo
      do i = 1,nsub
        sub_nums(i) = i
        sum = sum + time_subs(i)
      enddo
      target = sum/nproc

      call se_rsort(nsub,time_subs,sub_nums)

      ! Distribute 1 subsystem to each PE, biggest times first
      do i = 1,nproc
        isub = nsub-(i-1)
        total(i) = total(i) + time_subs(isub)
        ! make sure only pe getting the subsys updates it's list
        if (myid+1 .eq. i ) then
          my_numsubs = my_numsubs + 1
          my_subs(my_numsubs) = sub_nums(isub)
        endif
      enddo

      ! Place remaining subsystems in PE with least increase in total time
      do i = nsub-nproc,1,-1
      ! Place next subsystem in bin (PE) with smallest total time so far
        min_sum = total(1)
        min_pe = 1
        do j = 2,nproc
          if( total(j) .lt. min_sum) then
            min_sum = total(j)
            min_pe = j
          endif
        enddo
        total(min_pe) = total(min_pe) + time_subs(i)
        if (myid+1 .eq. min_pe ) then
          my_numsubs = my_numsubs + 1
          my_subs(my_numsubs) = sub_nums(i)
        endif
      enddo
      return
#endif
      end
