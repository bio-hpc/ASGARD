      subroutine se_bsort(n,ia,ib)
!
!     Heap sort routine for list of integer point numbers ia.
!     The integer values ib are carried along.  Sorts in order of
!     increasing value.
!
!     heap sort algorithm instead of bubble sort:
!     bubble sort = N*N, heap sort = N Log2(N)

      implicit none
    
      integer n
      integer ia(*), ib(*)

      integer l, ir, iia, iib, i, j

      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          iia=ia(l)
          iib=ib(l)
        else
          iia=ia(ir)
          iib=ib(ir)
          ia(ir)=ia(1)
          ib(ir)=ib(1)
          ir=ir-1
          if(ir.eq.1)then
            ia(1)=iia
            ib(1)=iib
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ia(j).lt.ia(j+1))j=j+1
          endif
          if(iia.lt.ia(j))then
            ia(i)=ia(j)
            ib(i)=ib(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ia(i)=iia
        ib(i)=iib
      goto 10
      end
!---------------------------------------------------------------------
      subroutine se_bovsort1(n2,ia,ierror)
!
!     returns list of unique and non-overlapping integers
!     used in rddmx.F to significantly decrease copy time
!     of density matrix later on
!
!     note that ia(i) contains the starting index,
!     .         ia(i+1) the ending index, i is odd
!
!     so array 5 7 1 3 8 12 13 13 will be overwritten by 1 3 5 13
!
!     ierror will be set to one in case there were overlapping /
!     nonuniqe numbers, else ierror = 0
!     AvdV

      implicit none

      integer n2, ierror, ia(*)
      integer n, l, l2, l1, ir, ir2, ir1, iia1, iia2
      integer i, j, i2, i1, j1, j2, k, k1, k2, ia1, ia2
      integer ileft1, iright1
      integer ileft2, iright2

      ierror = 0

!     A. sort stuff

      n = n2/2
      if (n.lt.2) return
      l=n/2+1
      l2=2*l
      l1=l2-1
      ir=n
      ir2=2*ir
      ir1=ir2-1
10    continue
        if(l.gt.1)then
          l=l-1
          l2=2*l
          l1=l2-1
          iia1=ia(l1)
          iia2=ia(l2)
        else
          iia1=ia(ir1)
          iia2=ia(ir2)
          ia(ir1)=ia(1)
          ia(ir2)=ia(2)
          ir=ir-1
          ir2=2*ir
          ir1=ir2-1
          if(ir.eq.1)then
            ia(1)=iia1
            ia(2)=iia2
            goto 100
          endif
        endif
        i=l
        i2=2*i
        i1=i2-1
        j=l+l
        j2=2*j
        j1=j2-1
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ia(j1).lt.ia(j1+2)) then
               j=j+1
               j2=2*j
               j1=j2-1
            endif
          endif
          if(iia1.lt.ia(j1))then
            ia(i1)=ia(j1)
            ia(i2)=ia(j2)
            i=j
            i2=2*i
            i1=i2-1
            j=j+j
            j2=2*j
            j1=j2-1
          else
            j=ir+1
            j2=2*j
            j1=j2-1
          endif
        goto 20
        endif
        ia(i1)=iia1
        ia(i2)=iia2
      goto 10

!     B. get unique pairs

 100  i=1
      i2=2*i
      i1=i2-1
      k=i
      k2=i2
      k1=i1
      ia1=ia(i1)
      ia2=ia(i2)

 110  i=i+1
      if (i.gt.n) goto 200
      i2=2*i
      i1=i2-1
      if (ia(i1).eq.ia1) then
         ierror = 1
         ia2 = max(ia2,ia(i2))
      else
         ia(k1) = ia1
         ia(k2) = ia2
         k = k+1
         k2=2*k
         k1=k2-1
         ia1=ia(i1)
         ia2=ia(i2)
      endif
      goto 110

 200  ia(k1) = ia1
      ia(k2) = ia2
      n = k

!     C. get overlap

      i=1
      i2=2*i
      i1=i2-1
      k=2
      k2=2*k
      k1=k2-1
      ileft1 = ia(1)
      iright1 = ia(2)

 310  if (k.le.n) then
         ileft2 = ia(k1)
         iright2 = ia(k2)
         if (iright1.ge.(ileft2-1)) then
            if (iright1.ge.ileft2) ierror = 1
            if (iright1.lt.iright2) then
               iright1 = iright2
            endif
         else
            ia(i1) = ileft1
            ia(i2) = iright1
            i=i+1
            i2=2*i
            i1=i2-1
            ileft1 = ia(k1)
            iright1 = ia(k2)
         endif
         k=k+1
         k2=2*k
         k1=k2-1
         goto 310
      endif
      ia(i1) = ileft1
      ia(i2) = iright1
      n2 = i2

      end
!---------------------------------------------------------------------
      subroutine se_rsort(n,ra,ia)
!
!     Heap sort routine for list of floating point numbers ra.
!     The integer values ia are carried along.  Sorts in order of
!     increasing value.
!
      implicit none

      integer n, ia(*)
      double precision ra(*)

      integer ir, l, iia, i, j
      double precision rra

      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1.0)then
          l=l-1
          rra=ra(l)
          iia=ia(l)
        else
          rra=ra(ir)
          iia=ia(ir)
          ra(ir)=ra(1)
          ia(ir)=ia(1)
          ir=ir-1
          if(ir.eq.1.0)then
            ra(1)=rra
            ia(1)=iia
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            ia(i)=ia(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ra(i)=rra
        ia(i)=iia
      goto 10
      end
