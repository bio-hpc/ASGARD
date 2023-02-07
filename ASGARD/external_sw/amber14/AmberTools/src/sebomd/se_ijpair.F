!------------------------------------------------------------------------
!     this file contains subroutines ijmake, ijfind, and binloc.  they
!     allow fast location of an atom pair (iatm,jatm) in the pairlist
!     ipair.  after the pairlist has been created, ijmake should be
!     called to create three arrays that allow fast access to the pairlist.
!     after calling ijmake, subroutine ijfind may be used to locate the
!     position of (iatm,jatm) in ipair.  binloc is a binary search
!     routine used by ijmake and ijfind.
!
!     name change variables in order to have the same variable names
!     everywhere:: oldname -> newname:
!     .            ipmod1  -> ip1old
!     .            ipack   -> ipold
!     .            ipaddr  -> ijold
!
      subroutine se_ijmake
!
!
!
!     creates three pairlist-size arrays (ip1old,ipold,ijold) that
!     allow near direct access of the address of any atom pair (iatm,jatm)
!     stored in the pairlist array ipair.
!
!     after calling this routine, then given any atom pair (iatm,jatm)
!     that is known to be in the pairlist, the following procedure
!     will locate the address of that pair in the array ipair:
!
!     let npairs = ip1(natoms+1)-1
!         ijpack = (iatm*(iatm-3))/2 + jatm + 1
!          ijmod = mod(ijpack,npairs)
!           imin = ip1old(ijmod)
!           imax = ip1old(ijmod+1)-1
!
!     search ipold(i) over the range i=imin,imax.
!
!     when ipold(i).eq.ijpack then ipmod(i) will contain the position
!     in ipair where (iatm,jatm) is stored.
!
!     the search from imin to imax will span at most nratio entries, where
!     nratio = ((natoms*(natoms-1)/2)/npairs
!            = ratio of potential pairs to actual pairs stored.
!
!     if a binary search is used, then the entry can be located in at
!     most log2(nratio)+1 steps.
!
      implicit none
#include "sebomd.dim"
#include "se_bpairs.h"
#include "se_fulsys.h"
#include "se_pairij.h"

      integer npairs, nmax, istore
      integer ijmod, nwrap, n, ijpack
      double precision discr, xiatm
      integer iatm, jatm, istart, istop, nii, iaddr
!
      npairs = ip1(natoms+1)-1
      nmax = (natoms*(natoms-1))/2 - 1
!
!     loop over the possible values for ijmod.
!
      istore = 1
      do 500 ijmod=1,npairs
!
!       ip1old(ijmod) marks the beginning of the ijmod block in the
!       arrays ipold and ijold.
!
        ip1old(ijmod) = istore
!
!       if all the pairs have been stored then we are done.
!
        if(istore.gt.npairs) return
!
!       determine which combinations of (iatm,jatm) will yield a value
!       of ijpack that gives rise to the current ijmod value.  nwrap
!       keeps track of the number of times the target ijpack value
!       wraps or cycles through npairs.
!
        nwrap = -1
 100    nwrap = nwrap + 1
        n = ijmod - 1 + nwrap*npairs
        if(n.gt.nmax) go to 500
!
!       target ijpack value:
!
        ijpack = n + 1
!
!       must solve the following  quadratic to determine the iatm value
!       that will yield the target ijpack value:
!
!       (iatm*(iatm-3))/2 + 2 = ijpack
!
        discr = 2*n + 0.25d0
        xiatm = dsqrt(discr) + 1.50000001d0
        iatm = int(xiatm)
!
!       skip iatm if it has no pairs in the actual pairlist.
!
        istart = ip1(iatm)
        istop = ip1(iatm+1)-1
        if(istop.lt.istart) go to 100
!
!       determine jatm based on the location of the (iatm,1) pair
!       in the packed storage scheme:
!
        nii = (iatm*(iatm-3))/2
        jatm = n - nii
!
!       do a binary search of the iatm section of the actual pairlist to
!       determine whether the (iatm,jatm) pair has been stored.  if the
!       pair has been stored, then iaddr will contain the corresponding
!       address in ipair.
!
        call se_binloc(ipair,istart,istop,jatm,iaddr)
        if(iaddr.ne.0)then
          ijold(istore) = iaddr
          ipold(istore) = ijpack
          istore = istore + 1
        endif
        go to 100
 500  continue
!
!     overflow pointer in case all possible pairs have been stored:
!
      ip1old(npairs+1) = npairs+1
      return
      end
!
!
!
      subroutine se_ijfind(npairs,iatm,jatm,ijaddr)
!
!     locates the atom pair (iatm,jatm) in the pairlist after the arrays
!     ip1old, ipold, and ijold have been created by subroutine ijmake.
!
!     input:
!
!     npairs = total number of stored pairs.
!
!     iatm,
!     jatm   = atom pair to be located (iatm > jatm).
!
!     returned:
!
!     ijaddr = position in ipair where (iatm,jatm) is stored.
!
!     if (iatm,jatm) is not stored, then ijaddr is returned with
!     a value of zero.
!
      implicit none
#include "sebomd.dim"
#include "se_pairij.h"
      integer npairs, iatm, jatm, ijaddr

      integer ijpack, ijmod, imin, imax, imatch
!
!     index for packed, lower triagular storage:
!
      ijpack = (iatm*(iatm-3))/2 + jatm + 1
!
!     remainder of ijpack/npairs:
!
      ijmod = mod(ijpack,npairs)
      if(ijmod.eq.0) ijmod = npairs
!
!     do a binary search in the array ipold to find a matching value
!     of ijpack.
!
      imin = ip1old(ijmod)
      imax = ip1old(ijmod+1) - 1
      call se_binloc(ipold,imin,imax,ijpack,imatch)
      if(imatch.ne.0)then
        ijaddr = ijold(imatch)
      else
        ijaddr = 0
      endif
      return
      end
!
!
!
      subroutine se_binloc(list,istart,istop,itargt,ifound)
!
!     uses a binary search to locate a target element in a sorted list.
!
!     list   = sorted list of integers to be searched.
!     istart = starting point of search in list.
!     istop  = stopping point of search in list.
!     itargt = value of target element sought in list.
!     ifound = returned location of itarget in list.  if not present, then
!              ifound is returned with a value of 0.
!
      implicit none
      integer istart, istop, itargt, ifound
      integer list(*)

      integer ibegin, iend, ihalf, lhalf
!
      ifound = 0
!
!     use binary search only if list has more than one value.
!
      if(istart.eq.istop)then
        if(list(istart).eq.itargt) ifound = istart
      else
        ibegin = istart
        iend = istop
 10     ihalf = ibegin + (iend-ibegin)/2
        lhalf = list(ihalf)
        if(lhalf.eq.itargt)then
          ifound = ihalf
          return
        elseif(lhalf.gt.itargt)then
          iend = ihalf - 1
        else
          ibegin = ihalf + 1
        endif
        if(ibegin.gt.iend)then
          return
        else
          go to 10
        endif
      endif
      return
      end
