
!************************************************************************
!                              AMBER                                   **
!                                                                      **
!                  Copyright (c) 1986, 1991, 1995                      **
!             Regents of the University of California                  **
!                       All Rights Reserved.                           **
!                                                                      **
!  This software provided pursuant to a license agreement containing   **
!  restrictions on its disclosure, duplication, and use. This software **
!  contains confidential and proprietary information, and may not be   **
!  extracted or distributed, in whole or in part, for any purpose      **
!  whatsoever, without the express written permission of the authors.  **
!  This notice, and the associated author list, must be attached to    **
!  all copies, or extracts, of this software. Any additional           **
!  restrictions set forth in the license agreement also apply to this  **
!  software.                                                           **
!************************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rdparm here]
subroutine rdparm(x, ix, ih, cg, maxdup)
   
   !     --- finishes reading the prmtop file, which was begun in nmode()
   
   implicit double precision (a-h, o-z)
   
   dimension x(*)
   dimension ix(*)
   character(len=4) ih(*)
   dimension cg(*)
   
#  include "pointer.h"
#  include "inpdat.h"
#  include "pol.h"
   
   character(len=80) fmt
   character(len=80) fmtin,ifmt,afmt,rfmt,type
   ifmt = '(12I6)'
   afmt = '(20A4)'
   rfmt = '(5E16.8)'
   nf = 51
   
   ntype = ntypes * ntypes
   nttyp = ntypes * (ntypes+1) / 2
   
   !     ----- read atomic and residue properties of the molecule(s)
   
   fmtin = afmt
   type = 'ATOM_NAME'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read (51,fmt) ( ih(i), i = mgraph, mgraph + natom - 1 )
   
   fmtin = rfmt
   type = 'CHARGE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=mchrg,mchrg+natom-1)
   
   fmtin = rfmt
   type = 'MASS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=mamass,mamass+natom-1)
   
   fmtin = ifmt
   type = 'ATOM_TYPE_INDEX'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (ix(i),i=miac,miac+natom-1)
   
   fmtin = ifmt
   type = 'NUMBER_EXCLUDED_ATOMS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (ix(i),i=miblo,miblo+natom-1)
   
   fmtin = ifmt
   type = 'NONBONDED_PARM_INDEX'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (ix(i),i=mico,mico+ntype-1)
   
   fmtin = afmt
   type = 'RESIDUE_LABEL'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read (51,fmt) ( ih(i), i = mlbres, mlbres + nres - 1 )
   
   fmtin = ifmt
   type = 'RESIDUE_POINTER'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (ix(i),i=mipres,mipres+nres-1)
   
   !     ----- the last ipres array element should indicate the
   !     ----- beginning of a residue beyond the last one as if
   !     ----- it were present
   
   ix(mipres + nres) = natom + 1
   
   !     ----- read parameters such as bond lengths, force constants
   
   fmtin = rfmt
   type = 'BOND_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=mrk,mrk+numbnd-1)
   
   fmtin = rfmt
   type = 'BOND_EQUIL_VALUE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=mreq,mreq+numbnd-1)
   
   fmtin = rfmt
   type = 'ANGLE_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=mtk,mtk+numang-1)
   
   fmtin = rfmt
   type = 'ANGLE_EQUIL_VALUE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=mteq,mteq+numang-1)
   
   fmtin = rfmt
   type = 'DIHEDRAL_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=mpk,mpk+nptra-1)
   
   fmtin = rfmt
   type = 'DIHEDRAL_PERIODICITY'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=mpn,mpn+nptra-1)
   
   fmtin = rfmt
   type = 'DIHEDRAL_PHASE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=mphase,mphase+nptra-1)
   
   fmtin = rfmt
   type = 'SOLTY'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=msolty,msolty+natyp-1)
   
   fmtin = rfmt
   type = 'LENNARD_JONES_ACOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=mcn1,mcn1+nttyp-1)
   
   fmtin = rfmt
   type = 'LENNARD_JONES_BCOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=mcn2,mcn2+nttyp-1)
   
   !     ----- read bonding information
   
   fmtin = ifmt
   type = 'BONDS_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt)      (ix(i+ mibh-1), &
         ix(i+ mjbh-1), &
         ix(i+micbh-1), &
         i=1,nbonh)
   
   fmtin = ifmt
   type = 'BONDS_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt)      (ix(i+ miba-1), &
         ix(i+ mjba-1), &
         ix(i+micba-1), &
         i=1,nbona)
   
   fmtin = ifmt
   type = 'ANGLES_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt)      (ix(i+ mith-1), &
         ix(i+ mjth-1), &
         ix(i+ mkth-1), &
         ix(i+micth-1), &
         i=1,ntheth)
   
   fmtin = ifmt
   type = 'ANGLES_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt)      (ix(i+ mita-1), &
         ix(i+ mjta-1), &
         ix(i+ mkta-1), &
         ix(i+micta-1), &
         i=1,ntheta)
   
   fmtin = ifmt
   type = 'DIHEDRALS_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt)      (ix(i+ miph-1), &
         ix(i+ mjph-1), &
         ix(i+ mkph-1), &
         ix(i+ mlph-1), &
         ix(i+micph-1), &
         i=1,nphih)
   
   fmtin = ifmt
   type = 'DIHEDRALS_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt)      (ix(i+ mipa-1), &
         ix(i+ mjpa-1), &
         ix(i+ mkpa-1), &
         ix(i+ mlpa-1), &
         ix(i+micpa-1), &
         i=1,nphia)
   
   fmtin = ifmt
   type = 'EXCLUDED_ATOMS_LIST'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (ix(i+minb-1),i=1,nnb)
   
   !     ----- read H-bond parameters
   
   fmtin = rfmt
   type = 'HBOND_ACOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=masol,masol+nphb-1)
   
   fmtin = rfmt
   type = 'HBOND_BCOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=mbsol,mbsol+nphb-1)
   
   fmtin = rfmt
   type = 'HBCUT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (x(i),i=mhbcut,mhbcut+nphb-1)
   
   !     ----- read the symbol, tree, join, irotat arrays:
   
   fmtin = afmt
   type = 'AMBER_ATOM_TYPE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read (51,fmt) ( ih(i), i = msymbl, msymbl + natom - 1 )
   
   fmtin = afmt
   type = 'TREE_CHAIN_CLASSIFICATION'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read (51,fmt) ( ih(i), i = mitree , mitree  + natom - 1 )
   
   fmtin = ifmt
   type = 'JOIN_ARRAY'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (ix(i),i=mjoin,mjoin+natom-1)
   
   fmtin = ifmt
   type = 'IROTAT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(51,fmt) (ix(i),i=mrotat,mrotat+natom-1)
   
   !  read in polariziabilites
   
   if(ipol /= 0) then
      fmtin = rfmt
      type = 'POLARIZABILITY'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(51,fmt) (x(i),i=mpol,mpol+natom-1)
   end if
   
   call pload(natom,ntypes,ix(miac),ix(mico))

   if(ipol /= 0) then
      ip14 = 0
      call load14(nphih,ix(miph),ix(mjph),ix(mkph),ix(mlph), &
            ix(micph),x(mpn),n14,ni14,ip)
      ip14 = ip14 + ip
      call load14(mphia,ix(mipa),ix(mjpa),ix(mkpa),ix(mlpa), &
            ix(micpa),x(mpn),n14,ni14,ip)
      ip14 = ip14 + ip
      
      write(6,3377)ip14
      3377 format(t2,'%MINMD/POL-I-14Pairs, ',i4,' pairs for pol 1-4.'/)
   end if
   
   !   ---duplicate the multiple-term torsions:
   
   idum = nphih
   call dihdup(nphih,idum,ix(miph),ix(mjph),ix(mkph),ix(mlph), &
         ix(micph),x(mpn),maxdup)
   call dihdup(mphia,nphia,ix(mipa),ix(mjpa),ix(mkpa),ix(mlpa), &
         ix(micpa),x(mpn),maxdup)

   
   !     ----- Scale the charges by the dielectric constatnt
   
   if (dielc /= 1.0) then
      sqdiel  = sqrt (dielc)
      do 10 i = 1,natom
         cg(i) = cg(i) / sqdiel
      10 continue
   end if
   
#ifdef CHARMM
   
   !   ---- read extra 1-4 nonbonded and urey-bradley parms
   
   read(51,'(5e16.8)') (x(i),i=mcn114,mcn114+nttyp-1)
   read(51,'(5e16.8)') (x(i),i=mcn214,mcn214+nttyp-1)
   
   read(51,'(5e16.8)') (x(i),i=mfk ,mfk +numang-1)
   read(51,'(5e16.8)') (x(i),i=mqeq,mqeq+numang-1)
#endif
   return
   1 format (12i6)
   2 format (20a4)
end subroutine rdparm 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine istuff here]
subroutine istuff(i,j,iarray,k)
   dimension iarray(15,*)
   iarray(i,j) = k
   return
end subroutine istuff 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine load14 here]
subroutine load14(nphi,ixi,ixj,ixk,ixl,ixm,x,n14,ni14,ip)
   implicit double precision (a-h,o-z)
   dimension n14(*),ni14(*),x(*)
   dimension ixi(*),ixj(*),ixk(*),ixl(*),ixm(*)
   ip = 0
   do 511 i = 1,nphi
      i1 = ixi(i)/3+1
      i3 = abs(ixk(i))/3+1
      i4 = abs(ixl(i))/3+1
      if(ixk(i) >= 0 .and. ixl(i) >= 0 )then
         n14(i1) = n14(i1) + 1
         if(n14(i1) > 15) then
            print *,'too many dihedrals '
            call mexit(6,1)
         end if
         call istuff(n14(i1),i1,ni14,i4)
         ip = ip+1
      end if
   511 continue
   return
end subroutine load14 
