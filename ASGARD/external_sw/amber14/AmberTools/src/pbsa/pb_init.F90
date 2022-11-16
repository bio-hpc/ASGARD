#include "copyright.h"
#include "../include/dprec.fh"
#include "pb_def.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PBMD cleanup routine
subroutine pb_free

   use poisson_boltzmann
   use solvent_accessibility

   implicit none
#  include "../include/md.h"
#ifdef SANDER
#  include "../sander/box.h"
#else
#  include "box.h"
#endif

   integer alloc_err(32)

   alloc_err(1:32) = 0

   if ( allocated(icrd   ) ) deallocate(   icrd, stat = alloc_err(1 ) )
   if ( allocated(grdcrg ) ) deallocate( grdcrg, stat = alloc_err(2 ) )
   if ( allocated(qgrdcrg) ) deallocate(qgrdcrg, stat = alloc_err(3 ) )

   if ( allocated(acrg) ) deallocate( acrg, stat = alloc_err(4 ) )
   if ( allocated(gcrg) ) deallocate( gcrg, stat = alloc_err(5 ) )

   if ( allocated(nshrt) ) deallocate(nshrt, stat = alloc_err(6 ) )
   if ( allocated(nex  ) ) deallocate(  nex, STAT = alloc_err(7 ) )
   if ( allocated(iex  ) ) deallocate(  iex, STAT = alloc_err(8 ) )
    
   if ( allocated(acrd) ) deallocate( acrd, stat = alloc_err(9 ) )
   if ( allocated(gcrd) ) deallocate( gcrd, stat = alloc_err(10) )

   if ( allocated(mdsig  ) ) deallocate(  mdsig, stat = alloc_err(11) )
   if ( allocated(rmin   ) ) deallocate(   rmin, stat = alloc_err(12) )
   if ( allocated(radi   ) ) deallocate(   radi, stat = alloc_err(13) )
   if ( allocated(radip  ) ) deallocate(  radip, stat = alloc_err(14) )
   if ( allocated(radip2 ) ) deallocate( radip2, stat = alloc_err(15) )
   if ( allocated(radip3 ) ) deallocate( radip3, stat = alloc_err(16) )
   if ( allocated(nzratm ) ) deallocate( nzratm, stat = alloc_err(17) )
   if ( allocated(nmax   ) ) deallocate(   nmax, stat = alloc_err(18) )
   if ( allocated(nexp   ) ) deallocate(   nexp, stat = alloc_err(19) )
   if ( allocated(sumnmax) ) deallocate(sumnmax, stat = alloc_err(20) )
   if ( allocated(sumnexp) ) deallocate(sumnexp, stat = alloc_err(21) )
   if ( allocated(avnmax ) ) deallocate( avnmax, stat = alloc_err(22) )
   if ( allocated(avnexp ) ) deallocate( avnexp, stat = alloc_err(23) )
   if ( allocated(scrd   ) ) deallocate(   scrd, stat = alloc_err(24) )
 
   if ( allocated(iar1pb ) ) deallocate( iar1pb, stat = alloc_err(25) )
   if ( allocated(iprshrt) ) deallocate(iprshrt, stat = alloc_err(26) )
   if ( allocated(cn1pb  ) ) deallocate(  cn1pb, stat = alloc_err(28) )
   if ( allocated(cn2pb  ) ) deallocate(  cn2pb, stat = alloc_err(29) )
   if ( allocated(cn3pb  ) ) deallocate(  cn3pb, stat = alloc_err(30) )

   if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3)+alloc_err( 4)+alloc_err( 5)+&
        alloc_err( 6)+alloc_err( 7)+alloc_err( 8)+alloc_err( 9)+alloc_err(10)+&
        alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+alloc_err(15)+&
        alloc_err(16)+alloc_err(17)+alloc_err(18)+alloc_err(19)+alloc_err(20)+&
        alloc_err(21)+alloc_err(22)+alloc_err(23)+alloc_err(24)+alloc_err(25)+&
        alloc_err(26)+alloc_err(27)+alloc_err(28)+alloc_err(29)+alloc_err(30)+&
        alloc_err(31)+alloc_err(32) /= 0 ) then
      write(6,'(a,i6)') 'PB Bomb in pb_free(): Deallocation aborted #1', alloc_err(1:32)
      call mexit(6, 1)
   end if

   alloc_err(1:32) = 0

   if ( allocated(phi    ) ) deallocate(    phi, stat = alloc_err(1 ) )
   if ( allocated(chgrd  ) ) deallocate(  chgrd, stat = alloc_err(2 ) )
   if ( allocated(epsx   ) ) deallocate(   epsx, stat = alloc_err(3 ) )
   if ( allocated(epsy   ) ) deallocate(   epsy, stat = alloc_err(4 ) )
   if ( allocated(epsz   ) ) deallocate(   epsz, stat = alloc_err(5 ) )
   if ( allocated(saltgrd) ) deallocate(saltgrd, stat = alloc_err(6 ) )

   if ( allocated(bv) ) deallocate( bv, stat = alloc_err(7 ) )

   if ( allocated(insas ) ) deallocate( insas, stat = alloc_err(8 ) )
   if ( allocated(atmsas) ) deallocate(atmsas, stat = alloc_err(9 ) )
   if ( allocated(lvlset) ) deallocate(lvlset, stat = alloc_err(10) )
   if ( allocated(zv    ) ) deallocate(    zv, stat = alloc_err(11) )

   if ( allocated(cphi) ) deallocate( cphi, stat = alloc_err(12) )

   if ( allocated(fedgex ) ) deallocate( fedgex, stat = alloc_err(13) )
   if ( allocated(fedgey ) ) deallocate( fedgey, stat = alloc_err(14) )
   if ( allocated(fedgez ) ) deallocate( fedgez, stat = alloc_err(15) )
   if ( allocated(iepsav ) ) deallocate( iepsav, stat = alloc_err(16) )
   if ( allocated(iepsavx) ) deallocate(iepsavx, stat = alloc_err(17) )
   if ( allocated(iepsavy) ) deallocate(iepsavy, stat = alloc_err(18) )
   if ( allocated(iepsavz) ) deallocate(iepsavz, stat = alloc_err(19) )

   if ( allocated(xs) ) deallocate( xs, stat = alloc_err(20) )

   if ( allocated(outflag    ) ) deallocate(    outflag, stat = alloc_err(21) )
   if ( allocated(outflagorig) ) deallocate(outflagorig, stat = alloc_err(22) )
   if ( allocated(mapout     ) ) deallocate(     mapout, stat = alloc_err(23) )

   if ( allocated(liveflag) ) deallocate( liveflag, stat = alloc_err(24) )
   if ( allocated(realflag) ) deallocate( realflag, stat = alloc_err(25) )

   if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3)+alloc_err( 4)+alloc_err( 5)+&
        alloc_err( 6)+alloc_err( 7)+alloc_err( 8)+alloc_err( 9)+alloc_err(10)+&
        alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+alloc_err(15)+&
        alloc_err(16)+alloc_err(17)+alloc_err(18)+alloc_err(19)+alloc_err(20)+&
        alloc_err(21)+alloc_err(22)+alloc_err(23)+alloc_err(24)+alloc_err(25)+&
        alloc_err(26)+alloc_err(27)+alloc_err(28)+alloc_err(29)+alloc_err(30)+&
        alloc_err(31)+alloc_err(32) /= 0 ) then
      write(6,'(a,i6)') 'PB Bomb in pb_free(): Deallocation aborted #2', alloc_err(1:32)
      call mexit(6, 1)
   end if

end subroutine pb_free
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PBMD initialization routine
subroutine pb_init(ifcap, natom, nres, ntypes, nbonh, nbona, ipres, iac, ico, &
                   numex, natex, ibh, jbh, iba, jba, ibel, lbres, igraph, &
                   isymbl, cg, rin)
     
   ! Module variables
     
   use parms, only : cn1, cn2
   use poisson_boltzmann
   use solvent_accessibility
   use dispersion_cavity

   implicit none
     
   ! Common variables

#  include "../include/md.h"
#  include "pb_md.h"
     
   ! Passed variables
      
   integer ifcap,natom, nres, ntypes, nbonh, nbona
   integer ipres(*), iac(*), ico(*), numex(*), natex(*), ibh(*), jbh(*), iba(*), jba(*), ibel(*)
   character (len=4) :: lbres(*), igraph(*), isymbl(*)
   _REAL_ cg(natom), rin(natom)
     
   ! Local variables
     
   integer ires, iatm, jatm, maxmax, ic, i, j, jp, idum
   integer alloc_err(64)
   _REAL_ maxnba_l
   _REAL_ rinchk, crgchk
   _REAL_ ucrgh(natom), ucrga(natom)
   character (len=4) :: residue, resid(natom)
   integer focusresidue(nres)
    
   ! begin of code
     
   alloc_err(1:64) = 0

   ! allocate topology informations
     
   allocate(   icrd( 3,  natom), stat = alloc_err(1 ) )
   allocate( grdcrg( 3,8*natom), stat = alloc_err(2 ) )
   allocate(qgrdcrg(   8*natom), stat = alloc_err(3 ) )
   allocate(   acrg(     natom), stat = alloc_err(4 ) )
   allocate(   gcrg( 8,  natom), stat = alloc_err(5 ) )
   allocate(  nshrt( 0:  natom), stat = alloc_err(6 ) )
   allocate(    nex(     natom), STAT = alloc_err(7 ) )
   allocate(    iex(64,  natom), STAT = alloc_err(8 ) )
    
   if ( ifcap == 0 .or. (ifcap >= 3 .and. ifcap <= 5) ) then

   allocate(   acrd( 3,  natom), stat = alloc_err(9 ) )
   allocate(   gcrd( 3,  natom), stat = alloc_err(10) )
    
   ! allocate sas informations
    
   allocate(  mdsig(  natom  ), stat = alloc_err(11) )
   allocate(   rmin(  natom  ), stat = alloc_err(11) )
   allocate(   radi(  natom  ), stat = alloc_err(12) )
   allocate(  radip(  natom  ), stat = alloc_err(13) )
   allocate( radip2(  natom  ), stat = alloc_err(14) )
   allocate( radip3(  natom  ), stat = alloc_err(15) )
   allocate( nzratm(  natom  ), stat = alloc_err(16) )
   allocate(   nmax(  natom  ), stat = alloc_err(17) )
   allocate(   nexp(  natom  ), stat = alloc_err(18) )
   allocate(sumnmax(  natom  ), stat = alloc_err(19) )
   allocate(sumnexp(  natom  ), stat = alloc_err(20) )
   allocate( avnmax(  natom  ), stat = alloc_err(21) )
   allocate( avnexp(  natom  ), stat = alloc_err(22) )
   allocate(   scrd(3,maxsph ), stat = alloc_err(23) )
 
   else
 
   allocate(   acrd( 3, 0:natom), stat = alloc_err(9 ) )
   allocate(   gcrd( 3, 0:natom), stat = alloc_err(10) )

   allocate(   radi(    0: 0   ), stat = alloc_err(12) )
   allocate( radip3(    1: 1   ), stat = alloc_err(15) )
   allocate( nzratm(    1: 1   ), stat = alloc_err(16) )

   end if
    
   ! allocate pb nblists
   
   maxnba_l = dble(natom) * ( sqrt( max(cutnb,cutsa,cutfd) ) )**3 / 3.0d0
   if ( natom >= 65536 ) then
      write(6,'(a)') "PB Warning: natom**2 exceeds integer limit (2147483647)."
      maxmax = 2147483647
   else
      maxmax = ceiling(dble(natom)/2*dble(natom))
   end if
   if ( maxnba_l > maxmax ) then
      maxnba = maxmax
   else
      maxnba = int(maxnba_l)
   end if

   allocate( iar1pb (4,0:natom), stat = alloc_err(24) )
   allocate( iprshrt(  maxnba ), stat = alloc_err(25) )
   allocate( cn1pb  (  maxnba ), stat = alloc_err(27) )
   allocate( cn2pb  (  maxnba ), stat = alloc_err(28) )
   allocate( cn3pb  (  maxnba ), stat = alloc_err(29) )

   ! allocate ibelly and outflag

   if ( ifcap == 3 .or. ifcap == 4 ) then
     allocate( ibelly(natom), stat = alloc_err(30) )
   end if
   allocate( outflag(  natom  ), stat = alloc_err(31) )
   allocate( outflagorig(  natom  ), stat = alloc_err(32) )
   allocate( mapout(  natom  ), stat = alloc_err(33) )

   if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3)+alloc_err( 4)+alloc_err( 5)+&
        alloc_err( 6)+alloc_err( 7)+alloc_err( 8)+alloc_err( 9)+alloc_err(10)+&
        alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+alloc_err(15)+&
        alloc_err(16)+alloc_err(17)+alloc_err(18)+alloc_err(19)+alloc_err(20)+&
        alloc_err(21)+alloc_err(22)+alloc_err(23)+alloc_err(24)+alloc_err(25)+&
        alloc_err(26)+alloc_err(27)+alloc_err(28)+alloc_err(29)+alloc_err(30)+&
        alloc_err(31)+alloc_err(32)+alloc_err(33) /= 0 ) then
      write(6,'(a,i6)') 'PB Bomb in pb_init(): Allocation aborted', alloc_err(1:31)
      call mexit(6, 1)
   end if 

   if ( pbverbose ) then
      write(6,'(a)')
      write(6,'(a)') '======== Implicit Solvent Initialization ========'
      write(6,'(a)')
      if ( pqropt == 0 ) write(6,'(5x,a,2i9)') 'Max Nonbonded Pairs:', maxnba, maxmax
      write(6,'(a)')
   end if

   ! the following setups are mostly for standard prmtop/inpcrd files

   if ( pqropt == 0 ) then

   ! getting some topology info into pb setup ...

   do ires = 1, nres 
      write(residue,'(a4)') lbres(ires)
      do iatm = ipres(ires), ipres(ires+1) - 1
         resid(iatm) = residue
      enddo
   enddo

   ! atomic and total charge in electron for printing only

   totcrgp = ZERO
   totcrgn = ZERO
   do iatm = 1, natom
       
      acrg(iatm) = cg(iatm)*INV_AMBER_ELECTROSTATIC
      if ( acrg(iatm) > ZERO) then
         totcrgp = totcrgp + acrg(iatm)
      else
         totcrgn = totcrgn + acrg(iatm)
      end if
       
      ! get info about belly atoms for GBSP

      if(ifcap == 3 .or. ifcap == 4) then
        ibelly(iatm) = ibel(iatm)
      end if

   end do
   totcrg = totcrgp + totcrgn

   ! for pure implicit solvent, set up radii arrays

   if ( ifcap == 0 .or. (ifcap >= 3 .and. ifcap <= 5) ) then

      ! set up group charges for cavity radii analysis

      ucrgh(1:natom) = acrg(1:natom)
      do idum = 1, nbonh
         iatm = ibh(idum)/3 + 1
         jatm = jbh(idum)/3 + 1
         if (isymbl(iatm)(1:1) == 'H' ) then
            ucrgh(jatm) = ucrgh(jatm) + acrg(iatm)
         else
            ucrgh(iatm) = ucrgh(iatm) + acrg(jatm)
         endif
      enddo
                                                                                                          
      ucrga(1:natom) = ucrgh(1:natom)
      do idum = 1, nbona
         iatm = iba(idum)/3 + 1
         jatm = jba(idum)/3 + 1
         ucrga(iatm) = ucrga(iatm) + ucrgh(jatm)
         ucrga(jatm) = ucrga(jatm) + ucrgh(iatm)
      enddo
                                                                                                       
      do idum = 1, natom
         if (isymbl(idum)(1:1) == 'H' ) then
            ucrga(idum) = 0.0d0
            ucrgh(idum) = 0.0d0
         endif
      enddo

      ! van der Waals sigma radii for nonpolar solvation
       
      do iatm = 1, natom
         ic = ico(ntypes*(iac(iatm)-1) + iac(iatm))
         if (cn2(ic) /= ZERO) then
            mdsig(iatm) = (cn1(ic)/cn2(ic))**(SIXTH)/2 ! this is sigma
            rmin(iatm) = mdsig(iatm)*(2.0d0**(SIXTH)) ! this is Rmin 
         else
            mdsig(iatm) = ZERO
            rmin(iatm) = ZERO
         endif
      end do

      ! cavity radii for polar solvation if not passing in from driver

      if ( radiopt == 0 ) then
         rinchk = ZERO
         do iatm = 1, natom
            rinchk = rinchk+rin(iatm)
         end do
         if (rinchk == ZERO) then
            write(6,'(a)') ' PB Bomb in pb_init(): Requested radi to be read in, but found none'
            !call mexit(6,1)
            if ( ipb > 0 ) then 
               write(6,'(a)') ' PB Bomb in pb_init(): No PB calculation done'
               write(6,'(a)') ' '
               ipb = 0
               igb = 100
            end if
            if ( inp == 1 ) then
               write(6,'(a)') ' PB Bomb in pb_init(): No NP calculation done'
               write(6,'(a)') ' '
               inp = 0
            end if
         end if
         radi = rin ! for pb
         mdsig = rin ! for np
      else if ( radiopt == 1 ) then
         call pb_aaradi( natom, nbonh, ibh, jbh, radi, acrg, ucrgh, ucrga, resid, igraph, isymbl, rin )
         radi = radi*radiscale
      else if ( radiopt == 2 ) then
         call phi_aaradi( natom, isymbl, radi )
      else
         write(6,'(a,i6)') 'PB Bomb in pb_init(): Unknown radi assigment option', radiopt
         !call mexit(6,1)
         if ( ipb > 0 ) then
            write(6,'(a)') 'PB Bomb in pb_init(): No PB calculation done'
            write(6,'(a)') ' '
            ipb = 0
            igb = 100
         end if
         if ( inp == 1 ) then 
            write(6,'(a)') 'PB Bomb in pb_init(): No NP calculation done'
            write(6,'(a)') ' '
            inp = 0
         end if
      end if

      ! check if atomic charges are correctly input from prmtop file

      crgchk = ZERO
      do iatm = 1, natom
         crgchk = crgchk + abs(acrg(iatm))
      end do
      if ( crgchk < 1.0D-6 ) then
         write(6,'(a)') ' PB Bomb in pb_init(): Requested charges to be read in, but found none'
         write(6,'(a)') ' PB Bomb in pb_init(): No PB calculation done'
         write(6,'(a)') ' '
         ipb = 0
         igb = 100
      end if

      ! let's do a summary of key parameters for pb before proceeding

      if ( pbverbose ) then
         write(6,'(a,i6)') ' no. of atoms processed in PB initialization:', natom
         write(6, '(a5,3a6,5a10)') 'NUM','RESI','NAME','TYPE',&
         'CHARGE', 'ATM CRG/H', 'GRP CRG', 'PB RADI', 'NP RADI'
         do iatm = 1, natom
            if ( use_rmin == 0 ) write(6, '(i5,3a6,5f10.6)') iatm,resid(iatm),igraph(iatm),isymbl(iatm),&
            real(acrg(iatm)),real(ucrgh(iatm)), real(ucrga(iatm)),&
            real(radi(iatm)), real(mdsig(iatm))
            if ( use_rmin == 1 ) write(6, '(i5,3a6,5f10.6)') iatm,resid(iatm),igraph(iatm),isymbl(iatm),&
            real(acrg(iatm)),real(ucrgh(iatm)), real(ucrga(iatm)),&
            real(radi(iatm)), real(rmin(iatm))
         end do
         write(6,*)
         write(6,'(a,3f14.4)') ' total system charges (+/-) for PB', totcrg, totcrgp, totcrgn
         write(6,'(a,f14.4,a,f14.4)') ' cavity_surften =', cavity_surften, ' cavity_offset =', cavity_offset
         write(6,*)
      end if

      ! initialization for sas surface

      if ( srsas ) then

         call sa_sphere(maxsph, scrd)
         if ( pbverbose ) write(6,'(2x,a,i6)') ' SAS Surface: surface dots generated: ', maxsph

         ! nmax and nexp accumulators

         sumnmax(1:natom) = ZERO
         sumnexp(1:natom) = ZERO

      ! initialization for for vdw surface

      else
         if ( pbverbose ) write(6,'(a)') ' VDW Surface: setting up working radii'
         radip3(1:natom) = radi(1:natom)
      end if

   else
      acrd(1:3,0:0) = ZERO; gcrd(1:3,0:0) = ZERO
   end if ! if ( ifcap == 0 .or. (ifcap >= 3 .and. ifcap <= 5) ) then

   ! assinging atom-based pointers to exclusion list

   nshrt(0) = 0
   do i = 1, natom
      nshrt(i) = nshrt(i-1) + numex(i)
   end do
   nex = 0
   do i = 1, natom-1
      do jp = nshrt(i-1) + 1, nshrt(i)
         j = natex(jp)
         if (j == 0)  cycle
         nex(i) = nex(i) + 1
         iex(nex(i),i) = j
         nex(j) = nex(j) + 1
         iex(nex(j),j) = i
      end do
   end do

   ! only the following setups are needed/can be done for the pqr file

   else

      ! apparently there is exclusion list

      nshrt = 0

      ! initialization for sas surface

      if ( srsas ) then

         call sa_sphere(maxsph, scrd)
         if ( pbverbose ) write(6,'(2x,a,i6)') ' SAS Surface: surface dots generated: ', maxsph

         ! nmax and nexp accumulators

         sumnmax(1:natom) = ZERO
         sumnexp(1:natom) = ZERO

      end if

   end if ! if ( pqropt == 0 ) then

   ! set up green 3-d array for pb_fdcoulomb

   call pb_green

   ! set up variables if using ligand or multiblock

   if ( ligand .or. multiblock ) then

      allocate( liveflag(natom), stat = alloc_err(1) )
      allocate( realflag(natom), stat = alloc_err(2) )
      if ( SUM(alloc_err(1:2)) /= 0 ) then
         write(6,'(a,i6)') ' PB Bomb in pb_init(): Allocation aborted', alloc_err(1:2)
         call mexit(6,1)
      end if
     
      ! Mengjuei -
      ! Basically liveflag HERE is a temporary switch storing the list of
      ! atoms labelled ligandmask. However if the geormetries of ligand
      ! box were specified, we don't need liveflag that early.
      ! Later liveflag will become the list of atom inside the block of
      ! rectangular bounding box.
     
      if ( multiblock ) then
         continue!noop
      else if ( xmax - xmin > ZERO .and. &
                ymax - ymin > ZERO .and. &
                zmax - zmin > ZERO       ) then
         ! If the geometries of bounding box are correctly specified,
         ! skip the ligandmask
         continue
      else if ( ligand .and. len_trim(ligandmask) > 0 ) then
         focusresidue = 0
             liveflag = 0
         call myresmask ( ligandmask, LEN_TRIM(ligandmask), focusresidue, nres )
         j = 0
         do i = 1, nres
            if ( focusresidue(i) == 0 ) cycle
            j = j + 1
            liveflag(ipres(i):ipres(i+1)-1)=1
         enddo
         write(6,'(a,a,a,i5)') ' Focusing Mask ', &
              ligandmask(1:LEN_TRIM(ligandmask)), ' matches residue', j
         xmax = ZERO; ymax = ZERO; zmax = ZERO
         xmin = ZERO; ymin = ZERO; zmin = ZERO
      else
         write(6,'(a)') 'wrong ligand box'
         call mexit(6,1)
      endif

   endif

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Setup table of precomputed FD Green's function.
subroutine pb_green

   implicit none

   ! Common variables

   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green
   
   _REAL_ rdsqr(0:20, 0:20, 0:20)
   common /blk_rdsqr/ rdsqr
   
   ! Local variables
   
   integer i, j, k
   _REAL_ green_data

   do i = 0, 20
      do j = 0, i
         do k = 0, j
            green_data = green(i, j, k)
            green(i, k, j) = green_data
            green(j, i, k) = green_data
            green(k, i, j) = green_data
            green(j, k, i) = green_data
            green(k, j, i) = green_data
         end do
      end do
   end do
   do i = 0, 20
      do j = 0, 20
         do k = 0, 20
            if ( i == 0 .and. j == 0 .and. k == 0 ) then
               green_data = 0.0d0
            else
               green_data = 1.0d0/sqrt( dble(i**2 + j**2 + k**2) )
            end if
            rdsqr(i, j, k) = green_data
         end do
      end do
   end do

end subroutine pb_green

end subroutine pb_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Data for precomputed FD Green's function
block data blkgreen

   implicit none
   
   ! Common variables
   
   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green
  
   data green( 0, 0, 0) /  3.1759115343998565      / 
   data green( 1, 0, 0) /  1.0815164320066608      /
   data green( 1, 1, 0) / 0.69355600923256988      /
   data green( 1, 1, 1) / 0.54762175090034226      /
   data green( 2, 0, 0) / 0.53896302070982693      /
   data green( 2, 1, 0) / 0.45152984479070624      /
   data green( 2, 1, 1) / 0.40168749256811487      /
   data green( 2, 2, 0) / 0.35255169528423141      /
   data green( 2, 2, 1) / 0.32872036340200106      /
   data green( 2, 2, 2) / 0.28464545919844825      /
   data green( 3, 0, 0) / 0.34614231308947585      /
   data green( 3, 1, 0) / 0.32073335838137879      /
   data green( 3, 1, 1) / 0.30200278812293213      /
   data green( 3, 2, 0) / 0.27740487765998695      /
   data green( 3, 2, 1) / 0.26587502039654837      /
   data green( 3, 2, 2) / 0.24057055499489527      /
   data green( 3, 3, 0) / 0.23517700551562060      /
   data green( 3, 3, 1) / 0.22833414279274872      /
   data green( 3, 3, 2) / 0.21176997182120025      /
   data green( 3, 3, 3) / 0.19125147413327950      /
   data green( 4, 0, 0) / 0.25495742430151219      /
   data green( 4, 1, 0) / 0.24531753850223970      /
   data green( 4, 1, 1) / 0.23711247861362333      /
   data green( 4, 2, 0) / 0.22421716598559402      /
   data green( 4, 2, 1) / 0.21821739540672588      /
   data green( 4, 2, 2) / 0.20348788633542600      /
   data green( 4, 3, 0) / 0.19979199609412604      /
   data green( 4, 3, 1) / 0.19565391931328732      /
   data green( 4, 3, 2) / 0.18494655200569138      /
   data green( 4, 3, 3) / 0.17073297644535870      /
   data green( 4, 4, 0) / 0.17649822095805937      /
   data green( 4, 4, 1) / 0.17368106314142176      /
   data green( 4, 4, 2) / 0.16610647203353304      /
   data green( 4, 4, 3) / 0.15557916272053965      /
   data green( 4, 4, 4) / 0.14383287016300686      /
   data green( 5, 0, 0) / 0.20233207871063805      /
   data green( 5, 1, 0) / 0.19777232511770629      /
   data green( 5, 1, 1) / 0.19360221574087669      /
   data green( 5, 2, 0) / 0.18635379284375955      /
   data green( 5, 2, 1) / 0.18295790179587609      /
   data green( 5, 2, 2) / 0.17402886819282612      /
   data green( 5, 3, 0) / 0.17155174547890764      /
   data green( 5, 3, 1) / 0.16895236643900982      /
   data green( 5, 3, 2) / 0.16192808608534290      /
   data green( 5, 3, 3) / 0.15209495508641044      /
   data green( 5, 4, 0) / 0.15602160363863021      /
   data green( 5, 4, 1) / 0.15408692361518175      /
   data green( 5, 4, 2) / 0.14874275116392699      /
   data green( 5, 4, 3) / 0.14103484061799032      /
   data green( 5, 4, 4) / 0.13208657760547401      /
   data green( 5, 5, 0) / 0.14126545291334486      /
   data green( 5, 5, 1) / 0.13983695654019662      /
   data green( 5, 5, 2) / 0.13582392009095975      /
   data green( 5, 5, 3) / 0.12988672367513968      /
   data green( 5, 5, 4) / 0.12279310610564238      /
   data green( 5, 5, 5) / 0.11521211763561941      /
   data green( 6, 0, 0) / 0.16794574749149105      /
   data green( 6, 1, 0) / 0.16542610916784664      /
   data green( 6, 1, 1) / 0.16304036200447203      /
   data green( 6, 2, 0) / 0.15866571688859704      /
   data green( 6, 2, 1) / 0.15659277215205825      /
   data green( 6, 2, 2) / 0.15091334705909260      /
   data green( 6, 3, 0) / 0.14923834741891018      /
   data green( 6, 3, 1) / 0.14753562234546314      /
   data green( 6, 3, 2) / 0.14280302362419248      /
   data green( 6, 3, 3) / 0.13591090066643738      /
   data green( 6, 4, 0) / 0.13864035525110568      /
   data green( 6, 4, 1) / 0.13728680076790511      /
   data green( 6, 4, 2) / 0.13347626454055395      /
   data green( 6, 4, 3) / 0.12781887345645104      /
   data green( 6, 4, 4) / 0.12103070202257171      /
   data green( 6, 5, 0) / 0.12793779856120766      /
   data green( 6, 5, 1) / 0.12687925950325574      /
   data green( 6, 5, 2) / 0.12386716900128411      /
   data green( 6, 5, 3) / 0.11931681730912760      /
   data green( 6, 5, 4) / 0.11374332005607347      /
   data green( 6, 5, 5) / 0.10763112916559645      /
   data green( 6, 6, 0) / 0.11775701708375737      /
   data green( 6, 6, 1) / 0.11693365296828663      /
   data green( 6, 6, 2) / 0.11457015905523310      /
   data green( 6, 6, 3) / 0.11094579188109956      /
   data green( 6, 6, 4) / 0.10642359184315756      /
   data green( 6, 6, 5) / 0.10136496364921296      /
   data green( 6, 6, 6) / 9.60759687686214559D-002 /
   data green( 7, 0, 0) / 0.14363796956692135      /
   data green( 7, 1, 0) / 0.14209214150034138      /
   data green( 7, 1, 1) / 0.14060219364614540      /
   data green( 7, 2, 0) / 0.13779050759694922      /
   data green( 7, 2, 1) / 0.13644368281884856      /
   data green( 7, 2, 2) / 0.13265962260922792      /
   data green( 7, 3, 0) / 0.13150102220392421      /
   data green( 7, 3, 1) / 0.13034042367070278      /
   data green( 7, 3, 2) / 0.12705392104826485      /
   data green( 7, 3, 3) / 0.12212665475092678      /
   data green( 7, 4, 0) / 0.12407078035207596      /
   data green( 7, 4, 1) / 0.12310237935187017      /
   data green( 7, 4, 2) / 0.12033896922956387      /
   data green( 7, 4, 3) / 0.11614371558202512      /
   data green( 7, 4, 4) / 0.11097324750490715      /
   data green( 7, 5, 0) / 0.11620544711252656      /
   data green( 7, 5, 1) / 0.11541317918065419      /
   data green( 7, 5, 2) / 0.11313659350857447      /
   data green( 7, 5, 3) / 0.10963902578471769      /
   data green( 7, 5, 4) / 0.10526457389034512      /
   data green( 7, 5, 5) / 0.10035808994738646      /
   data green( 7, 6, 0) / 0.10839959972177780      /
   data green( 7, 6, 1) / 0.10775811133210896      /
   data green( 7, 6, 2) / 0.10590358573972203      /
   data green( 7, 6, 3) / 0.10302368288497579      /
   data green( 7, 6, 4) / 9.93720777082429535D-002 /
   data green( 7, 6, 5) / 9.52139814761528286D-002 /
   data green( 7, 6, 6) / 9.07869738880299482D-002 /
   data green( 7, 7, 0) / 0.10095457662354870      /
   data green( 7, 7, 1) / 0.10043709528831846      /
   data green( 7, 7, 2) / 9.89333081671740916D-002 /
   data green( 7, 7, 3) / 9.65760097217587543D-002 /
   data green( 7, 7, 4) / 9.35500217771100540D-002 /
   data green( 7, 7, 5) / 9.00559653631305690D-002 /
   data green( 7, 7, 6) / 8.62815540718346874D-002 /
   data green( 7, 7, 7) / 8.23848148089189602D-002 /
   data green( 8, 0, 0) / 0.12551350390867147      /
   data green( 8, 1, 0) / 0.12449387537804016      /
   data green( 8, 1, 1) / 0.12350115123402043      /
   data green( 8, 2, 0) / 0.12159679935113553      /
   data green( 8, 2, 1) / 0.12067657723800752      /
   data green( 8, 2, 2) / 0.11804918086204812      /
   data green( 8, 3, 0) / 0.11722565051420447      /
   data green( 8, 3, 1) / 0.11640591425584559      /
   data green( 8, 3, 2) / 0.11405483240497505      /
   data green( 8, 3, 3) / 0.11045375457854335      /
   data green( 8, 4, 0) / 0.11187309884115869      /
   data green( 8, 4, 1) / 0.11116412291031905      /
   data green( 8, 4, 2) / 0.10912094134609471      /
   data green( 8, 4, 3) / 0.10596552276558401      /
   data green( 8, 4, 4) / 0.10199220406213064      /
   data green( 8, 5, 0) / 0.10599814567878951      /
   data green( 8, 5, 1) / 0.10539728427558927      /
   data green( 8, 5, 2) / 0.10365763211550470      /
   data green( 8, 5, 3) / 0.10094877153325783      /
   data green( 8, 5, 4) / 9.75016823407428290D-002 /
   data green( 8, 5, 5) / 9.35602997857264695D-002 /
   data green( 8, 6, 0) / 9.99643348466161402D-002 /
   data green( 8, 6, 1) / 9.94615550938944826D-002 /
   data green( 8, 6, 2) / 9.79996594902655810D-002 /
   data green( 8, 6, 3) / 9.57056064743137674D-002 /
   data green( 8, 6, 4) / 9.27566143777162483D-002 /
   data green( 8, 6, 5) / 8.93458183009139539D-002 /
   data green( 8, 6, 6) / 8.56548034635831873D-002 /
   data green( 8, 7, 0) / 9.40270348605498768D-002 /
   data green( 8, 7, 1) / 9.36092321374849545D-002 /
   data green( 8, 7, 2) / 9.23897862567616135D-002 /
   data green( 8, 7, 3) / 9.04626813081583908D-002 /
   data green( 8, 7, 4) / 8.79620000806425051D-002 /
   data green( 8, 7, 5) / 8.50381266887664938D-002 /
   data green( 8, 7, 6) / 8.18372982414493216D-002 /
   data green( 8, 7, 7) / 7.84880755460032331D-002 /
   data green( 8, 8, 0) / 8.83470984223495293D-002 /
   data green( 8, 8, 1) / 8.80008148673751822D-002 /
   data green( 8, 8, 2) / 8.69867315876210745D-002 /
   data green( 8, 8, 3) / 8.53741053967011637D-002 /
   data green( 8, 8, 4) / 8.32636437937793322D-002 /
   data green( 8, 8, 5) / 8.07712357260426672D-002 /
   data green( 8, 8, 6) / 7.80130771303589182D-002 /
   data green( 8, 8, 7) / 7.50949080072080694D-002 /
   data green( 8, 8, 8) / 7.21059876168646668D-002 /
   data green( 9, 0, 0) / 0.11146755237294673      /
   data green( 9, 1, 0) / 0.11075850504005154      /
   data green( 9, 1, 1) / 0.11006380852588185      /
   data green( 9, 2, 0) / 0.10871760814160422      /
   data green( 9, 2, 1) / 0.10806273490614675      /
   data green( 9, 2, 2) / 0.10617264327709565      /
   data green( 9, 3, 0) / 0.10557115417731718      /
   data green( 9, 3, 1) / 0.10497387879686464      /
   data green( 9, 3, 2) / 0.10324528233905353      /
   data green( 9, 3, 3) / 0.10055516237921519      /
   data green( 9, 4, 0) / 0.10161577068124410      /
   data green( 9, 4, 1) / 0.10108511939135567      /
   data green( 9, 4, 2) / 9.95445686506214039D-002 /
   data green( 9, 4, 3) / 9.71337494914521898D-002 /
   data green( 9, 4, 4) / 9.40455666552229980D-002 /
   data green( 9, 5, 0) / 9.71514247212569848D-002 /
   data green( 9, 5, 1) / 9.66890706743735034D-002 /
   data green( 9, 5, 2) / 9.53425425392463999D-002 /
   data green( 9, 5, 3) / 9.32231597186838623D-002 /
   data green( 9, 5, 4) / 9.04876303952805855D-002 /
   data green( 9, 5, 5) / 8.73087074836587651D-002 /
   data green( 9, 6, 0) / 9.24381186307905583D-002 /
   data green( 9, 6, 1) / 9.20407084813021359D-002 /
   data green( 9, 6, 2) / 9.08797912613967945D-002 /
   data green( 9, 6, 3) / 8.90422292515088282D-002 /
   data green( 9, 6, 4) / 8.66525013614414252D-002 /
   data green( 9, 6, 5) / 8.38510840135384267D-002 /
   data green( 9, 6, 6) / 8.07756138087426800D-002 /
   data green( 9, 7, 0) / 8.76777349958148544D-002 /
   data green( 9, 7, 1) / 8.73391064580100601D-002 /
   data green( 9, 7, 2) / 8.63471048498655053D-002 /
   data green( 9, 7, 3) / 8.47685799187723604D-002 /
   data green( 9, 7, 4) / 8.27009125383244142D-002 /
   data green( 9, 7, 5) / 8.02564424204199184D-002 /
   data green( 9, 7, 6) / 7.75481525481492295D-002 /
   data green( 9, 7, 7) / 7.46792259697855731D-002 /
   data green( 9, 8, 0) / 8.30134455391234732D-002 /
   data green( 9, 8, 1) / 8.27262974596552347D-002 /
   data green( 9, 8, 2) / 8.18829483740634301D-002 /
   data green( 9, 8, 3) / 8.05344471912449594D-002 /
   data green( 9, 8, 4) / 7.87562607393235553D-002 /
   data green( 9, 8, 5) / 7.66372200272922993D-002 /
   data green( 9, 8, 6) / 7.42688612830020994D-002 /
   data green( 9, 8, 7) / 7.17371161020092102D-002 /
   data green( 9, 8, 8) / 6.91170672265212366D-002 /
   data green( 9, 9, 0) / 7.85381639271570953D-002 /
   data green( 9, 9, 1) / 7.82951291375275499D-002 /
   data green( 9, 9, 2) / 7.75797037762221975D-002 /
   data green( 9, 9, 3) / 7.64307462261306220D-002 /
   data green( 9, 9, 4) / 7.49064582104172622D-002 /
   data green( 9, 9, 5) / 7.30766032903040347D-002 /
   data green( 9, 9, 6) / 7.10146806477706255D-002 /
   data green( 9, 9, 7) / 6.87914613338068415D-002 /
   data green( 9, 9, 8) / 6.64705444915650590D-002 /
   data green( 9, 9, 9) / 6.41059460647893098D-002 /
   data green(10, 0, 0) / 0.10025779016880276      /
   data green(10, 1, 0) / 9.97443772959542146D-002 /
   data green(10, 1, 1) / 9.92392200288740867D-002 /
   data green(10, 2, 0) / 9.82537204688275490D-002 /
   data green(10, 2, 1) / 9.77718934574265813D-002 /
   data green(10, 2, 2) / 9.63706443101250976D-002 /
   data green(10, 3, 0) / 9.59201381331208774D-002 /
   data green(10, 3, 1) / 9.54730677114690451D-002 /
   data green(10, 3, 2) / 9.41706085255491071D-002 /
   data green(10, 3, 3) / 9.21191560357362549D-002 /
   data green(10, 4, 0) / 9.29287075650203770D-002 /
   data green(10, 4, 1) / 9.25233046347111177D-002 /
   data green(10, 4, 2) / 9.13397767965258428D-002 /
   data green(10, 4, 3) / 8.94685167793856928D-002 /
   data green(10, 4, 4) / 8.70384360957417030D-002 /
   data green(10, 5, 0) / 8.94783719879706646D-002 /
   data green(10, 5, 1) / 8.91173446374905337D-002 /
   data green(10, 5, 2) / 8.80610328148980492D-002 /
   data green(10, 5, 3) / 8.63840351013573093D-002 /
   data green(10, 5, 4) / 8.41941648119335778D-002 /
   data green(10, 5, 5) / 8.16145162985881384D-002 /
   data green(10, 6, 0) / 8.57538002584509879D-002 /
   data green(10, 6, 1) / 8.54366087693473608D-002 /
   data green(10, 6, 2) / 8.45065029561923170D-002 /
   data green(10, 6, 3) / 8.30237367744446064D-002 /
   data green(10, 6, 4) / 8.10765375922800485D-002 /
   data green(10, 6, 5) / 7.87674207060537340D-002 /
   data green(10, 6, 6) / 7.62004062654974279D-002 /
   data green(10, 7, 0) / 8.19095980284050701D-002 /
   data green(10, 7, 1) / 8.16335608239375093D-002 /
   data green(10, 7, 2) / 8.08224168301887452D-002 /
   data green(10, 7, 3) / 7.95241043735320369D-002 /
   data green(10, 7, 4) / 7.78096907093466933D-002 /
   data green(10, 7, 5) / 7.57631587064486745D-002 /
   data green(10, 7, 6) / 7.34714735654957429D-002 /
   data green(10, 7, 7) / 7.10167429723932014D-002 /
   data green(10, 8, 0) / 7.80650809701087800D-002 /
   data green(10, 8, 1) / 7.78263403818316157D-002 /
   data green(10, 8, 2) / 7.71234053797715952D-002 /
   data green(10, 8, 3) / 7.59940424924785002D-002 /
   data green(10, 8, 4) / 7.44948826748830367D-002 /
   data green(10, 8, 5) / 7.26939167046615486D-002 /
   data green(10, 8, 6) / 7.06629212424322584D-002 /
   data green(10, 8, 7) / 6.84711727917312885D-002 /
   data green(10, 8, 8) / 6.61810945551142421D-002 /
   data green(10, 9, 0) / 7.43059171048202349D-002 /
   data green(10, 9, 1) / 7.41001561012377546D-002 /
   data green(10, 9, 2) / 7.34932252727740210D-002 /
   data green(10, 9, 3) / 7.25147104938271558D-002 /
   data green(10, 9, 4) / 7.12094391337109167D-002 /
   data green(10, 9, 5) / 6.96320204145258193D-002 /
   data green(10, 9, 6) / 6.78411483482542627D-002 /
   data green(10, 9, 7) / 6.58946553297434512D-002 /
   data green(10, 9, 8) / 6.38458625488757814D-002 /
   data green(10, 9, 9) / 6.17413476380135745D-002 /
   data green(10,10, 0) / 7.06892023220685034D-002 /
   data green(10,10, 1) / 7.05121138934122899D-002 /
   data green(10,10, 2) / 6.99889001464206600D-002 /
   data green(10,10, 3) / 6.91426840752861482D-002 /
   data green(10,10, 4) / 6.80088299964522608D-002 /
   data green(10,10, 5) / 6.66309932908370761D-002 /
   data green(10,10, 6) / 6.50568736152566801D-002 /
   data green(10,10, 7) / 6.33343767639574645D-002 /
   data green(10,10, 8) / 6.15086234019173728D-002 /
   data green(10,10, 9) / 5.96199586748402399D-002 /
   data green(10,10,10) / 5.77029003756548836D-002 /
   data green(11, 0, 0) / 9.11016794560529780D-002 /
   data green(11, 1, 0) / 9.07178080402950282D-002 /
   data green(11, 1, 1) / 9.03389701406010076D-002 /
   data green(11, 2, 0) / 8.95964123274329977D-002 /
   data green(11, 2, 1) / 8.92319733191168918D-002 /
   data green(11, 2, 2) / 8.81662186177034901D-002 /
   data green(11, 3, 0) / 8.78211111646219900D-002 /
   data green(11, 3, 1) / 8.74785827211419231D-002 /
   data green(11, 3, 2) / 8.64757239603848316D-002 /
   data green(11, 3, 3) / 8.48815232253326696D-002 /
   data green(11, 4, 0) / 8.55113553183643738D-002 /
   data green(11, 4, 1) / 8.51958117064051668D-002 /
   data green(11, 4, 2) / 8.42706293739895890D-002 /
   data green(11, 4, 3) / 8.27959471555007320D-002 /
   data green(11, 4, 4) / 8.08596867365886651D-002 /
   data green(11, 5, 0) / 8.28016101081145012D-002 /
   data green(11, 5, 1) / 8.25156789436423538D-002 /
   data green(11, 5, 2) / 8.16759948585757123D-002 /
   data green(11, 5, 3) / 8.03335997087980397D-002 /
   data green(11, 5, 4) / 7.85638333883536405D-002 /
   data green(11, 5, 5) / 7.64552192718954138D-002 /
   data green(11, 6, 0) / 7.98234953648447737D-002 /
   data green(11, 6, 1) / 7.95677354587106950D-002 /
   data green(11, 6, 2) / 7.88154312868782347D-002 /
   data green(11, 6, 3) / 7.76090113717969460D-002 /
   data green(11, 6, 4) / 7.60117111904602405D-002 /
   data green(11, 6, 5) / 7.40988213599695217D-002 /
   data green(11, 6, 6) / 7.19490352411429196D-002 /
   data green(11, 7, 0) / 7.66938502981806686D-002 /
   data green(11, 7, 1) / 7.64672944758421902D-002 /
   data green(11, 7, 2) / 7.57998225978334661D-002 /
   data green(11, 7, 3) / 7.47261595159612046D-002 /
   data green(11, 7, 4) / 7.32985483706119767D-002 /
   data green(11, 7, 5) / 7.15800081327143400D-002 /
   data green(11, 7, 6) / 6.96374596580536381D-002 /
   data green(11, 7, 7) / 6.75359391501195033D-002 /
   data green(11, 8, 0) / 7.35088443846406564D-002 /
   data green(11, 8, 1) / 7.33095415562787650D-002 /
   data green(11, 8, 2) / 7.27214589272932449D-002 /
   data green(11, 8, 3) / 7.17727048416121338D-002 /
   data green(11, 8, 4) / 7.05059462697771033D-002 /
   data green(11, 8, 5) / 6.89732971623870234D-002 /
   data green(11, 8, 6) / 6.72309547614485248D-002 /
   data green(11, 8, 7) / 6.53345065486952425D-002 /
   data green(11, 8, 8) / 6.33354294229500764D-002 /
   data green(11, 9, 0) / 7.03427432071114106D-002 /
   data green(11, 9, 1) / 7.01682108170607055D-002 /
   data green(11, 9, 2) / 6.96524757391646515D-002 /
   data green(11, 9, 3) / 6.88181257625826437D-002 /
   data green(11, 9, 4) / 6.76997330121598961D-002 /
   data green(11, 9, 5) / 6.63400217193870628D-002 /
   data green(11, 9, 6) / 6.47857388397967970D-002 /
   data green(11, 9, 7) / 6.30839101918350131D-002 /
   data green(11, 9, 8) / 6.12789098769009055D-002 /
   data green(11, 9, 9) / 5.94104973158599975D-002 /
   data green(11,10, 0) / 6.72495759679729577D-002 /
   data green(11,10, 1) / 6.70971343447545193D-002 /
   data green(11,10, 2) / 6.66460761821387054D-002 /
   data green(11,10, 3) / 6.59144766605948146D-002 /
   data green(11,10, 4) / 6.49302121725842118D-002 /
   data green(11,10, 5) / 6.37281076521309525D-002 /
   data green(11,10, 6) / 6.23467874701184102D-002 /
   data green(11,10, 7) / 6.08257264535419298D-002 /
   data green(11,10, 8) / 5.92028399374774500D-002 /
   data green(11,10, 9) / 5.75127664977210032D-002 /
   data green(11,10,10) / 5.57858420764695412D-002 /
   data green(11,11, 0) / 6.42662175213827985D-002 /
   data green(11,11, 1) / 6.41332129587133676D-002 /
   data green(11,11, 2) / 6.37391873183614649D-002 /
   data green(11,11, 3) / 6.30985790640633537D-002 /
   data green(11,11, 4) / 6.22338104909962372D-002 /
   data green(11,11, 5) / 6.11731723144518166D-002 /
   data green(11,11, 6) / 5.99484368421452640D-002 /
   data green(11,11, 7) / 5.85925565737739767D-002 /
   data green(11,11, 8) / 5.71377099321229778D-002 /
   data green(11,11, 9) / 5.56138330258502886D-002 /
   data green(11,11,10) / 5.40476651159870064D-002 /
   data green(11,11,11) / 5.24622568815832149D-002 /
   data green(12, 0, 0) / 8.34810544063349369D-002 /
   data green(12, 1, 0) / 8.31864388811280053D-002 /
   data green(12, 1, 1) / 8.28950380959079247D-002 /
   data green(12, 2, 0) / 8.23218876526195936D-002 /
   data green(12, 2, 1) / 8.20397626503952954D-002 /
   data green(12, 2, 2) / 8.12112728370923409D-002 /
   data green(12, 3, 0) / 8.09415957665298585D-002 /
   data green(12, 3, 1) / 8.06738084648534054D-002 /
   data green(12, 3, 2) / 7.98867812985920439D-002 /
   data green(12, 3, 3) / 7.86266410844886077D-002 /
   data green(12, 4, 0) / 7.91250796596190131D-002 /
   data green(12, 4, 1) / 7.88753192465816572D-002 /
   data green(12, 4, 2) / 7.81405217665451096D-002 /
   data green(12, 4, 3) / 7.69617271089097221D-002 /
   data green(12, 4, 4) / 7.54001232360815149D-002 /
   data green(12, 5, 0) / 7.69650800902222931D-002 /
   data green(12, 5, 1) / 7.67355768925574722D-002 /
   data green(12, 5, 2) / 7.60595970232479518D-002 /
   data green(12, 5, 3) / 7.49727763772037592D-002 /
   data green(12, 5, 4) / 7.35286186104458228D-002 /
   data green(12, 5, 5) / 7.17914898361379505D-002 /
   data green(12, 6, 0) / 7.45562406069010392D-002 /
   data green(12, 6, 1) / 7.43479039117091733D-002 /
   data green(12, 6, 2) / 7.37335204781603137D-002 /
   data green(12, 6, 3) / 7.27434297542392111D-002 /
   data green(12, 6, 4) / 7.14235150597491569D-002 /
   data green(12, 6, 5) / 6.98295336175505099D-002 /
   data green(12, 6, 6) / 6.80212431453136868D-002 /
   data green(12, 7, 0) / 7.19865750595089920D-002 /
   data green(12, 7, 1) / 7.17992561201120372D-002 /
   data green(12, 7, 2) / 7.12461745508370797D-002 /
   data green(12, 7, 3) / 7.03527655403806407D-002 /
   data green(12, 7, 4) / 6.91577744054122368D-002 /
   data green(12, 7, 5) / 6.77087635388151782D-002 /
   data green(12, 7, 6) / 6.60573470974006732D-002 /
   data green(12, 7, 7) / 6.42549595148261266D-002 /
   data green(12, 8, 0) / 6.93323087198854937D-002 /
   data green(12, 8, 1) / 6.91651003510040108D-002 /
   data green(12, 8, 2) / 6.86708034490989411D-002 /
   data green(12, 8, 3) / 6.78704960815799396D-002 /
   data green(12, 8, 4) / 6.67965115570084284D-002 /
   data green(12, 8, 5) / 6.54889353863334911D-002 /
   data green(12, 8, 6) / 6.39918051173260782D-002 /
   data green(12, 8, 7) / 6.23496329740869798D-002 /
   data green(12, 8, 8) / 6.06046491313939270D-002 /
   data green(12, 9, 0) / 6.66557001511131203D-002 /
   data green(12, 9, 1) / 6.65072139538171320D-002 /
   data green(12, 9, 2) / 6.60677574731385053D-002 /
   data green(12, 9, 3) / 6.53546538281371825D-002 /
   data green(12, 9, 4) / 6.43946530149174801D-002 /
   data green(12, 9, 5) / 6.32212332353218470D-002 /
   data green(12, 9, 6) / 6.18716105477375503D-002 /
   data green(12, 9, 7) / 6.03839241023316695D-002 /
   data green(12, 9, 8) / 5.87949198444070315D-002 /
   data green(12, 9, 9) / 5.71382835079026069D-002 /
   data green(12,10, 0) / 6.40050240677659116D-002 /
   data green(12,10, 1) / 6.38736162492291037D-002 /
   data green(12,10, 2) / 6.34842828835360251D-002 /
   data green(12,10, 3) / 6.28511827069137968D-002 /
   data green(12,10, 4) / 6.19963152231710890D-002 /
   data green(12,10, 5) / 6.09474589454071652D-002 /
   data green(12,10, 6) / 5.97358414178388100D-002 /
   data green(12,10, 7) / 5.83938877840891948D-002 /
   data green(12,10, 8) / 5.69533034626604970D-002 /
   data green(12,10, 9) / 5.54436279558284048D-002 /
   data green(12,10,10) / 5.38912888557463512D-002 /
   data green(12,11, 0) / 6.14158636374620079D-002 /
   data green(12,11, 1) / 6.12998021115134586D-002 /
   data green(12,11, 2) / 6.09555897615572731D-002 /
   data green(12,11, 3) / 6.03947616269163573D-002 /
   data green(12,11, 4) / 5.96353436111468765D-002 /
   data green(12,11, 5) / 5.87002856246536947D-002 /
   data green(12,11, 6) / 5.76156586122044678D-002 /
   data green(12,11, 7) / 5.64088698806458483D-002 /
   data green(12,11, 8) / 5.51070950590793057D-002 /
   data green(12,11, 9) / 5.37360450557748670D-002 /
   data green(12,11,10) / 5.23191083177746499D-002 /
   data green(12,11,11) / 5.08768486471794512D-002 /
   data green(12,12, 0) / 5.89130257485459066D-002 /
   data green(12,12, 1) / 5.88106036142706581D-002 /
   data green(12,12, 2) / 5.85065602182457054D-002 /
   data green(12,12, 3) / 5.80102801554029726D-002 /
   data green(12,12, 4) / 5.73365062169431394D-002 /
   data green(12,12, 5) / 5.65041519477287382D-002 /
   data green(12,12, 6) / 5.55349125538198915D-002 /
   data green(12,12, 7) / 5.44518593434600504D-002 /
   data green(12,12, 8) / 5.32781693284598082D-002 /
   data green(12,12, 9) / 5.20360886053427910D-002 /
   data green(12,12,10) / 5.07461728062902767D-002 /
   data green(12,12,11) / 4.94268022394494244D-002 /
   data green(12,12,12) / 4.80939386802734029D-002 /
   data green(13, 0, 0) / 7.70388914574446221D-002 /
   data green(13, 1, 0) / 7.68078069957024018D-002 /
   data green(13, 1, 1) / 7.65788553718000087D-002 /
   data green(13, 2, 0) / 7.61273536398360678D-002 /
   data green(13, 2, 1) / 7.59045955327814215D-002 /
   data green(13, 2, 2) / 7.52483305068758213D-002 /
   data green(13, 3, 0) / 7.50338791926117299D-002 /
   data green(13, 3, 1) / 7.48208091058795433D-002 /
   data green(13, 3, 2) / 7.41927196781878989D-002 /
   data green(13, 3, 3) / 7.31813064665954444D-002 /
   data green(13, 4, 0) / 7.35818082894341280D-002 /
   data green(13, 4, 1) / 7.33811169895098175D-002 /
   data green(13, 4, 2) / 7.27890765479496382D-002 /
   data green(13, 4, 3) / 7.18343530336385538D-002 /
   data green(13, 4, 4) / 7.05603612411894177D-002 /
   data green(13, 5, 0) / 7.18363963815843026D-002 /
   data green(13, 5, 1) / 7.16498821399413205D-002 /
   data green(13, 5, 2) / 7.10991917664451495D-002 /
   data green(13, 5, 3) / 7.02096860575817244D-002 /
   data green(13, 5, 4) / 6.90199737651488315D-002 /
   data green(13, 5, 5) / 6.75774152889397073D-002 /
   data green(13, 6, 0) / 6.98664853034117878D-002 /
   data green(13, 6, 1) / 6.96950939138134551D-002 /
   data green(13, 6, 2) / 6.91885863420501762D-002 /
   data green(13, 6, 3) / 6.83689896981443529D-002 /
   data green(13, 6, 4) / 6.72700227803867401D-002 /
   data green(13, 6, 5) / 6.59333687653175515D-002 /
   data green(13, 6, 6) / 6.44046622008368902D-002 /
   data green(13, 7, 0) / 6.77385384918626904D-002 /
   data green(13, 7, 1) / 6.75824883717706382D-002 /
   data green(13, 7, 2) / 6.71208791194369686D-002 /
   data green(13, 7, 3) / 6.63725589342540750D-002 /
   data green(13, 7, 4) / 6.53665423659080813D-002 /
   data green(13, 7, 5) / 6.41389825934796098D-002 /
   data green(13, 7, 6) / 6.27298516100692621D-002 /
   data green(13, 7, 7) / 6.11798577958619222D-002 /
   data green(13, 8, 0) / 6.55125320220421026D-002 /
   data green(13, 8, 1) / 6.53714783068316818D-002 /
   data green(13, 8, 2) / 6.49538333107407273D-002 /
   data green(13, 8, 3) / 6.42755372732422975D-002 /
   data green(13, 8, 4) / 6.33612641840302915D-002 /
   data green(13, 8, 5) / 6.22420017071424053D-002 /
   data green(13, 8, 6) / 6.09523499369491809D-002 /
   data green(13, 8, 7) / 5.95279534299488489D-002 /
   data green(13, 8, 8) / 5.80033597284254354D-002 /
   data green(13, 9, 0) / 6.32396970042815726D-002 /
   data green(13, 9, 1) / 6.31128986813571380D-002 /
   data green(13, 9, 2) / 6.27371149850770721D-002 /
   data green(13, 9, 3) / 6.21257079296907641D-002 /
   data green(13, 9, 4) / 6.12994712337064235D-002 /
   data green(13, 9, 5) / 6.02847197981482980D-002 /
   data green(13, 9, 6) / 5.91111205738101345D-002 /
   data green(13, 9, 7) / 5.78095832718342403D-002 /
   data green(13, 9, 8) / 5.64104489852525828D-002 /
   data green(13, 9, 9) / 5.49421081310847506D-002 /
   data green(13,10, 0) / 6.09617721515890582D-002 /
   data green(13,10, 1) / 6.08482401339875689D-002 /
   data green(13,10, 2) / 6.05114749282386552D-002 /
   data green(13,10, 3) / 5.99626060191272292D-002 /
   data green(13,10, 4) / 5.92190408880568855D-002 /
   data green(13,10, 5) / 5.83029705193266051D-002 /
   data green(13,10, 6) / 5.72396451474761134D-002 /
   data green(13,10, 7) / 5.60556613875163381D-002 /
   data green(13,10, 8) / 5.47774501950815534D-002 /
   data green(13,10, 9) / 5.34300803551650758D-002 /
   data green(13,10,10) / 5.20364185108024221D-002 /
   data green(13,11, 0) / 5.87113102640505341D-002 /
   data green(13,11, 1) / 5.86099264478482720D-002 /
   data green(13,11, 2) / 5.83089444107706828D-002 /
   data green(13,11, 3) / 5.78175944624138294D-002 /
   data green(13,11, 4) / 5.71503824842007485D-002 /
   data green(13,11, 5) / 5.63259283169830552D-002 /
   data green(13,11, 6) / 5.53656053541232568D-002 /
   data green(13,11, 7) / 5.42921619112680873D-002 /
   data green(13,11, 8) / 5.31284726948117667D-002 /
   data green(13,11, 9) / 5.18965173707737065D-002 /
   data green(13,11,10) / 5.06166294256698915D-002 /
   data green(13,11,11) / 4.93070138870452881D-002 /
   data green(13,12, 0) / 5.65126099939050122D-002 /
   data green(13,12, 1) / 5.64222157479026890D-002 /
   data green(13,12, 2) / 5.61536490083429998D-002 /
   data green(13,12, 3) / 5.57145456216981347D-002 /
   data green(13,12, 4) / 5.51169589881166308D-002 /
   data green(13,12, 5) / 5.43764608331509455D-002 /
   data green(13,12, 6) / 5.35110734036607916D-002 /
   data green(13,12, 7) / 5.25401672085944182D-002 /
   data green(13,12, 8) / 5.14834389518986912D-002 /
   data green(13,12, 9) / 5.03600496928784461D-002 /
   data green(13,12,10) / 4.91879646787000344D-002 /
   data green(13,12,11) / 4.79835023278869788D-002 /
   data green(13,12,12) / 4.67610751210973674D-002 /
   data green(13,13, 0) / 5.43829373564185745D-002 /
   data green(13,13, 1) / 5.43023946210955330D-002 /
   data green(13,13, 2) / 5.40629236535351276D-002 /
   data green(13,13, 3) / 5.36708343815004499D-002 /
   data green(13,13, 4) / 5.31361281185516870D-002 /
   data green(13,13, 5) / 5.24718029142820716D-002 /
   data green(13,13, 6) / 5.16930181050038670D-002 /
   data green(13,13, 7) / 5.08162171368132035D-002 /
   data green(13,13, 8) / 4.98582964212330657D-002 /
   data green(13,13, 9) / 4.88358850065051045D-002 /
   data green(13,13,10) / 4.77647725876775772D-002 /
   data green(13,13,11) / 4.66594980060228426D-002 /
   data green(13,13,12) / 4.55330906698840771D-002 /
   data green(13,13,13) / 4.43969445921931077D-002 /
   data green(14, 0, 0) / 7.15210663555232440D-002 /
   data green(14, 1, 0) / 7.13364472522056703D-002 /
   data green(14, 1, 1) / 7.11532890779243560D-002 /
   data green(14, 2, 0) / 7.07913569325197278D-002 /
   data green(14, 2, 1) / 7.06124619219017230D-002 /
   data green(14, 2, 2) / 7.00840797822240014D-002 /
   data green(14, 3, 0) / 6.99108992481112385D-002 /
   data green(14, 3, 1) / 6.97387347773329030D-002 /
   data green(14, 3, 2) / 6.92300141432349025D-002 /
   data green(14, 3, 3) / 6.84070522914311113D-002 /
   data green(14, 4, 0) / 6.87332605237700039D-002 /
   data green(14, 4, 1) / 6.85698066072725904D-002 /
   data green(14, 4, 2) / 6.80865560533712721D-002 /
   data green(14, 4, 3) / 6.73039607796052514D-002 /
   data green(14, 4, 4) / 6.62533906134802203D-002 /
   data green(14, 5, 0) / 6.73052403265549382D-002 /
   data green(14, 5, 1) / 6.71519168957376289D-002 /
   data green(14, 5, 2) / 6.66983224879000719D-002 /
   data green(14, 5, 3) / 6.59628317049096163D-002 /
   data green(14, 5, 4) / 6.49737386123495908D-002 /
   data green(14, 5, 5) / 6.37663168365674576D-002 /
   data green(14, 6, 0) / 6.56775485124957009D-002 /
   data green(14, 6, 1) / 6.55352174139975097D-002 /
   data green(14, 6, 2) / 6.51138430763007620D-002 /
   data green(14, 6, 3) / 6.44296543203541627D-002 /
   data green(14, 6, 4) / 6.35077470280523693D-002 /
   data green(14, 6, 5) / 6.23795961107118171D-002 /
   data green(14, 6, 6) / 6.10802893089339646D-002 /
   data green(14, 7, 0) / 6.39006618226719836D-002 /
   data green(14, 7, 1) / 6.37696842785668572D-002 /
   data green(14, 7, 2) / 6.33816332069690597D-002 /
   data green(14, 7, 3) / 6.27506396084120538D-002 /
   data green(14, 7, 4) / 6.18986512978855347D-002 /
   data green(14, 7, 5) / 6.08533675736250973D-002 /
   data green(14, 7, 6) / 5.96459100358872826D-002 /
   data green(14, 7, 7) / 5.83085771803091776D-002 /
   data green(14, 8, 0) / 6.20216913025595021D-002 /
   data green(14, 8, 1) / 6.19020171040754394D-002 /
   data green(14, 8, 2) / 6.15471867307574305D-002 /
   data green(14, 8, 3) / 6.09693631991578971D-002 /
   data green(14, 8, 4) / 6.01875209671740852D-002 /
   data green(14, 8, 5) / 5.92257583439134705D-002 /
   data green(14, 8, 6) / 5.81113671833983006D-002 /
   data green(14, 8, 7) / 5.68729368725352860D-002 /
   data green(14, 8, 8) / 5.55387044087558013D-002 /
   data green(14, 9, 0) / 6.00823803382308924D-002 /
   data green(14, 9, 1) / 5.99736477041478006D-002 /
   data green(14, 9, 2) / 5.96510175872965870D-002 /
   data green(14, 9, 3) / 5.91248642388543105D-002 /
   data green(14, 9, 4) / 5.84114415873947870D-002 /
   data green(14, 9, 5) / 5.75315215195823376D-002 /
   data green(14, 9, 6) / 5.65088147407153268D-002 /
   data green(14, 9, 7) / 5.53683911521457642D-002 /
   data green(14, 9, 8) / 5.41352727406824580D-002 /
   data green(14, 9, 9) / 5.28333065977705726D-002 /
   data green(14,10, 0) / 5.81181213054611512D-002 /
   data green(14,10, 1) / 5.80197523456631795D-002 /
   data green(14,10, 2) / 5.77276611369333184D-002 /
   data green(14,10, 3) / 5.72506351994494095D-002 /
   data green(14,10, 4) / 5.66024998488093081D-002 /
   data green(14,10, 5) / 5.58010300198879886D-002 /
   data green(14,10, 6) / 5.48666716322414941D-002 /
   data green(14,10, 7) / 5.38212400153487980D-002 /
   data green(14,10, 8) / 5.26867342850830048D-002 /
   data green(14,10, 9) / 5.14843599674195962D-002 /
   data green(14,10,10) / 5.02338026473981847D-002 /
   data green(14,11, 0) / 5.61577629056505476D-002 /
   data green(14,11, 1) / 5.60690460188646295D-002 /
   data green(14,11, 2) / 5.58054318562229978D-002 /
   data green(14,11, 3) / 5.53743266117698307D-002 /
   data green(14,11, 4) / 5.47874286384871304D-002 /
   data green(14,11, 5) / 5.40598650864430666D-002 /
   data green(14,11, 6) / 5.32091647331469422D-002 /
   data green(14,11, 7) / 5.22541949419168816D-002 /
   data green(14,11, 8) / 5.12141726807692002D-002 /
   data green(14,11, 9) / 5.01078270003420881D-002 /
   data green(14,11,10) / 4.89527537889232894D-002 /
   data green(14,11,11) / 4.77649711679783845D-002 /
   data green(14,12, 0) / 5.42239550986096802D-002 /
   data green(14,12, 1) / 5.41441108019536099D-002 /
   data green(14,12, 2) / 5.39067043979055902D-002 /
   data green(14,12, 3) / 5.35179567344118770D-002 /
   data green(14,12, 4) / 5.29877306541550741D-002 /
   data green(14,12, 5) / 5.23288494281343858D-002 /
   data green(14,12, 6) / 5.15562763672723010D-002 /
   data green(14,12, 7) / 5.06862525044656648D-002 /
   data green(14,12, 8) / 4.97354783654146218D-002 /
   data green(14,12, 9) / 4.87204035440502101D-002 /
   data green(14,12,10) / 4.76566612317969804D-002 /
   data green(14,12,11) / 4.65586600350068672D-002 /
   data green(14,12,12) / 4.54393260507686692D-002 /
   data green(14,13, 0) / 5.23338074542551643D-002 /
   data green(14,13, 1) / 5.22620376104070625D-002 /
   data green(14,13, 2) / 5.20485074509643639D-002 /
   data green(14,13, 3) / 5.16984316367597660D-002 /
   data green(14,13, 4) / 5.12201067196471835D-002 /
   data green(14,13, 5) / 5.06243747979174472D-002 /
   data green(14,13, 6) / 4.99239708858031581D-002 /
   data green(14,13, 7) / 4.91328269387267397D-002 /
   data green(14,13, 8) / 4.82653992401413309D-002 /
   data green(14,13, 9) / 4.73360708221815599D-002 /
   data green(14,13,10) / 4.63586615780687064D-002 /
   data green(14,13,11) / 4.53460600614007114D-002 /
   data green(14,13,12) / 4.43099755894468886D-002 /
   data green(14,13,13) / 4.32607985145020899D-002 /
   data green(14,14, 0) / 5.04996884219655121D-002 /
   data green(14,14, 1) / 5.04352107591626297D-002 /
   data green(14,14, 2) / 5.02432655691877131D-002 /
   data green(14,14, 3) / 4.99282202320272747D-002 /
   data green(14,14, 4) / 4.94970459279078234D-002 /
   data green(14,14, 5) / 4.89588962601307937D-002 /
   data green(14,14, 6) / 4.83245893390172945D-002 /
   data green(14,14, 7) / 4.76060480196575098D-002 /
   data green(14,14, 8) / 4.68157496760220199D-002 /
   data green(14,14, 9) / 4.59662269210508179D-002 /
   data green(14,14,10) / 4.50696472427775038D-002 /
   data green(14,14,11) / 4.41374856848144176D-002 /
   data green(14,14,12) / 4.31802927082936064D-002 /
   data green(14,14,13) / 4.22075504415200947D-002 /
   data green(14,14,14) / 4.12276049343636639D-002 /
   data green(15, 0, 0) / 6.67417176668720635D-002 /
   data green(15, 1, 0) / 6.65918750736398529D-002 /
   data green(15, 1, 1) / 6.64430607475313545D-002 /
   data green(15, 2, 0) / 6.61485176111618889D-002 /
   data green(15, 2, 1) / 6.60027154286278867D-002 /
   data green(15, 2, 2) / 6.55711960561948942D-002 /
   data green(15, 3, 0) / 6.54294292851001769D-002 /
   data green(15, 3, 1) / 6.52884176375973235D-002 /
   data green(15, 3, 2) / 6.48709422768622285D-002 /
   data green(15, 3, 3) / 6.41930574363108464D-002 /
   data green(15, 4, 0) / 6.44620020639746077D-002 /
   data green(15, 4, 1) / 6.43272544039138749D-002 /
   data green(15, 4, 2) / 6.39281557542651507D-002 /
   data green(15, 4, 3) / 6.32795809808007759D-002 /
   data green(15, 4, 4) / 6.24045836557821504D-002 /
   data green(15, 5, 0) / 6.32804027500043775D-002 /
   data green(15, 5, 1) / 6.31530323987592318D-002 /
   data green(15, 5, 2) / 6.27755954306358360D-002 /
   data green(15, 5, 3) / 6.21616279716667580D-002 /
   data green(15, 5, 4) / 6.13321717259390220D-002 /
   data green(15, 5, 5) / 6.03138162843421879D-002 /
   data green(15, 6, 0) / 6.19224687943404767D-002 /
   data green(15, 6, 1) / 6.18032178070706123D-002 /
   data green(15, 6, 2) / 6.14496446865335014D-002 /
   data green(15, 6, 3) / 6.08738748063057181D-002 /
   data green(15, 6, 4) / 6.00948190466262522D-002 /
   data green(15, 6, 5) / 5.91364871517743862D-002 /
   data green(15, 6, 6) / 5.80260613595687050D-002 /
   data green(15, 7, 0) / 6.04268240719801897D-002 /
   data green(15, 7, 1) / 6.03160877519164015D-002 /
   data green(15, 7, 2) / 5.99875664283402443D-002 /
   data green(15, 7, 3) / 5.94519766918515172D-002 /
   data green(15, 7, 4) / 5.87260902441414726D-002 /
   data green(15, 7, 5) / 5.78313070598728274D-002 /
   data green(15, 7, 6) / 5.67920073589878588D-002 /
   data green(15, 7, 7) / 5.56339114691479161D-002 /
   data green(15, 8, 0) / 5.88305394242611621D-002 /
   data green(15, 8, 1) / 5.87284143015893156D-002 /
   data green(15, 8, 2) / 5.84252559763048657D-002 /
   data green(15, 8, 3) / 5.79304303765072123D-002 /
   data green(15, 8, 4) / 5.72586471906624747D-002 /
   data green(15, 8, 5) / 5.64287711125585345D-002 /
   data green(15, 8, 6) / 5.54624331703892082D-002 /
   data green(15, 8, 7) / 5.43826278806538582D-002 /
   data green(15, 8, 8) / 5.32124474976739328D-002 /
   data green(15, 9, 0) / 5.71674770087873749D-002 /
   data green(15, 9, 1) / 5.70738201682635324D-002 /
   data green(15, 9, 2) / 5.67956307280094719D-002 /
   data green(15, 9, 3) / 5.63410199301363490D-002 /
   data green(15, 9, 4) / 5.57227717162422503D-002 /
   data green(15, 9, 5) / 5.49573646274341063D-002 /
   data green(15, 9, 6) / 5.40638163831138815D-002 /
   data green(15, 9, 7) / 5.30624992717584620D-002 /
   data green(15, 9, 8) / 5.19740510150869667D-002 /
   data green(15, 9, 9) / 5.08184660393346185D-002 /
   data green(15,10, 0) / 5.54673077459700153D-002 /
   data green(15,10, 1) / 5.53817977745845666D-002 /
   data green(15,10, 2) / 5.51276549047290745D-002 /
   data green(15,10, 3) / 5.47118533412023975D-002 /
   data green(15,10, 4) / 5.41454227595796131D-002 /
   data green(15,10, 5) / 5.34426515129250296D-002 /
   data green(15,10, 6) / 5.26201351368737746D-002 /
   data green(15,10, 7) / 5.16957866931893051D-002 /
   data green(15,10, 8) / 5.06879101111963812D-002 /
   data green(15,10, 9) / 4.96144089187585818D-002 /
   data green(15,10,10) / 4.84921698609009771D-002 /
   data green(15,11, 0) / 5.37550987280526124D-002 /
   data green(15,11, 1) / 5.36772917558491422D-002 /
   data green(15,11, 2) / 5.34459085610939008D-002 /
   data green(15,11, 3) / 5.30669127796336709D-002 /
   data green(15,11, 4) / 5.25497671455447540D-002 /
   data green(15,11, 5) / 5.19067893820188073D-002 /
   data green(15,11, 6) / 5.11523750168845839D-002 /
   data green(15,11, 7) / 5.03021778065025765D-002 /
   data green(15,11, 8) / 4.93723287970468314D-002 /
   data green(15,11, 9) / 4.83787546501165747D-002 /
   data green(15,11,10) / 4.73366312603540895D-002 /
   data green(15,11,11) / 4.62599854729646848D-002 /
   data green(15,12, 0) / 5.20513286339401582D-002 /
   data green(15,12, 1) / 5.19807059380319736D-002 /
   data green(15,12, 2) / 5.17705705355376994D-002 /
   data green(15,12, 3) / 5.14260014841828456D-002 /
   data green(15,12, 4) / 5.09550834161331331D-002 /
   data green(15,12, 5) / 5.03683888298675220D-002 /
   data green(15,12, 6) / 4.96783472484228147D-002 /
   data green(15,12, 7) / 4.88985712048689986D-002 /
   data green(15,12, 8) / 4.80432032711625850D-002 /
   data green(15,12, 9) / 4.71263341516875364D-002 /
   data green(15,12,10) / 4.61615237660327679D-002 /
   data green(15,12,11) / 4.51614393702094510D-002 /
   data green(15,12,12) / 4.41376099346071638D-002 /
   data green(15,13, 0) / 5.03721886277229830D-002 /
   data green(15,13, 1) / 5.03081945750110740D-002 /
   data green(15,13, 2) / 5.01176818379909031D-002 /
   data green(15,13, 3) / 4.98049643020074678D-002 /
   data green(15,13, 4) / 4.93769291825912754D-002 /
   data green(15,13, 5) / 4.88426225795071112D-002 /
   data green(15,13, 6) / 4.82127397668812435D-002 /
   data green(15,13, 7) / 4.74990738454795086D-002 /
   data green(15,13, 8) / 4.67139732172699509D-002 /
   data green(15,13, 9) / 4.58698486432731270D-002 /
   data green(15,13,10) / 4.49787575225778224D-002 /
   data green(15,13,11) / 4.40520794750444977D-002 /
   data green(15,13,12) / 4.31002855318321290D-002 /
   data green(15,13,13) / 4.21327944328854578D-002 /
   data green(15,14, 0) / 4.87300470524787352D-002 /
   data green(15,14, 1) / 4.86721176715041934D-002 /
   data green(15,14, 2) / 4.84995737610037433D-002 /
   data green(15,14, 3) / 4.82160733107742587D-002 /
   data green(15,14, 4) / 4.78274728179972283D-002 /
   data green(15,14, 5) / 4.73414963490123367D-002 /
   data green(15,14, 6) / 4.67673249913545669D-002 /
   data green(15,14, 7) / 4.61151476127261220D-002 /
   data green(15,14, 8) / 4.53957123175705510D-002 /
   data green(15,14, 9) / 4.46199114815710868D-002 /
   data green(15,14,10) / 4.37984238473311735D-002 /
   data green(15,14,11) / 4.29414270175069793D-002 /
   data green(15,14,12) / 4.20583844722666259D-002 /
   data green(15,14,13) / 4.11579039887295522D-002 /
   data green(15,14,14) / 4.02476594272072052D-002 /
   data green(15,15, 0) / 4.71339851868556800D-002 /
   data green(15,15, 1) / 4.70815684943477267D-002 /
   data green(15,15, 2) / 4.69253713128683150D-002 /
   data green(15,15, 3) / 4.66684932346309589D-002 /
   data green(15,15, 4) / 4.63159098164122873D-002 /
   data green(15,15, 5) / 4.58742086648960545D-002 /
   data green(15,15, 6) / 4.53512593822024898D-002 /
   data green(15,15, 7) / 4.47558484979013946D-002 /
   data green(15,15, 8) / 4.40973099877395996D-002 /
   data green(15,15, 9) / 4.33851776895595453D-002 /
   data green(15,15,10) / 4.26288792728814966D-002 /
   data green(15,15,11) / 4.18374839067754134D-002 /
   data green(15,15,12) / 4.10195086649216001D-002 /
   data green(15,15,13) / 4.01827828458965752D-002 /
   data green(15,15,14) / 3.93343651597620539D-002 /
   data green(15,15,15) / 3.84805061597912057D-002 /
   data green(16, 0, 0) / 6.25617393511497255D-002 /
   data green(16, 1, 0) / 6.24384464165367165D-002 /
   data green(16, 1, 1) / 6.23158944027282571D-002 /
   data green(16, 2, 0) / 6.20730135184558857D-002 /
   data green(16, 2, 1) / 6.19526385973800456D-002 /
   data green(16, 2, 2) / 6.15957811439651545D-002 /
   data green(16, 3, 0) / 6.14783215121586860D-002 /
   data green(16, 3, 1) / 6.13614296537467810D-002 /
   data green(16, 3, 2) / 6.10148126335701285D-002 /
   data green(16, 3, 3) / 6.04502458111079582D-002 /
   data green(16, 4, 0) / 6.06744110171452480D-002 /
   data green(16, 4, 1) / 6.05621119616142692D-002 /
   data green(16, 4, 2) / 6.02290053800068614D-002 /
   data green(16, 4, 3) / 5.96861002871744428D-002 /
   data green(16, 4, 4) / 5.89506059077331626D-002 /
   data green(16, 5, 0) / 5.96866405176377512D-002 /
   data green(16, 5, 1) / 5.95798071051930195D-002 /
   data green(16, 5, 2) / 5.92627892846903231D-002 /
   data green(16, 5, 3) / 5.87457131814095937D-002 /
   data green(16, 5, 4) / 5.80444447848670886D-002 /
   data green(16, 5, 5) / 5.71792631140588398D-002 /
   data green(16, 6, 0) / 5.85436018174212910D-002 /
   data green(16, 6, 1) / 5.84428557968764695D-002 /
   data green(16, 6, 2) / 5.81437705705478286D-002 /
   data green(16, 6, 3) / 5.76555261208020614D-002 /
   data green(16, 6, 4) / 5.69925433235445589D-002 /
   data green(16, 6, 5) / 5.61733230495244859D-002 /
   data green(16, 6, 6) / 5.52190898269537614D-002 /
   data green(16, 7, 0) / 5.72750988867747404D-002 /
   data green(16, 7, 1) / 5.71808196239511207D-002 /
   data green(16, 7, 2) / 5.69008002564661483D-002 /
   data green(16, 7, 3) / 5.64432586874023118D-002 /
   data green(16, 7, 4) / 5.58211401779501876D-002 /
   data green(16, 7, 5) / 5.50511189181496080D-002 /
   data green(16, 7, 6) / 5.41524210590611510D-002 /
   data green(16, 7, 7) / 5.31456211552948504D-002 /
   data green(16, 8, 0) / 5.59104155590611565D-002 /
   data green(16, 8, 1) / 5.58227653847144858D-002 /
   data green(16, 8, 2) / 5.55623072926253875D-002 /
   data green(16, 8, 3) / 5.51363192709301356D-002 /
   data green(16, 8, 4) / 5.45562987273512240D-002 /
   data green(16, 8, 5) / 5.38371162830790576D-002 /
   data green(16, 8, 6) / 5.29960091036227876D-002 /
   data green(16, 8, 7) / 5.20515390024182609D-002 /
   data green(16, 8, 8) / 5.10226227858061182D-002 /
   data green(16, 9, 0) / 5.44769942077351008D-002 /
   data green(16, 9, 1) / 5.43959534924626439D-002 /
   data green(16, 9, 2) / 5.41550158013263741D-002 /
   data green(16, 9, 3) / 5.37605691800023683D-002 /
   data green(16, 9, 4) / 5.32227342022461164D-002 /
   data green(16, 9, 5) / 5.25546555201825971D-002 /
   data green(16, 9, 6) / 5.17716513515123627D-002 /
   data green(16, 9, 7) / 5.08903225063609410D-002 /
   data green(16, 9, 8) / 4.99277104298759614D-002 /
   data green(16, 9, 9) / 4.89005697705460068D-002 /
   data green(16,10, 0) / 5.29995538843498551D-002 /
   data green(16,10, 1) / 5.29249597270324071D-002 /
   data green(16,10, 2) / 5.27030778865507851D-002 /
   data green(16,10, 3) / 5.23394744736862330D-002 /
   data green(16,10, 4) / 5.18429929927538835D-002 /
   data green(16,10, 5) / 5.12251671517558460D-002 /
   data green(16,10, 6) / 5.04995095828883606D-002 /
   data green(16,10, 7) / 4.96807578174557826D-002 /
   data green(16,10, 8) / 4.87841509580135904D-002 /
   data green(16,10, 9) / 4.78247928835833294D-002 /
   data green(16,10,10) / 4.68171361597822661D-002 /
   data green(16,11, 0) / 5.14996095710567314D-002 /
   data green(16,11, 1) / 5.14311935144671009D-002 /
   data green(16,11, 2) / 5.12275895345907575D-002 /
   data green(16,11, 3) / 5.08936195340082692D-002 /
   data green(16,11, 4) / 5.04369658974161070D-002 /
   data green(16,11, 5) / 4.98676887004478875D-002 /
   data green(16,11, 6) / 4.91976357943425119D-002 /
   data green(16,11, 7) / 4.84398101851088170D-002 /
   data green(16,11, 8) / 4.76077542625336358D-002 /
   data green(16,11, 9) / 4.67149977725102791D-002 /
   data green(16,11,10) / 4.57746000231862432D-002 /
   data green(16,11,11) / 4.47988004086825875D-002 /
   data green(16,12, 0) / 4.99953174731916916D-002 /
   data green(16,12, 1) / 4.99327393259001023D-002 /
   data green(16,12, 2) / 4.97464209940209137D-002 /
   data green(16,12, 3) / 4.94405211373731904D-002 /
   data green(16,12, 4) / 4.90216832004572584D-002 /
   data green(16,12, 5) / 4.84986409249887965D-002 /
   data green(16,12, 6) / 4.78817323047621771D-002 /
   data green(16,12, 7) / 4.71823725531808624D-002 /
   data green(16,12, 8) / 4.64125338906875573D-002 /
   data green(16,12, 9) / 4.55842710354899192D-002 /
   data green(16,12,10) / 4.47093190595706652D-002 /
   data green(16,12,11) / 4.37987775376006902D-002 /
   data green(16,12,12) / 4.28628837527910908D-002 /
   data green(16,13, 0) / 4.85015594756416923D-002 /
   data green(16,13, 1) / 4.84444357644093493D-002 /
   data green(16,13, 2) / 4.82742804034210771D-002 /
   data green(16,13, 3) / 4.79946683597457513D-002 /
   data green(16,13, 4) / 4.76113252602554937D-002 /
   data green(16,13, 5) / 4.71318065507727521D-002 /
   data green(16,13, 6) / 4.65650990507203086D-002 /
   data green(16,13, 7) / 4.59211843324039068D-002 /
   data green(16,13, 8) / 4.52106019859925615D-002 /
   data green(16,13, 9) / 4.44440446643507570D-002 /
   data green(16,13,10) / 4.36320078257165858D-002 /
   data green(16,13,11) / 4.27845073467398862D-002 /
   data green(16,13,12) / 4.19108692867421473D-002 /
   data green(16,13,13) / 4.10195890416873016D-002 /
   data green(16,14, 0) / 4.70301847353198285D-002 /
   data green(16,14, 1) / 4.69781113870211681D-002 /
   data green(16,14, 2) / 4.68229328636970560D-002 /
   data green(16,14, 3) / 4.65677155169788720D-002 /
   data green(16,14, 4) / 4.62173823212853330D-002 /
   data green(16,14, 5) / 4.57784527801882446D-002 /
   data green(16,14, 6) / 4.52587174982879356D-002 /
   data green(16,14, 7) / 4.46668780043931457D-002 /
   data green(16,14, 8) / 4.40121819300945197D-002 /
   data green(16,14, 9) / 4.33040794706412782D-002 /
   data green(16,14,10) / 4.25519205466721590D-002 /
   data green(16,14,11) / 4.17647047188096782D-002 /
   data green(16,14,12) / 4.09508889223158812D-002 /
   data green(16,14,13) / 4.01182523126012852D-002 /
   data green(16,14,14) / 3.92738133318963761D-002 /
   data green(16,15, 0) / 4.55903400137406059D-002 /
   data green(16,15, 1) / 4.55429095616769614D-002 /
   data green(16,15, 2) / 4.54015093131118347D-002 /
   data green(16,15, 3) / 4.51687658284782995D-002 /
   data green(16,15, 4) / 4.48489056814760853D-002 /
   data green(16,15, 5) / 4.44475450463684210D-002 /
   data green(16,15, 6) / 4.39714245738541676D-002 /
   data green(16,15, 7) / 4.34281131960070413D-002 /
   data green(16,15, 8) / 4.28257045519177537D-002 /
   data green(16,15, 9) / 4.21725269567969976D-002 /
   data green(16,15,10) / 4.14768831731457954D-002 /
   data green(16,15,11) / 4.07468307339176883D-002 /
   data green(16,15,12) / 3.99900081461621593D-002 /
   data green(16,15,13) / 3.92135076366183014D-002 /
   data green(16,15,14) / 3.84237915492350729D-002 /
   data green(16,15,15) / 3.76266471598203506D-002 /
   data green(16,16, 0) / 4.41888368876139151D-002 /
   data green(16,16, 1) / 4.41456508775730436D-002 /
   data green(16,16, 2) / 4.40168548646909727D-002 /
   data green(16,16, 3) / 4.38046973587843322D-002 /
   data green(16,16, 4) / 4.35128043883580493D-002 /
   data green(16,16, 5) / 4.31460094573079406D-002 /
   data green(16,16, 6) / 4.27101378723707195D-002 /
   data green(16,16, 7) / 4.22117636716594150D-002 /
   data green(16,16, 8) / 4.16579577315638883D-002 /
   data green(16,16, 9) / 4.10560438393440877D-002 /
   data green(16,16,10) / 4.04133762123440252D-002 /
   data green(16,16,11) / 3.97371478638243666D-002 /
   data green(16,16,12) / 3.90342350506857408D-002 /
   data green(16,16,13) / 3.83110793161428875D-002 /
   data green(16,16,14) / 3.75736056781821592D-002 /
   data green(16,16,15) / 3.68271734326269598D-002 /
   data green(16,16,16) / 3.60765548058803037D-002 /
   data green(17, 0, 0) / 5.88749327738794237D-002 /
   data green(17, 1, 0) / 5.87722617505182579D-002 /
   data green(17, 1, 1) / 5.86701356410046435D-002 /
   data green(17, 2, 0) / 5.84675183761178621D-002 /
   data green(17, 2, 1) / 5.83669974367562669D-002 /
   data green(17, 2, 2) / 5.80685883456957125D-002 /
   data green(17, 3, 0) / 5.79702159447572155D-002 /
   data green(17, 3, 1) / 5.78722755801602748D-002 /
   data green(17, 3, 2) / 5.75814715357317528D-002 /
   data green(17, 3, 3) / 5.71065915888477604D-002 /
   data green(17, 4, 0) / 5.72952780858718425D-002 /
   data green(17, 4, 1) / 5.72007642096797747D-002 /
   data green(17, 4, 2) / 5.69200623587268473D-002 /
   data green(17, 4, 3) / 5.64614504619883259D-002 /
   data green(17, 4, 4) / 5.58379616465337625D-002 /
   data green(17, 5, 0) / 5.64618133108695172D-002 /
   data green(17, 5, 1) / 5.63714126715800373D-002 /
   data green(17, 5, 2) / 5.61028440403486470D-002 /
   data green(17, 5, 3) / 5.56637906392568466D-002 /
   data green(17, 5, 4) / 5.50663714565172366D-002 /
   data green(17, 5, 5) / 5.43262267312277711D-002 /
   data green(17, 6, 0) / 5.54916911120218112D-002 /
   data green(17, 6, 1) / 5.54059178570748961D-002 /
   data green(17, 6, 2) / 5.51510072779183846D-002 /
   data green(17, 6, 3) / 5.47339961556022880D-002 /
   data green(17, 6, 4) / 5.41660067614972152D-002 /
   data green(17, 6, 5) / 5.34614359626657681D-002 /
   data green(17, 6, 6) / 5.26369893849825340D-002 /
   data green(17, 7, 0) / 5.44081126242835639D-002 /
   data green(17, 7, 1) / 5.43273096669584579D-002 /
   data green(17, 7, 2) / 5.40870789359298582D-002 /
   data green(17, 7, 3) / 5.36937896064137304D-002 /
   data green(17, 7, 4) / 5.31575311671118808D-002 /
   data green(17, 7, 5) / 5.24914058794098484D-002 /
   data green(17, 7, 6) / 5.17106799913579773D-002 /
   data green(17, 7, 7) / 5.08318953396623002D-002 /
   data green(17, 8, 0) / 5.32343300661667976D-002 /
   data green(17, 8, 1) / 5.31586820385972769D-002 /
   data green(17, 8, 2) / 5.29336870660102737D-002 /
   data green(17, 8, 3) / 5.25650513616922957D-002 /
   data green(17, 8, 4) / 5.20618352392393027D-002 /
   data green(17, 8, 5) / 5.14358443166094831D-002 /
   data green(17, 8, 6) / 5.07008937552766367D-002 /
   data green(17, 8, 7) / 4.98720305827710444D-002 /
   data green(17, 8, 8) / 4.89647903525743111D-002 /
   data green(17, 9, 0) / 5.19926118092869233D-002 /
   data green(17, 9, 1) / 5.19221656657038941D-002 /
   data green(17, 9, 2) / 5.17125562283075602D-002 /
   data green(17, 9, 3) / 5.13688514016889738D-002 /
   data green(17, 9, 4) / 5.08991170769443682D-002 /
   data green(17, 9, 5) / 5.03138995050680382D-002 /
   data green(17, 9, 6) / 4.96255950129056225D-002 /
   data green(17, 9, 7) / 4.88477771651447817D-002 /
   data green(17, 9, 8) / 4.79945455434421103D-002 /
   data green(17, 9, 9) / 4.70799459570228826D-002 /
   data green(17,10, 0) / 5.07034923272723925D-002 /
   data green(17,10, 1) / 5.06381818097794492D-002 /
   data green(17,10, 2) / 5.04437728779398367D-002 /
   data green(17,10, 3) / 5.01247339075995765D-002 /
   data green(17,10, 4) / 4.96881934718393439D-002 /
   data green(17,10, 5) / 4.91435046013372762D-002 /
   data green(17,10, 6) / 4.85017102453897744D-002 /
   data green(17,10, 7) / 4.77749669791736747D-002 /
   data green(17,10, 8) / 4.69759802434364451D-002 /
   data green(17,10, 9) / 4.61174937218891828D-002 /
   data green(17,10,10) / 4.52118612842534254D-002 /
   data green(17,11, 0) / 4.93853003118119929D-002 /
   data green(17,11, 1) / 4.93249711723733816D-002 /
   data green(17,11, 2) / 4.91453167174035060D-002 /
   data green(17,11, 3) / 4.88502533813495871D-002 /
   data green(17,11, 4) / 4.84460438112845337D-002 /
   data green(17,11, 5) / 4.79409330521652494D-002 /
   data green(17,11, 6) / 4.73446989759632036D-002 /
   data green(17,11, 7) / 4.66681628766375395D-002 /
   data green(17,11, 8) / 4.59227039718347260D-002 /
   data green(17,11, 9) / 4.51198137801519447D-002 /
   data green(17,11,10) / 4.42707154782175441D-002 /
   data green(17,11,11) / 4.33860618575569248D-002 /
   data green(17,12, 0) / 4.80539285067113423D-002 /
   data green(17,12, 1) / 4.79983622712795363D-002 /
   data green(17,12, 2) / 4.78328250273026764D-002 /
   data green(17,12, 3) / 4.75607332518240486D-002 /
   data green(17,12, 4) / 4.71875625665768433D-002 /
   data green(17,12, 5) / 4.67205459636252027D-002 /
   data green(17,12, 6) / 4.61682982569177683D-002 /
   data green(17,12, 7) / 4.55404034012536554D-002 /
   data green(17,12, 8) / 4.48470002357656897D-002 /
   data green(17,12, 9) / 4.40983966741327199D-002 /
   data green(17,12,10) / 4.33047341693977642D-002 /
   data green(17,12,11) / 4.24757152876104049D-002 /
   data green(17,12,12) / 4.16203989334536714D-002 /
   data green(17,13, 0) / 4.67227944887969593D-002 /
   data green(17,13, 1) / 4.66717294194609542D-002 /
   data green(17,13, 2) / 4.65195426006623988D-002 /
   data green(17,13, 3) / 4.62692035384383860D-002 /
   data green(17,13, 4) / 4.59254819466805159D-002 /
   data green(17,13, 5) / 4.54946987089765095D-002 /
   data green(17,13, 6) / 4.49844138512138225D-002 /
   data green(17,13, 7) / 4.44030805546570401D-002 /
   data green(17,13, 8) / 4.37596938811486355D-002 /
   data green(17,13, 9) / 4.30634590249910423D-002 /
   data green(17,13,10) / 4.23234978143882112D-002 /
   data green(17,13,11) / 4.15486052365256556D-002 /
   data green(17,13,12) / 4.07470611250865880D-002 /
   data green(17,13,13) / 3.99264966185514311D-002 /
   data green(17,14, 0) / 4.54029390960156293D-002 /
   data green(17,14, 1) / 4.53560877255195852D-002 /
   data green(17,14, 2) / 4.52164068006456060D-002 /
   data green(17,14, 3) / 4.49864704178924435D-002 /
   data green(17,14, 4) / 4.46704218708160322D-002 /
   data green(17,14, 5) / 4.42737689154026753D-002 /
   data green(17,14, 6) / 4.38031255892171456D-002 /
   data green(17,14, 7) / 4.32659234568393003D-002 /
   data green(17,14, 8) / 4.26701152500518072D-002 /
   data green(17,14, 9) / 4.20238912443620868D-002 /
   data green(17,14,10) / 4.13354242443884079D-002 /
   data green(17,14,11) / 4.06126537457054543D-002 /
   data green(17,14,12) / 3.98631145973133427D-002 /
   data green(17,14,13) / 3.90938109543602710D-002 /
   data green(17,14,14) / 3.83111328404982798D-002 /
   data green(17,15, 0) / 4.41032141493002683D-002 /
   data green(17,15, 1) / 4.40602772842672782D-002 /
   data green(17,15, 2) / 4.39322214472593758D-002 /
   data green(17,15, 3) / 4.37212738658876721D-002 /
   data green(17,15, 4) / 4.34310266879540732D-002 /
   data green(17,15, 5) / 4.30662691204879711D-002 /
   data green(17,15, 6) / 4.26327744478883774D-002 /
   data green(17,15, 7) / 4.21370598763163434D-002 /
   data green(17,15, 8) / 4.15861375093045035D-002 /
   data green(17,15, 9) / 4.09872730161734630D-002 /
   data green(17,15,10) / 4.03477653162623712D-002 /
   data green(17,15,11) / 3.96747565947886335D-002 /
   data green(17,15,12) / 3.89750778685136884D-002 /
   data green(17,15,13) / 3.82551316496718349D-002 /
   data green(17,15,14) / 3.75208103291311407D-002 /
   data green(17,15,15) / 3.67774468354068046D-002 /
   data green(17,16, 0) / 4.28305197715280750D-002 /
   data green(17,16, 1) / 4.27911971948897393D-002 /
   data green(17,16, 2) / 4.26738811627823539D-002 /
   data green(17,16, 3) / 4.24804966213501514D-002 /
   data green(17,16, 4) / 4.22141540755519087D-002 /
   data green(17,16, 5) / 4.18790121951910232D-002 /
   data green(17,16, 6) / 4.14801024787742889D-002 /
   data green(17,16, 7) / 4.10231300170038754D-002 /
   data green(17,16, 8) / 4.05142648872721389D-002 /
   data green(17,16, 9) / 3.99599375892813191D-002 /
   data green(17,16,10) / 3.93666496123020182D-002 /
   data green(17,16,11) / 3.87408072260405217D-002 /
   data green(17,16,12) / 3.80885834159114187D-002 /
   data green(17,16,13) / 3.74158099473763522D-002 /
   data green(17,16,14) / 3.67278991109264916D-002 /
   data green(17,16,15) / 3.60297928960292835D-002 /
   data green(17,16,16) / 3.53259361791336407D-002 /
   data green(17,17, 0) / 4.15900612598764305D-002 /
   data green(17,17, 1) / 4.15540595443157912D-002 /
   data green(17,17, 2) / 4.14466168589230355D-002 /
   data green(17,17, 3) / 4.12693960673060425D-002 /
   data green(17,17, 4) / 4.10250889832989019D-002 /
   data green(17,17, 5) / 4.07173040107817047D-002 /
   data green(17,17, 6) / 4.03504219727055524D-002 /
   data green(17,17, 7) / 3.99294310957238746D-002 /
   data green(17,17, 8) / 3.94597526538444376D-002 /
   data green(17,17, 9) / 3.89470680828795934D-002 /
   data green(17,17,10) / 3.83971567350649204D-002 /
   data green(17,17,11) / 3.78157512182099448D-002 /
   data green(17,17,12) / 3.72084148302860435D-002 /
   data green(17,17,13) / 3.65804432742990043D-002 /
   data green(17,17,14) / 3.59367908418726997D-002 /
   data green(17,17,15) / 3.52820197081586961D-002 /
   data green(17,17,16) / 3.46202699177196871D-002 /
   data green(17,17,17) / 3.39552470276695897D-002 /
   data green(18, 0, 0) / 5.55988102900537851D-002 /
   data green(18, 1, 0) / 5.55124016545661750D-002 /
   data green(18, 1, 1) / 5.54264010687504985D-002 /
   data green(18, 2, 0) / 5.52556241694632869D-002 /
   data green(18, 2, 1) / 5.51708280801790560D-002 /
   data green(18, 2, 2) / 5.49188109852329423D-002 /
   data green(18, 3, 0) / 5.48356265340743526D-002 /
   data green(18, 3, 1) / 5.47527747002897955D-002 /
   data green(18, 3, 2) / 5.45064987073898208D-002 /
   data green(18, 3, 3) / 5.41034597265383499D-002 /
   data green(18, 4, 0) / 5.42636998230995038D-002 /
   data green(18, 4, 1) / 5.41834446001253356D-002 /
   data green(18, 4, 2) / 5.39448385246057149D-002 /
   data green(18, 4, 3) / 5.35541962513902334D-002 /
   data green(18, 4, 4) / 5.30215201344582596D-002 /
   data green(18, 5, 0) / 5.35544448065255058D-002 /
   data green(18, 5, 1) / 5.34773295063143278D-002 /
   data green(18, 5, 2) / 5.32480020099193943D-002 /
   data green(18, 5, 3) / 5.28723685396749465D-002 /
   data green(18, 5, 4) / 5.23597981757207495D-002 /
   data green(18, 5, 5) / 5.17224824349417633D-002 /
   data green(18, 6, 0) / 5.27247832054066404D-002 /
   data green(18, 6, 1) / 5.26512306170941677D-002 /
   data green(18, 6, 2) / 5.24324361080067483D-002 /
   data green(18, 6, 3) / 5.20738565277253998D-002 /
   data green(18, 6, 4) / 5.15841625035415446D-002 /
   data green(18, 6, 5) / 5.09746639693526016D-002 /
   data green(18, 6, 6) / 5.02586145748939936D-002 /
   data green(18, 7, 0) / 5.17929363468210835D-002 /
   data green(18, 7, 1) / 5.17232469219139832D-002 /
   data green(18, 7, 2) / 5.15158797418120984D-002 /
   data green(18, 7, 3) / 5.11758213307437135D-002 /
   data green(18, 7, 4) / 5.07110093381609381D-002 /
   data green(18, 7, 5) / 5.01318249205643662D-002 /
   data green(18, 7, 6) / 4.94504745297553727D-002 /
   data green(18, 7, 7) / 4.86803297344208866D-002 /
   data green(18, 8, 0) / 5.07774763271745186D-002 /
   data green(18, 8, 1) / 5.07118343820297385D-002 /
   data green(18, 8, 2) / 5.05164465389092432D-002 /
   data green(18, 8, 3) / 5.01958255858712815D-002 /
   data green(18, 8, 4) / 4.97571687857265160D-002 /
   data green(18, 8, 5) / 4.92099152375839663D-002 /
   data green(18, 8, 6) / 4.85652035243928082D-002 /
   data green(18, 8, 7) / 4.78352878815498511D-002 /
   data green(18, 8, 8) / 4.70329670772134528D-002 /
   data green(18, 9, 0) / 4.96965229231394884D-002 /
   data green(18, 9, 1) / 4.96350086157894138D-002 /
   data green(18, 9, 2) / 4.94518445571759946D-002 /
   data green(18, 9, 3) / 4.91510806555876673D-002 /
   data green(18, 9, 4) / 4.87391886415844203D-002 /
   data green(18, 9, 5) / 4.82246805024288752D-002 /
   data green(18, 9, 6) / 4.76176380550421066D-002 /
   data green(18, 9, 7) / 4.69292023662153113D-002 /
   data green(18, 9, 8) / 4.61710691125982176D-002 /
   data green(18, 9, 9) / 4.53550274409287163D-002 /
   data green(18,10, 0) / 4.85671243386267132D-002 /
   data green(18,10, 1) / 4.85097290883547971D-002 /
   data green(18,10, 2) / 4.83387707179981294D-002 /
   data green(18,10, 3) / 4.80578578390934708D-002 /
   data green(18,10, 4) / 4.76727684411163907D-002 /
   data green(18,10, 5) / 4.71911241818053426D-002 /
   data green(18,10, 6) / 4.66219863200704673D-002 /
   data green(18,10, 7) / 4.59754135269777320D-002 /
   data green(18,10, 8) / 4.52620202862653101D-002 /
   data green(18,10, 9) / 4.44925681828869934D-002 /
   data green(18,10,10) / 4.36776131455247907D-002 /
   data green(18,11, 0) / 4.74048291210847073D-002 /
   data green(18,11, 1) / 4.73514724094986558D-002 /
   data green(18,11, 2) / 4.71924883108647553D-002 /
   data green(18,11, 3) / 4.69310730659775119D-002 /
   data green(18,11, 4) / 4.65723544983600302D-002 /
   data green(18,11, 5) / 4.61231162603334063D-002 /
   data green(18,11, 6) / 4.55914536303262952D-002 /
   data green(18,11, 7) / 4.49863937464911601D-002 /
   data green(18,11, 8) / 4.43175124324830458D-002 /
   data green(18,11, 9) / 4.35945750623272163D-002 /
   data green(18,11,10) / 4.28272217547589762D-002 /
   data green(18,11,11) / 4.20247092050030011D-002 /
   data green(18,12, 0) / 4.62234342239082815D-002 /
   data green(18,12, 1) / 4.61739801759287402D-002 /
   data green(18,12, 2) / 4.60265743286255927D-002 /
   data green(18,12, 3) / 4.57840338599035876D-002 /
   data green(18,12, 4) / 4.54508872255894800D-002 /
   data green(18,12, 5) / 4.50331422721259866D-002 /
   data green(18,12, 6) / 4.45379950446884720D-002 /
   data green(18,12, 7) / 4.39735059303630466D-002 /
   data green(18,12, 8) / 4.33482695955368438D-002 /
   data green(18,12, 9) / 4.26711017989998900D-002 /
   data green(18,12,10) / 4.19507607024670187D-002 /
   data green(18,12,11) / 4.11957139911276679D-002 /
   data green(18,12,12) / 4.04139570225369657D-002 /
   data green(18,13, 0) / 4.50348810154911489D-002 /
   data green(18,13, 1) / 4.49891536660978614D-002 /
   data green(18,13, 2) / 4.48528104147056586D-002 /
   data green(18,13, 3) / 4.46283246538251716D-002 /
   data green(18,13, 4) / 4.43196797350197957D-002 /
   data green(18,13, 5) / 4.39321750261640678D-002 /
   data green(18,13, 6) / 4.34721809467941214D-002 /
   data green(18,13, 7) / 4.29468644050828158D-002 /
   data green(18,13, 8) / 4.23639062354336376D-002 /
   data green(18,13, 9) / 4.17312298715637880D-002 /
   data green(18,13,10) / 4.10567563853097561D-002 /
   data green(18,13,11) / 4.03481960996234235D-002 /
   data green(18,13,12) / 3.96128820779332591D-002 /
   data green(18,13,13) / 3.88576465107275115D-002 /
   data green(18,14, 0) / 4.38492657516375353D-002 /
   data green(18,14, 1) / 4.38070623657068128D-002 /
   data green(18,14, 2) / 4.36811857488427421D-002 /
   data green(18,14, 3) / 4.34738009145880580D-002 /
   data green(18,14, 4) / 4.31884009356811524D-002 /
   data green(18,14, 5) / 4.28296454227301213D-002 /
   data green(18,14, 6) / 4.24031553656707277D-002 /
   data green(18,14, 7) / 4.19152814664002921D-002 /
   data green(18,14, 8) / 4.13728634785617763D-002 /
   data green(18,14, 9) / 4.07829964599264738D-002 /
   data green(18,14,10) / 4.01528167989401857D-002 /
   data green(18,14,11) / 3.94893170824069661D-002 /
   data green(18,14,12) / 3.87991949678981526D-002 /
   data green(18,14,13) / 3.80887377075254593D-002 /
   data green(18,14,14) / 3.73637411441104308D-002 /
   data green(18,15, 0) / 4.26749314459827084D-002 /
   data green(18,15, 1) / 4.26360336269577811D-002 /
   data green(18,15, 2) / 4.25199802568614685D-002 /
   data green(18,15, 3) / 4.23286621923916681D-002 /
   data green(18,15, 4) / 4.20651355135047075D-002 /
   data green(18,15, 5) / 4.17334874301232214D-002 /
   data green(18,15, 6) / 4.13386650486802995D-002 /
   data green(18,15, 7) / 4.08862806308549279D-002 /
   data green(18,15, 8) / 4.03824074740954941D-002 /
   data green(18,15, 9) / 3.98333794810334721D-002 /
   data green(18,15,10) / 3.92456052567758881D-002 /
   data green(18,15,11) / 3.86254046782920843D-002 /
   data green(18,15,12) / 3.79788728072346302D-002 /
   data green(18,15,13) / 3.73117731620312074D-002 /
   data green(18,15,14) / 3.66294599890483949D-002 /
   data green(18,15,15) / 3.59368274022995940D-002 /
   data green(18,16, 0) / 4.15186119425983435D-002 /
   data green(18,16, 1) / 4.14827945288718247D-002 /
   data green(18,16, 2) / 4.13758999895808346D-002 /
   data green(18,16, 3) / 4.11995771977885714D-002 /
   data green(18,16, 4) / 4.09564955771591771D-002 /
   data green(18,16, 5) / 4.06502340282422978D-002 /
   data green(18,16, 6) / 4.02851383674862340D-002 /
   data green(18,16, 7) / 3.98661580922771708D-002 /
   data green(18,16, 8) / 3.93986738226348165D-002 /
   data green(18,16, 9) / 3.88883260977166134D-002 /
   data green(18,16,10) / 3.83408545948189725D-002 /
   data green(18,16,11) / 3.77619546512067417D-002 /
   data green(18,16,12) / 3.71571555725661237D-002 /
   data green(18,16,13) / 3.65317229173064903D-002 /
   data green(18,16,14) / 3.58905849729673071D-002 /
   data green(18,16,15) / 3.52382821099231011D-002 /
   data green(18,16,16) / 3.45789366414236132D-002 /
   data green(18,17, 0) / 4.03856044637854183D-002 /
   data green(18,17, 1) / 4.03526423786578978D-002 /
   data green(18,17, 2) / 4.02542416081758114D-002 /
   data green(18,17, 3) / 4.00918386594569728D-002 /
   data green(18,17, 4) / 3.98677628353009336D-002 /
   data green(18,17, 5) / 3.95851443591518359D-002 /
   data green(18,17, 6) / 3.92477958860895718D-002 /
   data green(18,17, 7) / 3.88600759568927326D-002 /
   data green(18,17, 8) / 3.84267434849594192D-002 /
   data green(18,17, 9) / 3.79528119649027648D-002 /
   data green(18,17,10) / 3.74434109423479253D-002 /
   data green(18,17,11) / 3.69036606459137856D-002 /
   data green(18,17,12) / 3.63385638286922302D-002 /
   data green(18,17,13) / 3.57529170394412579D-002 /
   data green(18,17,14) / 3.51512419234627019D-002 /
   data green(18,17,15) / 3.45377358486505873D-002 /
   data green(18,17,16) / 3.39162402061112464D-002 /
   data green(18,17,17) / 3.32902241376194785D-002 /
   data green(18,18, 0) / 3.92799526725313569D-002 /
   data green(18,18, 1) / 3.92496259834842973D-002 /
   data green(18,18, 2) / 3.91590683763494221D-002 /
   data green(18,18, 3) / 3.90095308181411804D-002 /
   data green(18,18, 4) / 3.88030448008398338D-002 /
   data green(18,18, 5) / 3.85423464009018377D-002 /
   data green(18,18, 6) / 3.82307778943979173D-002 /
   data green(18,18, 7) / 3.78721736895126920D-002 /
   data green(18,18, 8) / 3.74707378396375093D-002 /
   data green(18,18, 9) / 3.70309201838776997D-002 /
   data green(18,18,10) / 3.65572973525220488D-002 /
   data green(18,18,11) / 3.60544636573612784D-002 /
   data green(18,18,12) / 3.55269354633846873D-002 /
   data green(18,18,13) / 3.49790711950733330D-002 /
   data green(18,18,14) / 3.44150078127137896D-002 /
   data green(18,18,15) / 3.38386134969194743D-002 /
   data green(18,18,16) / 3.32534554484503475D-002 /
   data green(18,18,17) / 3.26627811499748677D-002 /
   data green(18,18,18) / 3.20695111218862669D-002 /
   data green(19, 0, 0) / 5.26683223481785520D-002 /
   data green(19, 1, 0) / 5.25949115798606953D-002 /
   data green(19, 1, 1) / 5.25218113020078786D-002 /
   data green(19, 2, 0) / 5.23765422916631848D-002 /
   data green(19, 2, 1) / 5.23043601205815112D-002 /
   data green(19, 2, 2) / 5.20896239905642017D-002 /
   data green(19, 3, 0) / 5.20186698665465533D-002 /
   data green(19, 3, 1) / 5.19479746998099054D-002 /
   data green(19, 3, 2) / 5.17376367719403416D-002 /
   data green(19, 3, 3) / 5.13927768528222029D-002 /
   data green(19, 4, 0) / 5.15299603118746299D-002 /
   data green(19, 4, 1) / 5.14612608367627997D-002 /
   data green(19, 4, 2) / 5.12568272200825958D-002 /
   data green(19, 4, 3) / 5.09215401210757063D-002 /
   data green(19, 4, 4) / 5.04631703059938083D-002 /
   data green(19, 5, 0) / 5.09217134871487037D-002 /
   data green(19, 5, 1) / 5.08554423326414359D-002 /
   data green(19, 5, 2) / 5.06581953405659396D-002 /
   data green(19, 5, 3) / 5.03345676340370279D-002 /
   data green(19, 5, 4) / 4.98918839851906284D-002 /
   data green(19, 5, 5) / 4.93397435882760305D-002 /
   data green(19, 6, 0) / 5.02071657328830370D-002 /
   data green(19, 6, 1) / 5.01436701038484520D-002 /
   data green(19, 6, 2) / 4.99546404735710381D-002 /
   data green(19, 6, 3) / 4.96443545287831370D-002 /
   data green(19, 6, 4) / 4.92196402487923498D-002 /
   data green(19, 6, 5) / 4.86894634195081388D-002 /
   data green(19, 6, 6) / 4.80644210661654234D-002 /
   data green(19, 7, 0) / 4.94007520802337285D-002 /
   data green(19, 7, 1) / 4.93402907767683946D-002 /
   data green(19, 7, 2) / 4.91602486153690302D-002 /
   data green(19, 7, 3) / 4.88645671844788326D-002 /
   data green(19, 7, 4) / 4.84595473212775241D-002 /
   data green(19, 7, 5) / 4.79534805691233870D-002 /
   data green(19, 7, 6) / 4.73561944329021212D-002 /
   data green(19, 7, 7) / 4.66785582442524605D-002 /
   data green(19, 8, 0) / 4.85173998628601819D-002 /
   data green(19, 8, 1) / 4.84601458497939327D-002 /
   data green(19, 8, 2) / 4.82896079005559961D-002 /
   data green(19, 8, 3) / 4.80093848425682188D-002 /
   data green(19, 8, 4) / 4.76252386719192011D-002 /
   data green(19, 8, 5) / 4.71447693757817282D-002 /
   data green(19, 8, 6) / 4.65770116871488879D-002 /
   data green(19, 8, 7) / 4.59319940042855965D-002 /
   data green(19, 8, 8) / 4.52202981224102199D-002 /
   data green(19, 9, 0) / 4.75719078321698643D-002 /
   data green(19, 9, 1) / 4.75179550783325078D-002 /
   data green(19, 9, 2) / 4.73572045864638497D-002 /
   data green(19, 9, 3) / 4.70929159081118764D-002 /
   data green(19, 9, 4) / 4.67303163877026556D-002 /
   data green(19, 9, 5) / 4.62763173934893635D-002 /
   data green(19, 9, 6) / 4.57391606042395205D-002 /
   data green(19, 9, 7) / 4.51280284559791162D-002 /
   data green(19, 9, 8) / 4.44526519615243698D-002 /
   data green(19, 9, 9) / 4.37229440975790001D-002 /
   data green(19,10, 0) / 4.65784434835540412D-002 /
   data green(19,10, 1) / 4.65278166384363243D-002 /
   data green(19,10, 2) / 4.63769316345598454D-002 /
   data green(19,10, 3) / 4.61287202462815002D-002 /
   data green(19,10, 4) / 4.57878920140156881D-002 /
   data green(19,10, 5) / 4.53606889655456397D-002 /
   data green(19,10, 6) / 4.48545782808815463D-002 /
   data green(19,10, 7) / 4.42779114636504406D-002 /
   data green(19,10, 8) / 4.36395782192093573D-002 /
   data green(19,10, 9) / 4.29486794403867098D-002 /
   data green(19,10,10) / 4.22142377136032615D-002 /
   data green(19,11, 0) / 4.55501710331639234D-002 /
   data green(19,11, 1) / 4.55028365883855257D-002 /
   data green(19,11, 2) / 4.53617226256850734D-002 /
   data green(19,11, 3) / 4.51294505062936405D-002 /
   data green(19,11, 4) / 4.48102381858589277D-002 /
   data green(19,11, 5) / 4.44096899272174578D-002 /
   data green(19,11, 6) / 4.39345314344110754D-002 /
   data green(19,11, 7) / 4.33923140821592182D-002 /
   data green(19,11, 8) / 4.27911119324429767D-002 /
   data green(19,11, 9) / 4.21392324246824268D-002 /
   data green(19,11,10) / 4.14449569350142516D-002 /
   data green(19,11,11) / 4.07163218806878072D-002 /
   data green(19,12, 0) / 4.44990063483050660D-002 /
   data green(19,12, 1) / 4.44548841561624233D-002 /
   data green(19,12, 2) / 4.43233081830481310D-002 /
   data green(19,12, 3) / 4.41066106335796376D-002 /
   data green(19,12, 4) / 4.38085504215505669D-002 /
   data green(19,12, 5) / 4.34341341123552702D-002 /
   data green(19,12, 6) / 4.29893892316035864D-002 /
   data green(19,12, 7) / 4.24811093891252628D-002 /
   data green(19,12, 8) / 4.19165909401757117D-002 /
   data green(19,12, 9) / 4.13033788879716979D-002 /
   data green(19,12,10) / 4.06490361152080648D-002 /
   data green(19,12,11) / 3.99609456295251658D-002 /
   data green(19,12,12) / 3.92461510636462271D-002 /
   data green(19,13, 0) / 4.34354842964083318D-002 /
   data green(19,13, 1) / 4.33944586052938119D-002 /
   data green(19,13, 2) / 4.32720814901801504D-002 /
   data green(19,13, 3) / 4.30704194602955226D-002 /
   data green(19,13, 4) / 4.27928086221783729D-002 /
   data green(19,13, 5) / 4.24437030713377614D-002 /
   data green(19,13, 6) / 4.20284819879447880D-002 /
   data green(19,13, 7) / 4.15532312968487638D-002 /
   data green(19,13, 8) / 4.10245161807079525D-002 /
   data green(19,13, 9) / 4.04491593247219006D-002 /
   data green(19,13,10) / 3.98340370248758541D-002 /
   data green(19,13,11) / 3.91859018244371943D-002 /
   data green(19,13,12) / 3.85112367417268994D-002 /
   data green(19,13,13) / 3.78161428748962081D-002 /
   data green(19,14, 0) / 4.23687182209220164D-002 /
   data green(19,14, 1) / 4.23306476751853025D-002 /
   data green(19,14, 2) / 4.22170537405488835D-002 /
   data green(19,14, 3) / 4.20297615388951773D-002 /
   data green(19,14, 4) / 4.17717221574281719D-002 /
   data green(19,14, 5) / 4.14468848633388412D-002 /
   data green(19,14, 6) / 4.10600337202023449D-002 /
   data green(19,14, 7) / 4.06166014613921281D-002 /
   data green(19,14, 8) / 4.01224739854629114D-002 /
   data green(19,14, 9) / 3.95837978850975478D-002 /
   data green(19,14,10) / 3.90068013648335946D-002 /
   data green(19,14,11) / 3.83976362039824687D-002 /
   data green(19,14,12) / 3.77622455349752514D-002 /
   data green(19,14,13) / 3.71062595060251618D-002 /
   data green(19,14,14) / 3.64349186310164996D-002 /
   data green(19,15, 0) / 4.13064295784444854D-002 /
   data green(19,15, 1) / 4.12711558800565453D-002 /
   data green(19,15, 2) / 4.11658785361363119D-002 /
   data green(19,15, 3) / 4.09922054057195104D-002 /
   data green(19,15, 4) / 4.07527402577189390D-002 /
   data green(19,15, 5) / 4.04509754470938135D-002 /
   data green(19,15, 6) / 4.00911540500582461D-002 /
   data green(19,15, 7) / 3.96781118273598910D-002 /
   data green(19,15, 8) / 3.92171099221834477D-002 /
   data green(19,15, 9) / 3.87136685815129139D-002 /
   data green(19,15,10) / 3.81734106713082014D-002 /
   data green(19,15,11) / 3.76019216773396531D-002 /
   data green(19,15,12) / 3.70046305941065037D-002 /
   data green(19,15,13) / 3.63867139014003443D-002 /
   data green(19,15,14) / 3.57530229237507033D-002 /
   data green(19,15,15) / 3.51080333804477326D-002 /
   data green(19,16, 0) / 4.02550269165501337D-002 /
   data green(19,16, 1) / 4.02223820405463381D-002 /
   data green(19,16, 2) / 4.01249251830049641D-002 /
   data green(19,16, 3) / 3.99640701467925755D-002 /
   data green(19,16, 4) / 3.97421098125666297D-002 /
   data green(19,16, 5) / 3.94621262403422116D-002 /
   data green(19,16, 6) / 3.91278746708537406D-002 /
   data green(19,16, 7) / 3.87436497587904313D-002 /
   data green(19,16, 8) / 3.83141428994880767D-002 /
   data green(19,16, 9) / 3.78442991336282938D-002 /
   data green(19,16,10) / 3.73391810085645925D-002 /
   data green(19,16,11) / 3.68038451896089483D-002 /
   data green(19,16,12) / 3.62432358150452036D-002 /
   data green(19,16,13) / 3.56620968094566726D-002 /
   data green(19,16,14) / 3.50649037871366004D-002 /
   data green(19,16,15) / 3.44558148981681592D-002 /
   data green(19,16,16) / 3.38386390373393367D-002 /
   data green(19,17, 0) / 3.92197161503906530D-002 /
   data green(19,17, 1) / 3.91895281433141257D-002 /
   data green(19,17, 2) / 3.90993833860866641D-002 /
   data green(19,17, 3) / 3.89505234300292214D-002 /
   data green(19,17, 4) / 3.87449646318988455D-002 /
   data green(19,17, 5) / 3.84854229935946077D-002 /
   data green(19,17, 6) / 3.81752167659030750D-002 /
   data green(19,17, 7) / 3.78181534927936949D-002 /
   data green(19,17, 8) / 3.74184086718442055D-002 /
   data green(19,17, 9) / 3.69804029976353102D-002 /
   data green(19,17,10) / 3.65086843608651013D-002 /
   data green(19,17,11) / 3.60078195776646001D-002 /
   data green(19,17,12) / 3.54822994205614348D-002 /
   data green(19,17,13) / 3.49364590978137463D-002 /
   data green(19,17,14) / 3.43744150251305836D-002 /
   data green(19,17,15) / 3.38000176473282551D-002 /
   data green(19,17,16) / 3.32168192428037301D-002 /
   data green(19,17,17) / 3.26280550858750601D-002 /
   data green(19,18, 0) / 3.82046275703243274D-002 /
   data green(19,18, 1) / 3.81767250473545941D-002 /
   data green(19,18, 2) / 3.80933851200596849D-002 /
   data green(19,18, 3) / 3.79556972063719231D-002 /
   data green(19,18, 4) / 3.77654329576970274D-002 /
   data green(19,18, 5) / 3.75249834959348119D-002 /
   data green(19,18, 6) / 3.72372777518969120D-002 /
   data green(19,18, 7) / 3.69056872446276094D-002 /
   data green(19,18, 8) / 3.65339230972579199D-002 /
   data green(19,18, 9) / 3.61259309906505413D-002 /
   data green(19,18,10) / 3.56857891945987146D-002 /
   data green(19,18,11) / 3.52176139182166642D-002 /
   data green(19,18,12) / 3.47254751352445296D-002 /
   data green(19,18,13) / 3.42133249077294851D-002 /
   data green(19,18,14) / 3.36849391686822389D-002 /
   data green(19,18,15) / 3.31438730115257707D-002 /
   data green(19,18,16) / 3.25934288157926147D-002 /
   data green(19,18,17) / 3.20366360271368139D-002 /
   data green(19,18,18) / 3.14762410937976592D-002 /
   data green(19,19, 0) / 3.72129485359347367D-002 /
   data green(19,19, 1) / 3.71871640385223051D-002 /
   data green(19,19, 2) / 3.71101328172145500D-002 /
   data green(19,19, 3) / 3.69828104558257917D-002 /
   data green(19,19, 4) / 3.68067529349050057D-002 /
   data green(19,19, 5) / 3.65840642526880644D-002 /
   data green(19,19, 6) / 3.63173280119562930D-002 /
   data green(19,19, 7) / 3.60095272379835976D-002 /
   data green(19,19, 8) / 3.56639571003730463D-002 /
   data green(19,19, 9) / 3.52841351922956598D-002 /
   data green(19,19,10) / 3.48737136304644660D-002 /
   data green(19,19,11) / 3.44363965714308273D-002 /
   data green(19,19,12) / 3.39758659048149114D-002 /
   data green(19,19,13) / 3.34957169902523957D-002 /
   data green(19,19,14) / 3.29994054442860271D-002 /
   data green(19,19,15) / 3.24902052238520145D-002 /
   data green(19,19,16) / 3.19711776353601984D-002 /
   data green(19,19,17) / 3.14451504390538331D-002 /
   data green(19,19,18) / 3.09147059141995847D-002 /
   data green(19,19,19) / 3.03821765850236132D-002 /
   data green(20, 0, 0) / 5.00314774795746903D-002 /
   data green(20, 1, 0) / 4.99685805807404890D-002 /
   data green(20, 1, 1) / 4.99059233424123186D-002 /
   data green(20, 2, 0) / 4.97813278929454817D-002 /
   data green(20, 2, 1) / 4.97193803592648129D-002 /
   data green(20, 2, 2) / 4.95349391731085206D-002 /
   data green(20, 3, 0) / 4.94739406620472999D-002 /
   data green(20, 3, 1) / 4.94131459027383410D-002 /
   data green(20, 3, 2) / 4.92321191609733505D-002 /
   data green(20, 3, 3) / 4.89348476043627231D-002 /
   data green(20, 4, 0) / 4.90531570209274195D-002 /
   data green(20, 4, 1) / 4.89939158560428745D-002 /
   data green(20, 4, 2) / 4.88174917255450036D-002 /
   data green(20, 4, 3) / 4.85277024621283831D-002 /
   data green(20, 4, 4) / 4.81306534889718932D-002 /
   data green(20, 5, 0) / 4.85278254063260736D-002 /
   data green(20, 5, 1) / 4.84704847212083859D-002 /
   data green(20, 5, 2) / 4.82996923731441111D-002 /
   data green(20, 5, 3) / 4.80190632889317473D-002 /
   data green(20, 5, 4) / 4.76343839583237977D-002 /
   data green(20, 5, 5) / 4.71532842853168710D-002 /
   data green(20, 6, 0) / 4.79084054168122594D-002 /
   data green(20, 6, 1) / 4.78532506901325344D-002 /
   data green(20, 6, 2) / 4.76889381448528937D-002 /
   data green(20, 6, 3) / 4.74188551040941042D-002 /
   data green(20, 6, 4) / 4.70484297344530705D-002 /
   data green(20, 6, 5) / 4.65848310753389297D-002 /
   data green(20, 6, 6) / 4.60365961172779919D-002 /
   data green(20, 7, 0) / 4.72064289853012867D-002 /
   data green(20, 7, 1) / 4.71536810894511854D-002 /
   data green(20, 7, 2) / 4.69965056150278421D-002 /
   data green(20, 7, 3) / 4.67380464681312682D-002 /
   data green(20, 7, 4) / 4.63833479151904357D-002 /
   data green(20, 7, 5) / 4.59390839447064644D-002 /
   data green(20, 7, 6) / 4.54132205009671400D-002 /
   data green(20, 7, 7) / 4.48146428567184132D-002 /
   data green(20, 8, 0) / 4.64339712379950448D-002 /
   data green(20, 8, 1) / 4.63837870982167499D-002 /
   data green(20, 8, 2) / 4.62342169702316674D-002 /
   data green(20, 8, 3) / 4.59881538044720831D-002 /
   data green(20, 8, 4) / 4.56502453184584947D-002 /
   data green(20, 8, 5) / 4.52266526954254802D-002 /
   data green(20, 8, 6) / 4.47247481812914974D-002 /
   data green(20, 8, 7) / 4.41527796343730225D-002 /
   data green(20, 8, 8) / 4.35195297256279062D-002 /
   data green(20, 9, 0) / 4.56031705668003828D-002 /
   data green(20, 9, 1) / 4.55556469473416270D-002 /
   data green(20, 9, 2) / 4.54139724400468292D-002 /
   data green(20, 9, 3) / 4.51807887300673672D-002 /
   data green(20, 9, 4) / 4.48603456970953288D-002 /
   data green(20, 9, 5) / 4.44582885252377269D-002 /
   data green(20, 9, 6) / 4.39813897528960954D-002 /
   data green(20, 9, 7) / 4.34372503359594445D-002 /
   data green(20, 9, 8) / 4.28339937613703078D-002 /
   data green(20, 9, 9) / 4.21799743407231115D-002 /
   data green(20,10, 0) / 4.47258244204910838D-002 /
   data green(20,10, 1) / 4.46810039574311868D-002 /
   data green(20,10, 2) / 4.45473549924940843D-002 /
   data green(20,10, 3) / 4.43272735756144523D-002 /
   data green(20,10, 4) / 4.40246198575889522D-002 /
   data green(20,10, 5) / 4.36445319958643912D-002 /
   data green(20,10, 6) / 4.31931908973721135D-002 /
   data green(20,10, 7) / 4.26775562166956182D-002 /
   data green(20,10, 8) / 4.21050942309863160D-002 /
   data green(20,10, 9) / 4.14835160043591639D-002 /
   data green(20,10,10) / 4.08205403852928209D-002 /
   data green(20,11, 0) / 4.38130740692686399D-002 /
   data green(20,11, 1) / 4.37709526673668026D-002 /
   data green(20,11, 2) / 4.36453205309585077D-002 /
   data green(20,11, 3) / 4.34383382803791712D-002 /
   data green(20,11, 4) / 4.31534917477160856D-002 /
   data green(20,11, 5) / 4.27954306048003857D-002 /
   data green(20,11, 6) / 4.23697634542783488D-002 /
   data green(20,11, 7) / 4.18828265268343936D-002 /
   data green(20,11, 8) / 4.13414434959480173D-002 /
   data green(20,11, 9) / 4.07526922899516669D-002 /
   data green(20,11,10) / 4.01236917211449384D-002 /
   data green(20,11,11) / 3.94614169500449863D-002 /
   data green(20,12, 0) / 4.28751802240250057D-002 /
   data green(20,12, 1) / 4.28357150360132236D-002 /
   data green(20,12, 2) / 4.27179758640558668D-002 /
   data green(20,12, 3) / 4.25239013703863078D-002 /
   data green(20,12, 4) / 4.22566237497416786D-002 /
   data green(20,12, 5) / 4.19203297502962133D-002 /
   data green(20,12, 6) / 4.15200834210965944D-002 /
   data green(20,12, 7) / 4.10616248536012432D-002 /
   data green(20,12, 8) / 4.05511596552695155D-002 /
   data green(20,12, 9) / 3.99951527240421656D-002 /
   data green(20,12,10) / 3.94001375113943522D-002 /
   data green(20,12,11) / 3.87725489020439920D-002 /
   data green(20,12,12) / 3.81185846168361972D-002 /
   data green(20,13, 0) / 4.19213829831441287D-002 /
   data green(20,13, 1) / 4.18845003477287398D-002 /
   data green(20,13, 2) / 4.17744385371888252D-002 /
   data green(20,13, 3) / 4.15929298231146050D-002 /
   data green(20,13, 4) / 4.13427768874383775D-002 /
   data green(20,13, 5) / 4.10277338160452004D-002 /
   data green(20,13, 6) / 4.06523536608821012D-002 /
   data green(20,13, 7) / 4.02218143568395733D-002 /
   data green(20,13, 8) / 3.97417353016047414D-002 /
   data green(20,13, 9) / 3.92179960981144729D-002 /
   data green(20,13,10) / 3.86565671347445933D-002 /
   data green(20,13,11) / 3.80633592468893056D-002 /
   data green(20,13,12) / 3.74440970744732288D-002 /
   data green(20,13,13) / 3.68042182431456005D-002 /
   data green(20,14, 0) / 4.09598343486711688D-002 /
   data green(20,14, 1) / 4.09254372385836826D-002 /
   data green(20,14, 2) / 4.08227674540535684D-002 /
   data green(20,14, 3) / 4.06533675547909243D-002 /
   data green(20,14, 4) / 4.04197367267565139D-002 /
   data green(20,14, 5) / 4.01252293612408065D-002 /
   data green(20,14, 6) / 3.97739245928093246D-002 /
   data green(20,14, 7) / 3.93704764720785860D-002 /
   data green(20,14, 8) / 3.89199549848345885D-002 /
   data green(20,14, 9) / 3.84276875941274992D-002 /
   data green(20,14,10) / 3.78991096047973028D-002 /
   data green(20,14,11) / 3.73396297399021040D-002 /
   data green(20,14,12) / 3.67545151961122876D-002 /
   data green(20,14,13) / 3.61487983863371595D-002 /
   data green(20,14,14) / 3.55272057824368992D-002 /
   data green(20,15, 0) / 3.99975891270989217D-002 /
   data green(20,15, 1) / 3.99655638230690460D-002 /
   data green(20,15, 2) / 3.98699507506264580D-002 /
   data green(20,15, 3) / 3.97121197623823349D-002 /
   data green(20,15, 4) / 3.94942932100007388D-002 /
   data green(20,15, 5) / 3.92194598409814563D-002 /
   data green(20,15, 6) / 3.88912635861593733D-002 /
   data green(20,15, 7) / 3.85138751408801161D-002 /
   data green(20,15, 8) / 3.80918547651813780D-002 /
   data green(20,15, 9) / 3.76300143958264305D-002 /
   data green(20,15,10) / 3.71332861388225105D-002 /
   data green(20,15,11) / 3.66066027267396496D-002 /
   data green(20,15,12) / 3.60547938286439329D-002 /
   data green(20,15,13) / 3.54825004130318034D-002 /
   data green(20,15,14) / 3.48941078534545507D-002 /
   data green(20,15,15) / 3.42936972365490977D-002 /
   data green(20,16, 0) / 3.90406397467745886D-002 /
   data green(20,16, 1) / 3.90108615914803725D-002 /
   data green(20,16, 2) / 3.89219369988870467D-002 /
   data green(20,16, 3) / 3.87750798516465350D-002 /
   data green(20,16, 4) / 3.85722620214880296D-002 /
   data green(20,16, 5) / 3.83161404897021732D-002 /
   data green(20,16, 6) / 3.80099628425421901D-002 /
   data green(20,16, 7) / 3.76574575699700487D-002 /
   data green(20,16, 8) / 3.72627160878472100D-002 /
   data green(20,16, 9) / 3.68300732168522216D-002 /
   data green(20,16,10) / 3.63639921011580380D-002 /
   data green(20,16,11) / 3.58689584078328572D-002 /
   data green(20,16,12) / 3.53493873039715800D-002 /
   data green(20,16,13) / 3.48095453380376299D-002 /
   data green(20,16,14) / 3.42534880933461697D-002 /
   data green(20,16,15) / 3.36850134268339430D-002 /
   data green(20,16,16) / 3.31076293006686007D-002 /
   data green(20,17, 0) / 3.80939816650557386D-002 /
   data green(20,17, 1) / 3.80663198568486072D-002 /
   data green(20,17, 2) / 3.79836968319361007D-002 /
   data green(20,17, 3) / 3.78471865495683055D-002 /
   data green(20,17, 4) / 3.76585357622046535D-002 /
   data green(20,17, 5) / 3.74201024683368383D-002 /
   data green(20,17, 6) / 3.71347758001899370D-002 /
   data green(20,17, 7) / 3.68058825587040464D-002 /
   data green(20,17, 8) / 3.64370860589307843D-002 /
   data green(20,17, 9) / 3.60322828639209614D-002 /
   data green(20,17,10) / 3.55955024443793955D-002 /
   data green(20,17,11) / 3.51308139308216041D-002 /
   data green(20,17,12) / 3.46422430689082780D-002 /
   data green(20,17,13) / 3.41337013845629950D-002 /
   data green(20,17,14) / 3.36089285263599519D-002 /
   data green(20,17,15) / 3.30714478576906307D-002 /
   data green(20,17,16) / 3.25245346643758537D-002 /
   data green(20,17,17) / 3.19711958377497871D-002 /
   data green(20,18, 0) / 3.71616979683799878D-002 /
   data green(20,18, 1) / 3.71360194284227757D-002 /
   data green(20,18, 2) / 3.70593038869809488D-002 /
   data green(20,18, 3) / 3.69325004564785425D-002 /
   data green(20,18, 4) / 3.67571546762317652D-002 /
   data green(20,18, 5) / 3.65353566188303808D-002 /
   data green(20,18, 6) / 3.62696730985617652D-002 /
   data green(20,18, 7) / 3.59630681983208259D-002 /
   data green(20,18, 8) / 3.56188167364145172D-002 /
   data green(20,18, 9) / 3.52404152782379645D-002 /
   data green(20,18,10) / 3.48314949148735079D-002 /
   data green(20,18,11) / 3.43957393729999866D-002 /
   data green(20,18,12) / 3.39368111967599043D-002 /
   data green(20,18,13) / 3.34582878593106745D-002 /
   data green(20,18,14) / 3.29636088107077424D-002 /
   data green(20,18,15) / 3.24560337165800197D-002 /
   data green(20,18,16) / 3.19386115294787443D-002 /
   data green(20,18,17) / 3.14141595783268140D-002 /
   data green(20,18,18) / 3.08852515582268977D-002 /
   data green(20,19, 0) / 3.62470539989575846D-002 /
   data green(20,19, 1) / 3.62232263916376432D-002 /
   data green(20,19, 2) / 3.61520260844099064D-002 /
   data green(20,19, 3) / 3.60342912850456776D-002 /
   data green(20,19, 4) / 3.58713884927610374D-002 /
   data green(20,19, 5) / 3.56651687886987215D-002 /
   data green(20,19, 6) / 3.54179105386361257D-002 /
   data green(20,19, 7) / 3.51322519131585381D-002 /
   data green(20,19, 8) / 3.48111169887215799D-002 /
   data green(20,19, 9) / 3.44576392208176610D-002 /
   data green(20,19,10) / 3.40750858149314259D-002 /
   data green(20,19,11) / 3.36667860284361325D-002 /
   data green(20,19,12) / 3.32360657983585897D-002 /
   data green(20,19,13) / 3.27861903884772016D-002 /
   data green(20,19,14) / 3.23203160571236062D-002 /
   data green(20,19,15) / 3.18414511202071426D-002 /
   data green(20,19,16) / 3.13524262588350464D-002 /
   data green(20,19,17) / 3.08558735152447763D-002 /
   data green(20,19,18) / 3.03542131367623563D-002 /
   data green(20,19,19) / 2.98496472558476555D-002 /
   data green(20,20, 0) / 3.53525950410734924D-002 /
   data green(20,20, 1) / 3.53304890519835074D-002 /
   data green(20,20, 2) / 3.52644203743242568D-002 /
   data green(20,20, 3) / 3.51551290426699650D-002 /
   data green(20,20, 4) / 3.50038228185520009D-002 /
   data green(20,20, 5) / 3.48121403931247395D-002 /
   data green(20,20, 6) / 3.45821029926544829D-002 /
   data green(20,20, 7) / 3.43160571342382564D-002 /
   data green(20,20, 8) / 3.40166115920775142D-002 /
   data green(20,20, 9) / 3.36865716890995190D-002 /
   data green(20,20,10) / 3.33288738493244877D-002 /
   data green(20,20,11) / 3.29465229807156396D-002 /
   data green(20,20,12) / 3.25425347668048953D-002 /
   data green(20,20,13) / 3.21198843913341331D-002 /
   data green(20,20,14) / 3.16814626608974345D-002 /
   data green(20,20,15) / 3.12300399720686003D-002 /
   data green(20,20,16) / 3.07682381235762047D-002 /
   data green(20,20,17) / 3.02985096183565686D-002 /
   data green(20,20,18) / 2.98231238399709539D-002 /
   data green(20,20,19) / 2.93441593185186853D-002 /
   data green(20,20,20) / 2.88635012119702050D-002 /

end block data blkgreen
