! This file is meant to be included, not compiled
#define ERROR1 12345
#define ERROR2 12346
#define ERROR3 12347

#ifndef DISABLE_NCSU
   use ncsu_sander_hooks, only : &
      ncsu_on_sander_init => on_sander_init, &
      ncsu_on_sander_exit => on_sander_exit
#endif /* DISABLE_NCSU */
   use lmod_driver
   use constants, only : INV_AMBER_ELECTROSTATIC
   ! The main qmmm_struct contains all the QMMM variables and arrays
   use qmmm_module, only : qmmm_nml, qmmm_struct, qmmm_mpi, qm2_params, &
                           qm2_struct, qmmm_vsolv, deallocate_qmmm
   use qmmm_read_and_alloc, only : read_qmmm_nm_and_alloc
   use qmmm_vsolv_module, only: qmmm_vsolv_store_parameters, new
   use qm2_extern_module, only: qm2_extern_finalize
   use sebomd_module, only : sebomd_obj, &
                  sebomd_open_files, sebomd_close_files, sebomd_setup
   use sebomd_arrays, only : init_sebomd_arrays, cleanup_sebomd_arrays
#ifdef LES
   use genbornles
#else
   use genborn
#endif
   use decomp, only : allocate_int_decomp, allocate_real_decomp, nat, nrs, &
                      deallocate_int_decomp, deallocate_real_decomp
   use fastwt
   use les_data, only : temp0les
   use relax_mat
   use nmr, only: nmrrad, impnum
   use ew_recip, only: deallocate_m1m2m3,first_pme
   use parms
   use molecule, only : mol_info, allocate_molecule
   use nblist, only: first_list_flag
   use stack
   use amoeba_runmd, only : AM_RUNMD_get_coords,AM_RUNMD
   use amoeba_mdin, only : beeman_integrator,iamoeba,am_nbead
   use amoeba_interface, only : AMOEBA_deallocate,AMOEBA_readparm
#ifdef RISMSANDER
   use sander_rism_interface, only: rism_setparam, rism_init
#endif /* RISMSANDER */
#ifdef PUPIL_SUPPORT
   use pupildata
#endif /* PUPIL */
#ifdef APBS
   use apbs
#endif /* APBS */
   use xray_interface_module, only: xray_init, xray_read_parm, xray_init_globals
   ! for LIE calculations
   use linear_response, only: ilrt, setup_linear_response, &
                              cleanup_linear_response
#  define rem 0
!RCW+MJW CHARMM SUPPORT
   use charmm_mod, only : charmm_filter_out_qm_atoms
   use memory_module, only: x, ix, ih, memory_init, memory_free

! Self-Guided molecular/Langevin Dynamics (SGLD)
   use sgld, only : isgld,psgld
   
   use nbips, only: ipssys,ips

   use crg_reloc, only: ifcr, cr_backup_charge, cr_allocate, &
                        cr_read_input, cr_check_input

   use emap,only: temap,pemap,qemap

   use file_io_dat
   use constantph, only : cnstph_finalize
   use barostats, only : mcbar_setup
   use random, only: amrset

!AMD
   use amd_mod
!scaledMD
   use scaledMD_mod

   use abfqmmm_module
#ifndef USE_PRMTOP_FILE
   use prmtop_type, only : prmtop_struct
#endif

   implicit none

   ! Input parameters
#ifdef USE_PRMTOP_FILE
   character(len=*) :: prmname
#else
   type(prmtop_struct), intent(in) :: parmdata
#endif
   integer, intent(out) :: ierr
   double precision, dimension(6), intent(in) :: inbox
   double precision, dimension(*), intent(in) :: coordinates
   type(sander_input) :: input_options
   type(qmmm_input_options), optional :: qmmm_options

   logical belly, erstop
#  include "../include/memory.h"
#  include "nmr.h"
#  include "box.h"
#  include "../include/md.h"
#  include "extra.h"
#  include "tgtmd.h"
#  include "multitmd.h"

#  include "parallel.h"
#  include "ew_pme_recip.h"
#  include "ew_frc.h"
#  include "ew_erfc_spline.h"
#ifdef MPI
#  include "ew_parallel.h"
#endif
#  include "ew_mpole.h"
#  include "ew_cntrl.h"
#  include "def_time.h"

   integer native,nr3,nr,ier

   ! nmrcal vars
   _REAL_ f,enmr,devdis,devang,devtor,devplpt,devpln,devgendis,ag,bg,cg
   integer numphi,nhb

   integer, dimension(:), allocatable :: dummy
   integer i

   _REAL_ :: box_center(3)

   ! Make sure that QM/MM system have QM/MM options passed in
   if (input_options%ifqnt == 1 .and. .not. present(qmmm_options)) then
      write(0, *) 'Error: QM/MM setups require QM/MM options to be passed in!'
      stop 1
   end if

   ierr = 0

   ! Make sure the prmtop/inpcrd parsing routines are silenced
   call nxtsec_silence
   call nxtsec_crd_silence


#ifdef USE_PRMTOP_FILE
   parm = trim(prmname)
#endif
   mtmd = 'mtmd' ! needed to prevent sander from thinking mtmd is active
   call xray_init_globals

   ! ==== Flag to tell list builder to print size of list on first call =======
   first_list_flag = .true.
   ! ==== Flag to tell recip space routines to allocate on first call =======
   first_pme = .true.

   ! ==== Initialise first_call flags for QMMM ====
   qmmm_struct%qm_mm_first_call = .true.
   qmmm_struct%fock_first_call = .true.
   qmmm_struct%fock2_2atm_first_call = .true.
   qmmm_struct%qm2_allocate_e_repul_first_call = .true.
   qmmm_struct%qm2_calc_rij_eqns_first_call = .true.
   qmmm_struct%qm2_scf_first_call = .true.
   qmmm_struct%zero_link_charges_first_call = .true.
   qmmm_struct%adj_mm_link_pair_crd_first_call = .true.
   qmmm_struct%num_qmmm_calls = 0

   ! In the single-threaded version, the one process is master
   master = .true.
   erstop = .false.
!  num_direct = 1

   ! --- generic packing scheme ---

   nwdvar = 1
   native = 32
#ifdef ISTAR2

   ! --- Int*2 packing scheme ---

   nwdvar = 2
#endif  /*ISTAR2*/
   numpk = nwdvar
   nbit = native/numpk

   ! ----- Only the master node (only node when single-process)
   !       performs the initial setup and reading/writing -----

   call timer_start(TIME_TOTAL)

   call abfqmmm_init_param()

   masterwork: if (master) then

   if (abfqmmm_param%abfqmmm == 0) then

      ! ---- first, initial reads to determine memory sizes:
      call fill_ucell(inbox(1),inbox(2),inbox(3),inbox(4),inbox(5),inbox(6))
      call api_mdread1(input_options, ierr)
      if (ierr /= 0) return
#ifdef USE_PRMTOP_FILE
      call amopen(8,parm,'O','F','R')
      call rdparm1(8, ierr)
#else
      call cpparm1(parmdata,ierr)
#endif
      if (ierr /= 0) return

      ! Check for illegal input combinations that _would_ cause a fatal error
      ! later in the code. This way we don't quit fatally and give a helpful
      ! error message
      if (ifbox > 0 .and. ntb == 0) then
         write(0, '(a)') 'ntb must be 1 for a periodic system'
         ierr = 1
         return
      end if
      if (ifbox == 0 .and. ntb > 0) then
         write(0,'(a)') 'ntb must be 0 for a non-periodic system'
         ierr = 1
         return
      end if

      if (mtmd /= 'mtmd' .or. itgtmd == 2) call mtmdlx(natom)
      ! --- now, we can allocate memory:

      call locmem()

      ! --- dynamic memory allocation:

      ! GMS: 
      ! Allocate space for module molecule
      ! in the master node
      mol_info%natom = natom
      mol_info%nres  = nres
      call allocate_molecule()

      ! Allocate all global arrays
      if (allocated(ipairs)) deallocate(ipairs)
      allocate( x(lastr), ix(lasti), ipairs(lastpr), ih(lasth), stat = ierr )
      if (ierr /= 0) then
         ierr = 1
         return
      end if
      ix(1:lasti) = 0

      ! This sets up pointer arrays in MEMORY_MODULE to match array-offsets into
      ! the shared X, IX, and IH arrays. Eventually, LOCMEM code should be
      ! merged with MEMORY_MODULE to allocate individual allocatable arrays, but
      ! that will also require updating the MPI code to handle individual
      ! arrays.
      call memory_init()

      ! Allocate the parm arrays
      call allocate_parms()

      if ((igb /= 0 .and. igb /= 10 .and. ipb == 0) &
                    .or.hybridgb>0.or.icnstph.gt.1) &
         call allocate_gb( natom, ncopy )

      if( idecomp > 0 ) then
         nat = natom
         nrs = nres
         call allocate_int_decomp(natom)
      end if

      ! --- finish reading the prmtop file and other user input:
#ifdef USE_PRMTOP_FILE
      call rdparm2(x,ix,ih,8,ierr)
#else
      call cpparm2(x,ix,ih,parmdata,ierr)
#endif

      if (ierr /= 0) goto ERROR1

      call AMOEBA_readparm(8,ntf,ntc,natom,x(lmass))! ntf,ntc get reset if amoeba prmtop
      call xray_read_parm(8,6)

   end if

   if (qmmm_nml%ifqnt .or. abfqmmm_param%abfqmmm == 1) then
      if(abfqmmm_param%abfqmmm == 0) then
         call sebomd_setup
         call read_qmmm_nm_and_alloc(igb, ih, ix, x, cut, use_pme, ntb, 0, &
                                     dummy, 0, .false., qmmm_options)
         if (qmmm_nml%qmtheory%SEBOMD) then
            ! don't do QM/MM
            qmmm_nml%ifqnt= .false.
            sebomd_obj%do_sebomd = .true.
         end if
      end if
      if(qmmm_struct%abfqmmm == 1 .and. abfqmmm_param%abfqmmm == 0) then
         call abfqmmm_setup(natom,nres,ix(i02),ih(m04),ih(m02),x(lmass), &
                            nbonh,nbona,ix(iibh),ix(ijbh),ix(iiba),ix(ijba))
         nr=natom
         x(lcrd:lcrd+natom*3-1) = coordinates(1:natom*3)
         x(lvel:lvel+natom*3-1) = 0.d0
         abfqmmm_param%maxqmstep = nstlim
      end if
      if(abfqmmm_param%abfqmmm == 1) then
         if(abfqmmm_param%system == 1) then
            call abfqmmm_update_qmatoms(x(lcrd))
            if(abfqmmm_param%ntwpdb < 0) then
               call abfqmmm_write_pdb(x(lcrd),ix(i70))
               close(6)
               call mexit(6,1)
            end if
         end if
         call abfqmmm_select_system_qmatoms(natom)
         if(qmmm_nml%ifqnt) then
            call read_qmmm_nm_and_alloc(igb,ih,ix,x,cut,use_pme,ntb,&
                                        abfqmmm_param%qmstep, &
                                        abfqmmm_param%isqm, &
                                        abfqmmm_param%abfcharge,.false., &
                                        qmmm_options)
         endif
      endif
   endif

   if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
      call api_mdread2(x,ix,ih, ierr)
      if (ierr /= 0) goto ERROR2
   endif

#if defined(RISMSANDER) 
      call rism_setparam(mdin,&
           commsander,&
           natom,ntypes,x(L15:L15+natom-1),&
           x(LMASS:LMASS+natom-1),cn1,cn2,&
           ix(i04:i04+ntypes**2-1), ix(i06:i06+natom-1))
#endif /*RISMSANDER*/

      if ( ifcr /= 0 ) then
         call cr_read_input(natom)
         call cr_check_input( ips )
         call cr_backup_charge( x(l15), natom )
      end if

      ! --- alloc memory for decomp module that needs info from mdread2
      if (idecomp == 1 .or. idecomp == 2) then
         call allocate_real_decomp(nrs)
      else if( idecomp == 3 .or. idecomp == 4 ) then
         call allocate_real_decomp(npdec*npdec)
      end if

      ! ----- EVALUATE SOME CONSTANTS FROM MDREAD SETTINGS -----

      nr = nrp
      nr3 = 3*nr
      belly = ibelly > 0

      ! ========================= PUPIL INTERFACE =========================
#ifdef PUPIL_SUPPORT

      ! I moved the PUPIL interface down here so that write() statements work
      ! as advertised. BPR 9/7/09

      ! Initialise the CORBA interface
      puperror = 0
      call fixport()
      call inicorbaintfcmd(puperror)
      if (puperror .ne. 0) then
         write(6,*) 'Error creating PUPIL CORBA interface.'
         call mexit(6,1)
      end if
      pupactive = .true.

      ! Allocation of memory and initialization
      pupStep  = 0
      puperror = 0
      allocate (qcell   (12     ),stat=puperror)
      allocate (pupmask (natom  ),stat=puperror)
      allocate (pupqlist(natom  ),stat=puperror)
      allocate (pupatm  (natom  ),stat=puperror)
      allocate (pupchg  (natom  ),stat=puperror)
      allocate (qfpup   (natom*3),stat=puperror)
      allocate (qcdata  (natom*9),stat=puperror)
      allocate (keyMM   (natom  ),stat=puperror)
      allocate (pupres  (nres   ),stat=puperror)
      allocate (keyres  (nres   ),stat=puperror)

      if (puperror /= 0) then
         write(6,*) 'Error allocating PUPIL interface memory.'
         call mexit(6,1)
      end if

      ! Initialise the "atomic numbers" and "quantum forces" vectors        
      pupqatoms = 0
      iresPup   = 1
      pupres(1) = 1
      do iPup=1,natom
         bs1  = (iPup-1)*3
         call get_atomic_number_pupil(ih(iPup+m06-1),x(lmass+iPup-1),pupatm(iPup))
         if (iresPup .lt. nres) then
            if (iPup .ge. ix(iresPup+i02)) then
               iresPup = iresPup + 1
               pupres(iresPup) = iPup
            end if
         end if
         write (strAux,"(A4,'.',A4)") trim(ih(iresPup+m02-1)),adjustl(ih(iPup+m04-1))
         keyres(iresPup) = trim(ih(iresPup+m02-1))
         keyMM(iPup)     = trim(strAux)

         ! Retrieve the initial charges
         pupchg(iPup) = x(L15+iPup-1)

         do jPup=1,3
            qfpup(bs1+jPup) = 0.0d0
         end do
      end do

      ! Initialise the PUPIL cell
      do iPup=1,12
         qcell(iPup) = 0.0d0
      end do

      ! Submit the KeyMM particles and their respective atomic numbers to PUPIL
      puperror = 0
      call putatomtypes(natom,puperror,pupatm,keyMM)
      if (puperror .ne. 0) then
         write(6,*) 'Error sending MM atom types to PUPIL.'
         call mexit(6,1)
      end if

      puperror = 0
      call putresiduetypes(nres,puperror,pupres,keyres)
      if (puperror .ne. 0) then
         write(6,*) 'Error sending MM residue types to PUPIL.'
         call mexit(6,1)
      end if

#endif /* PUPIL_SUPPORT */
      ! ========================= PUPIL INTERFACE =========================

      ! --- seed the random number generator ---

      ! DAN ROE: Note master node only here
      if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
         call amrset(ig)
         if (ntp > 0.and.iabs(ntb) /= 2) then
            write(6,*) 'Input of NTP/NTB inconsistent'
            goto ERROR3
         end if
      end if

      ! ----- READ COORDINATES AND VELOCITIES -----

      if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then   ! <<< QMSTEP=1 BLOCK >>>

         call timer_start(TIME_RDCRD)
         x(lcrd:lcrd+natom*3-1) = coordinates(1:natom*3)
         x(lvel:lvel+natom*3-1) = 0.d0
         if (iamoeba > 0) then
            natom = natom*am_nbead
            nrp   = nrp*am_nbead
            nr    = nr*am_nbead
            nr3   = nr3*am_nbead
            ncopy = am_nbead
         end if

         ! M-WJ
         !if( igb == 0 .and. induced == 1 ) call get_dips(x,nr)
! WJM  if is a polarizable model, reading input dipole information
         if (igb == 0 .and. ipb == 0 .and. induced > 0) call get_dips(x,nr)

#ifdef APBS
         ! APBS initialization
         if (mdin_apbs) then
            ! in: natom, coords, charge and radii (from prmtop)
            ! out: pb charges and pb radii (via apbs_vars module)
            call apbs_init(natom, x(lcrd), x(l15), x(l97))
         end if
#endif /* APBS */

         call xray_init()

      ! ----- SET THE INITIAL VELOCITIES -----

         if (ntx <= 3) then
            call setvel(nr,x(lvel),x(lwinv),tempi,iscale,scalm)
            ! random numbers may have been "used up" in setting the intial
            ! velocities; re-set the generator so that all nodes are back in
            ! sync

            ! DAN ROE: Note master node only here
            call amrset(ig)

         end if
         if (belly) call bellyf(natom,ix(ibellygp),x(lvel))
         call timer_stop(TIME_RDCRD)

         if(abfqmmm_param%abfqmmm == 1 .and. ntb > 0) then
            call iwrap2(abfqmmm_param%n_user_qm, abfqmmm_param%user_qm, &
                        x(lcrd), box_center)
         end if

         ! --- If we are reading NMR restraints/weight changes,
         !     read them now:

         if (nmropt >= 1) then
            call nmrcal(x(lcrd),f,ih(m04),ih(m02),ix(i02),x(lwinv),enmr, &
                  devdis,devang,devtor,devplpt,devpln,devgendis,temp0,tautp,&
                  cut,ntb,x(lnmr01),ix(inmr02),x(l95),5,6,rk,tk,pk,cn1,cn2, &
                  ag,bg,cg,numbnd,numang,numphi,nimprp, &
                  nhb,natom,natom,ntypes,nres,rad,wel,radhb, &
                  welhb,rwell,isftrp,tgtrmsd,temp0les,-1,'READ')
            ! Updated 9/2007 by Matthew Seetin to enable plane-point and
            ! plane-plane restraints

            ! --- Determine how many of the torsional parameters
            !     are impropers
            call impnum(ix(i46),ix(i56),ix(i48),ix(i58),nphih,nphia, &
                  0,nptra,nimprp)
         end if

         ! -- Set up info related to weight changes for the non-bonds:

         call nmrrad(rad,wel,cn1,cn2,ntypes,0,0.0d0)
         call decnvh(asol,bsol,nphb,radhb,welhb)

         if (iredir(4) > 0) call noeread(x,ix,ih)
         if (iredir(8) > 0) call alignread(natom, x(lcrd))
         if (iredir(9) > 0) call csaread

      end if ! <<< QMSTEP=1 BLOCK >>>

      !---------------------------------------------------------------
      ! --- Call FASTWAT, which will tag those bonds which are part
      !     of 3-point water molecules. Constraints will be effected
      !     for these waters using a fast analytic routine -- dap.

      call timer_start(TIME_FASTWT)

      call fastwat(ih(m04),nres,ix(i02),ih(m02), &
            nbonh,nbona,ix(iibh),ix(ijbh),ibelly,ix(ibellygp), &
            iwtnm,iowtnm,ihwtnm,jfastw,ix(iifstwt), &
            ix(iifstwr),ibgwat,ienwat,ibgion,ienion,iorwat, &
            6,natom)
      call timer_stop(TIME_FASTWT)

      call getwds(ih(m04), nres        , ix(i02)     , ih(m02)     , &
             nbonh       , nbona       , 0           , ix(iibh)    , &
             ix(ijbh)    , iwtnm       , iowtnm      , ihwtnm      , &
             jfastw      , ix(iicbh)   , req         , x(lwinv)    , &
             rbtarg      , ibelly      , ix(ibellygp), 6)

      ! Assign link atoms between quantum mechanical and molecular mechanical
      ! atoms if quantum atoms are present.
      ! After assigning the link atoms, delete all connectivity between the
      ! QM atoms.
      if(qmmm_nml%ifqnt) then

         call identify_link_atoms(nbona,ix(iiba),ix(ijba))

         ! Variable QM solvent:
         ! Store the original bond parameters since we will need to rebuild
         ! the QM region (delete bonded terms etc) repeatedly
         if ( qmmm_nml%vsolv > 0 ) then
            call new(qmmm_vsolv, nbonh, nbona, ntheth, ntheta, nphih, nphia)
            call qmmm_vsolv_store_parameters(qmmm_vsolv, numbnd, &
                 ix(iibh), ix(ijbh), ix(iicbh), &
                 ix(iiba), ix(ijba), ix(iicba), &
                 ix(i24), ix(i26), ix(i28), ix(i30), &
                 ix(i32), ix(i34), ix(i36), ix(i38), &
                 ix(i40), ix(i42), ix(i44), ix(i46), ix(i48), &
                 ix(i50), ix(i52), ix(i54), ix(i56), ix(i58))
         end if

         if( abfqmmm_param%abfqmmm == 1 ) then
            if(abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
               call abfqmmm_allocate_arrays_of_parameters(numbnd, nbonh, &
                              nbona, ntheth, ntheta, nphih, nphia)
               call abfqmmm_store_parameters(&
                                 ix(iibh), ix(ijbh), ix(iicbh), &
                                 ix(iiba), ix(ijba), ix(iicba), &
                                 ix(i24), ix(i26), ix(i28), ix(i30), &
                                 ix(i32), ix(i34), ix(i36), ix(i38), &
                                 ix(i40), ix(i42), ix(i44), ix(i46), ix(i48), &
                                 ix(i50), ix(i52), ix(i54), ix(i56), ix(i58), &
                                 x(l15), rk, req)
            else
               call abfqmmm_set_parameters(&
                                 numbnd, nbonh, nbona, ntheth, ntheta, nphih, &
                                 nphia, ix(iibh), ix(ijbh), ix(iicbh), &
                                 ix(iiba), ix(ijba), ix(iicba), &
                                 ix(i24), ix(i26), ix(i28), ix(i30), &
                                 ix(i32), ix(i34), ix(i36), ix(i38), &
                                 ix(i40), ix(i42), ix(i44), ix(i46), ix(i48), &
                                 ix(i50), ix(i52), ix(i54), ix(i56), ix(i58), &
                                 x(l15), rk, req)

               call init_extra_pts(ix(iibh),ix(ijbh),ix(iicbh), &
                                 ix(iiba),ix(ijba),ix(iicba), &
                                 ix(i24),ix(i26),ix(i28),ix(i30), &
                                 ix(i32),ix(i34),ix(i36),ix(i38), &
                                 ix(i40),ix(i42),ix(i44),ix(i46),ix(i48), &
                                 ix(i50),ix(i52),ix(i54),ix(i56),ix(i58), &
                                 ih(m06),ix,x,ix(i08),ix(i10), &
                                 nspm,ix(i70),x(l75),tmass,tmassinv,&
                                 x(lmass),x(lwinv),req)

            end if
         end if

         ! Remove bonds between QM atoms from list (Hydrogen)
         if (nbonh .gt. 0) &
            call setbon(nbonh,ix(iibh),ix(ijbh),ix(iicbh),ix(ibellygp))

         ! Remove bonds between QM atoms from list (Heavy)
         if (nbona .gt. 0) &
            call setbon(nbona,ix(iiba),ix(ijba),ix(iicba),ix(ibellygp))

         ! Remove angles between QM atoms from list (Hydrogen)
         if (ntheth .gt. 0) &
            call setang(ntheth,ix(i24),ix(i26),ix(i28),ix(i30),ix(ibellygp))

         ! Remove angles between QM atoms from list (Heavy)
         if (ntheta .gt. 0) &
            call setang(ntheta,ix(i32),ix(i34),ix(i36),ix(i38),ix(ibellygp))

         ! Remove dihedrals between QM atoms from list (Hydrogen)
         if (nphih .gt. 0) &
            call setdih(nphih,ix(i40),ix(i42),ix(i44),ix(i46),&
                        ix(i48),ix(ibellygp))

         ! Remove dihedrals between QM atoms from list (Heavy)
         if (nphia .gt. 0) &
            call setdih(nphia,ix(i50),ix(i52),ix(i54),ix(i56),&
                        ix(i58),ix(ibellygp))

         ! Remove CHARMM energy terms from QM region
         call charmm_filter_out_qm_atoms()

         ! Now we should work out the type of each quantum atom present. 
         ! This is used for our arrays of pre-computed parameters. It is 
         ! essentially a re-basing of the atomic numbers and is done to save 
         ! memory. Note: qm_assign_atom_types will allocate the qm_atom_type 
         ! array for us. Only the master calls this routine. All other 
         ! threads get this allocated and broadcast to them by the mpi setup 
         ! routine.
         call qm_assign_atom_types

         ! Set default QMMM MPI parameters - for single cpu operation.
         ! These will get overwritten by qmmm_mpi_setup if MPI is on.
         qmmm_mpi%commqmmm_master = master
         qmmm_mpi%numthreads = 1
         qmmm_mpi%mytaskid = 0
         qmmm_mpi%natom_start = 1
         qmmm_mpi%natom_end = natom
         qmmm_mpi%nquant_nlink_start = 1
         qmmm_mpi%nquant_nlink_end = qmmm_struct%nquant_nlink

         ! Now we know how many link atoms we can allocate the scf_mchg array...
         allocate(qm2_struct%scf_mchg(qmmm_struct%nquant_nlink), stat = ier)
         REQUIRE(ier == 0) ! Deallocated in deallocate qmmm

         ! We can also allocate ewald_memory
         if (qmmm_nml%qm_ewald > 0 ) then
            call allocate_qmewald(natom)
         end if
         if (qmmm_nml%qmgb == 2 ) then
            call allocate_qmgb(qmmm_struct%nquant_nlink)
         end if

         allocate( qmmm_struct%dxyzqm(3, qmmm_struct%nquant_nlink), stat = ier )
         REQUIRE(ier == 0) ! Deallocated in deallocate qmmm

      else if(abfqmmm_param%abfqmmm == 1) then
      
         call abfqmmm_set_parameters(&
                        numbnd, nbonh, nbona, ntheth, ntheta, nphih, nphia, &
                        ix(iibh), ix(ijbh), ix(iicbh), &
                        ix(iiba), ix(ijba), ix(iicba), &
                        ix(i24), ix(i26), ix(i28), ix(i30), &
                        ix(i32), ix(i34), ix(i36), ix(i38), &
                        ix(i40), ix(i42), ix(i44), ix(i46), ix(i48), &
                        ix(i50), ix(i52), ix(i54), ix(i56), ix(i58), &
                        x(l15), rk, req)

         call init_extra_pts(ix(iibh),ix(ijbh),ix(iicbh), &
                           ix(iiba),ix(ijba),ix(iicba), &
                           ix(i24),ix(i26),ix(i28),ix(i30), &
                           ix(i32),ix(i34),ix(i36),ix(i38), &
                           ix(i40),ix(i42),ix(i44),ix(i46),ix(i48), &
                           ix(i50),ix(i52),ix(i54),ix(i56),ix(i58), &
                           ih(m06),ix,x,ix(i08),ix(i10), &
                           nspm,ix(i70),x(l75),tmass,tmassinv,&
                           x(lmass),x(lwinv),req)

      end if !if (qmmm_nml%ifqnt)

      ! --- Open the data dumping files and position it depending
      !     on the type of run:

   ! --- end of master process setup ---
   end if masterwork ! (master)

#  if defined(RISMSANDER)
   call rism_init(commsander)
#  endif /* RISMSANDER */

   if(abfqmmm_param%abfqmmm == 1) then
      if(abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
         allocate(abfqmmm_param%v(3*natom+iscale), stat=ier)
         REQUIRE(ier==0)
         allocate(abfqmmm_param%f(3*natom+iscale), stat=ier)
         REQUIRE(ier==0)
         allocate(abfqmmm_param%f1(3*natom+iscale), stat=ier)
         REQUIRE(ier==0)
         allocate(abfqmmm_param%f2(3*natom+iscale), stat=ier)
         REQUIRE(ier==0)
      end if
   end if

   !   debug needs to copy charges at start and they can't change later
   !   ---------------- Check system is neutral and print warning message ------
   !   ---------------- adjust charges for roundoff error.                ------
   if( igb == 0 .and. ipb == 0 .and. iyammp == 0 ) call check_neutral(x(l15),natom)

   if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
      call amrset(ig+1)
   end if

   if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
      call stack_setup()
   else
      call deallocate_stacks             
      call stack_setup()
   end if

   if (sebomd_obj%do_sebomd) then
      ! open necessary files
      call sebomd_open_files
      ! initialize SEBOMD arrays
      call init_sebomd_arrays(natom)
   endif

#ifdef OPENMP

   ! If -openmp was specified to configure_amber then -DOPENMP is defined and the 
   ! threaded version of MKL will have been linked in. It is important here that
   ! we set the default number of openmp threads for MKL to be 1 to stop conflicts
   ! with threaded vectorization routines when running in parallel etc.
   ! Individual calls to MKL from routines that know what they are doing - e.g.
   ! QMMM calls to diagonalizers etc can increase this limit as long as they
   ! put it back afterwards.
   call omp_set_num_threads(1)

   ! If we are using openmp for matrix diagonalization print some information.
   if (qmmm_nml%ifqnt .and. master) call qm_print_omp_info()
#endif

   ! allocate memory for crg relocation
   if (ifcr /= 0) call cr_allocate( master, natom )

   ! initialize LIE module if used
   if ( ilrt /= 0 ) then
      call setup_linear_response(natom,nres,ih(m04),ih(m06),ix(i02),ih(m02),x(lcrd),x(l15), &
                                 ntypes, ix(i04), ix(i06), cn1, cn2, master)
   end if

   if (igb == 7 .or. igb == 8 ) &
      call igb7_init(natom, x(l97)) !x(l97) is rborn()
     !Hai Nguyen: add igb ==8 here

   if (qmmm_nml%ifqnt) then
      ! Apply charge correction if required.
      if (qmmm_nml%adjust_q>0) then
         call qmmm_adjust_q(qmmm_nml%adjust_q, natom, qmmm_struct%nquant, &
               qmmm_struct%nquant_nlink, qmmm_struct%nlink, x(L15), &
               qmmm_struct%iqmatoms, qmmm_nml%qmcharge, qmmm_struct%atom_mask, &
               qmmm_struct%mm_link_mask, master,x(LCRD), qmmm_nml%vsolv)
      end if
      ! At this point we can also fill the qmmm_struct%scaled_mm_charges
      ! array - we only need to do this once as the charges are constant
      ! during a run. Having a separate array of scaled charges saves us
      ! having to do it on every qmmm routine call. Do this BEFORE zeroing
      ! the QM charges since that routine take care of these values as well.
      do i = 1, natom
         qmmm_struct%scaled_mm_charges(i) = x(L15+(i-1)) * &
                        INV_AMBER_ELECTROSTATIC * qmmm_nml%chg_lambda
                        ! charge scaling factor for FEP
      end do

      ! Zeroing of QM charges MUST be done AFTER call to check_neutral.
      ! Zero out the charges on the quantum mechanical atoms.
      call qm_zero_charges(x(L15),qmmm_struct%scaled_mm_charges,.true.)

      if (qmmm_struct%nlink > 0) then
         ! We need to exclude all electrostatic
         ! interactions with MM link pairs, both QM-MM and MM-MM. Do this by
         ! zeroing the MM link pair charges in the main charge array.
         ! These charges are stored in qmmm_struct%mm_link_pair_resp_charges in case
         ! they are later needed.
         call qm_zero_mm_link_pair_main_chg(qmmm_struct%nlink,&
                     qmmm_struct%link_pairs,x(L15), &
                     qmmm_struct%scaled_mm_charges,.true.)
      end if

   end if

  ! Prepare for SGLD simulation
   if (isgld > 0) call psgld(natom,x(lmass),x(lvel), 0)

  ! Prepare for EMAP constraints
   if (temap) call pemap(dt,temp0,x,ix,ih)

      ! Prepare for Isotropic periodic sum of nonbonded interaction
   if (ips .gt. 0) call ipssys(natom,ntypes,ntb,x(l15), &
            cut,cn1,cn2,ix(i04),ix(i06),x(lcrd))

   ! Set up the MC barostat if requested
   if (ntp > 0 .and. barostat == 2) call mcbar_setup(ig)

   is_setup_ = .true.

   return

ERROR3 continue
   if (idecomp > 0) &
      call deallocate_real_decomp

ERROR2 continue
   if(qmmm_nml%ifqnt) then
      call deallocate_qmmm(qmmm_nml, qmmm_struct, qmmm_vsolv, qm2_params)
      call get_qm2_forces_reset
   end if

ERROR1 continue
   ! Make sure we deallocate what we've allocated so far and bail
   ierr = 1
   is_setup_ = .false.
   call memory_free
   call clean_parms
   if ((igb /= 0 .and. igb /= 10 .and. ipb == 0) &
                 .or. hybridgb>0 .or. icnstph>1) &
      call deallocate_gb
   if (idecomp > 0) &
      call deallocate_int_decomp
   return

#undef ERROR1
#undef ERROR2
#undef ERROR3
