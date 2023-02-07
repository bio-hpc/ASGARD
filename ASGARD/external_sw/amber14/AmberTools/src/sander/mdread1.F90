! This file is _not_ meant to be compiled.  It is meant to be included.
#ifdef API
#  define FATAL_ERROR ierr = 1; return
#else
#  define FATAL_ERROR call mexit(6, 1)
#endif /* API */

   use file_io_dat
   use lmod_driver, only : read_lmod_namelist
   use qmmm_module, only : qmmm_nml, qm_gb
   use constants, only : RETIRED_INPUT_OPTION, zero, one, two, three, seven, &
                         eight, NO_INPUT_VALUE_FLOAT, NO_INPUT_VALUE
   use constantph, only : mccycles
   use amoeba_mdin, only: AMOEBA_read_mdin, iamoeba
   use nose_hoover_module, only: nchain  ! APJ
   use lscivr_vars, only: ilscivr, icorf_lsc
   use les_data, only : temp0les
   use pimd_vars, only: ipimd,itimass
   use neb_vars, only: ineb
   use cmd_vars, only: restart_cmd, eq_cmd, adiab_param
   use stack, only: lastist,lastrst
   use nmr, only: echoin
   use crg_reloc, only: ifcr, cropt, crcut, crskin, crin, crprintcharges
   use sgld, only : isgld, isgsta,isgend,fixcom, &
                    tsgavg,sgft,sgff,sgfd,tempsg,treflf,tsgavp
   use amd_mod, only: iamd,iamdlag,EthreshD,alphaD,EthreshP,alphaP, &
        w_amd,EthreshD_w,alphaD_w,EthreshP_w,alphaP_w
   use scaledMD_mod, only: scaledMD,scaledMD_lambda
   use nbips, only: ips,teips,tvips,teaips,tvaips,raips,mipsx,mipsy,mipsz, &
                    mipso,gridips,dvbips
   use emap,only: temap,gammamap
#ifdef DSSP
   use dssp, only: idssp
#endif /* DSSP */

   use emil_mod,          only : emil_do_calc
   use mdin_emil_dat_mod, only : error_hdr

#if !defined(DISABLE_NCSU) && defined(MPI)
   use ncsu_sander_hooks, only : ncsu_on_mdread1 => on_mdread1
#endif /* ! DISABLE_NCSU && MPI */
#ifndef API
   use xray_interface_module, only: xray_active, xray_read_mdin
#endif /* API */
#ifdef MPI /* SOFT CORE */
   use softcore, only : scalpha,scbeta,ifsc,scmask,logdvdl,dvdl_norest,dynlmb, &
                        sceeorder, tishake, emil_sc
   use mbar, only : ifmbar, bar_intervall, bar_l_min, bar_l_max, bar_l_incr
   use remd, only  : rem
#endif /* MPI */
   ! Parameter for LIE module
   use linear_response, only: ilrt, lrt_interval, lrtmask
#ifdef RISMSANDER
#  ifndef API
   use sander_rism_interface, only: xvvfile, guvfile, huvfile, cuvfile,&
        uuvfile, asympfile, quvFile, chgDistFile, exchemfile, solvenefile, &
        entropyfile, exchemGFfile, solveneGFfile, entropyGFfile, exchemUCfile, &
        solveneUCfile, entropyUCfile, potUVfile
#  endif /* API */
   use sander_rism_interface, only: rismprm
#endif /*RISMSANDER*/
#ifdef APBS
   use apbs
#endif /* APBS */
   use sebomd_module, only: read_sebomd_namelist, sebomd_namelist_default
   implicit none
#  include "box.h"
#  include "def_time.h"
#  include "ew_cntrl.h"
#  include "ew_pme_recip.h"
#  include "../include/md.h"
#  include "../include/memory.h"
#  include "nmr.h"
#  include "tgtmd.h"
#  include "multitmd.h"
#  include "ew_erfc_spline.h"

   character(len=4) watdef(4),watnam,owtnm,hwtnm1,hwtnm2

   _REAL_      dele
   integer     ierr
   integer     imcdo
   integer     itotst
   integer     inerr
   logical     mdin_cntrl, mdin_lmod, mdin_qmmm  ! true if namelists are in mdin
   logical     mdin_sebomd
   integer :: ifqnt    ! local here --> put into qmmm_nml%ifqnt after read here
   integer     mxgrp
   integer     iemap
   _REAL_      dtemp  ! retired 
   _REAL_      dxm  ! retired 
   _REAL_      heat  ! retired 
   _REAL_      timlim ! retired

#ifdef API
   ! Input options passed in to the setup API routine
   type(sander_input), intent(in) :: input_options
#else
   integer     ifind
   character(len=8) date
   character(len=10) time
   character(len=512) :: char_tmp_512
#endif /* API */

#ifdef RISMSANDER
   integer irism
#endif /*RISMSANDER*/

   namelist /cntrl/ irest,ibelly, &
         ntx,ntxo,ntcx,ig,tempi, &
         ntb,ntt,nchain,temp0,tautp, &
         ntp,pres0,comp,taup,barostat,mcbarint, &
         nscm,nstlim,t,dt, &
         ntc,ntcc,nconp,tol,ntf,ntn,nsnb, &
         cut,dielc, &
         ntpr,ntwx,ntwv,ntwe,ntwf,ntave,ntpp,ioutfm, &
         ntr,nrc,ntrx,taur,nmropt, &
         ivcap,cutcap,xcap,ycap,zcap,fcap, &
         xlorth,ylorth,zlorth,xorth,yorth,zorth,forth, &
         imin,drms,dele,dx0, &
         pencut,ipnlty,iscale,scalm,noeskp, &
         maxcyc,ncyc,ntmin,vlimit, &
         mxsub,ipol,jfastw,watnam,owtnm,hwtnm1,hwtnm2, iesp, &
         skmin, skmax, vv,vfac, tmode, ips, &
         mipsx,mipsy,mipsz,mipso,gridips,raips,dvbips, &
         iamd,iamdlag,EthreshD,alphaD,EthreshP,alphaP, &
         w_amd,EthreshD_w,alphaD_w,EthreshP_w,alphaP_w, &
         scaledMD,scaledMD_lambda, &
         iemap,gammamap, &
         isgld,isgsta,isgend,fixcom,tsgavg,sgft,sgff,sgfd,tempsg,treflf,tsgavp,&
         jar, iamoeba, &
         numexchg, repcrd, numwatkeep, hybridgb, &
         ntwprt,tausw, &
         ntwr,iyammp,imcdo, &
         plumed,plumedfile, &
         igb,alpb,Arad,rgbmax,saltcon,offset,gbsa,vrand, &
         surften,iwrap,nrespa,nrespai,gamma_ln,extdiel,intdiel, &
         cut_inner,icfe,clambda,klambda, rbornstat,lastrst,lastist,  &
         itgtmd,tgtrmsd,tgtmdfrc,tgtfitmask,tgtrmsmask, dec_verbose, &
         idecomp,temp0les,restraintmask,restraint_wt,bellymask, &
         noshakemask,crgmask, iwrap_mask, &
         rdt,icnstph,solvph,ntcnstph,ntrelax, mccycles, &
         ifqnt,ievb, ipimd, itimass, ineb,profile_mpi, ilscivr, icorf_lsc, &
         ipb, inp, nkija, idistr, &
         gbneckscale, & 
         gbalphaH,gbbetaH,gbgammaH, &
         gbalphaC,gbbetaC,gbgammaC, &
         gbalphaN,gbbetaN,gbgammaN, &
         gbalphaOS,gbbetaOS,gbgammaOS, &
         gbalphaP,gbbetaP,gbgammaP, &
         Sh,Sc,Sn,So,Ss,Sp, &
         lj1264, &
         ifcr, cropt, crcut, crskin, crin, crprintcharges, &
         csurften, ninterface, gamma_ten, &
#ifdef MPI /* SOFT CORE */
         scalpha, scbeta, ifsc, scmask, logdvdl, dvdl_norest, dynlmb, &
         sceeorder, &
         ifmbar, bar_intervall, bar_l_min, bar_l_max, bar_l_incr, tishake, &
         emil_sc, &
#endif
         ilrt, lrt_interval, lrtmask, &
#ifdef DSSP
         idssp, &
#endif
#ifdef RISMSANDER
         irism,&
#endif /*RISMSANDER*/
         emil_do_calc, & 
         restart_cmd, eq_cmd, adiab_param,  &
         vdwmodel, & ! mjhsieh - the model used for van der Waals
         dtemp, heat, timlim  !all retired 

   ! Define default water residue name and the names of water oxygen & hydrogens
   
   data watdef/'WAT ','O   ','H1  ','H2  '/
   
   !     ----- READ THE CONTROL DATA AND OPEN DIFFERENT FILES -----
   
#ifndef API /* If this is NOT the API */
   if (mdout /= "stdout" ) &
         call amopen(6,mdout,owrite,'F','W')
   call amopen(5,mdin,'O','F','R')
   write(6,9308)
   call date_and_time( DATE=date, TIME=time )
   write(6,'(12(a))') '| Run on ', date(5:6), '/', date(7:8), '/',  &
        date(1:4), ' at ', time(1:2), ':', time(3:4), ':', time(5:6)

   ! Write the path of the current executable and working directory
   call get_command_argument(0, char_tmp_512)
   write(6,'(/,a,a)') '|   Executable path: ', trim(char_tmp_512)
   call getcwd(char_tmp_512)
   write(6,'(a,a)') '| Working directory: ', trim(char_tmp_512)
! Write the hostname if we can get it from environment variable
! Note: get_environment_variable is part of the F2003 standard but seems
!       to be supported by GNU, Intel, IBM and Portland (2010+) compilers
   call get_environment_variable("HOSTNAME", char_tmp_512, inerr)
   if (inerr .eq. 0) then
     write(6,'(a,a,/)') '|          Hostname: Unknown'
   else
     write(6,'(a,a,/)') '|          Hostname: ', trim(char_tmp_512)
   end if

   if (owrite /= 'N') write(6, '(2x,a)') '[-O]verwriting output'

   ! Echo the file assignments to the user:
   
   write(6,9700) 'MDIN'   ,mdin(1:70)  , 'MDOUT' ,mdout(1:70) , &
         'INPCRD' ,inpcrd(1:70), 'PARM'  ,parm(1:70)  , &
         'RESTRT',restrt(1:70) , 'REFC'  ,refc(1:70)  , &
         'MDVEL' ,mdvel(1:70)  , 'MDFRC' ,mdfrc(1:70) , &
         'MDEN'   ,mden(1:70)  , &
         'MDCRD' ,mdcrd(1:70)  , 'MDINFO' ,mdinfo(1:70), &
         'MTMD'  ,mtmd(1:70)   , 'INPDIP', inpdip(1:70), &
         'RSTDIP', rstdip(1:70), 'INPTRAJ', inptraj(1:70)
#  ifdef MPI
   write(6,9702) 'REMLOG',     trim(remlog), &
                 'REMTYPE',    trim(remtype), &
                 'REMSTRIP',   trim(remstripcoord), &
                 'SAVEENE',    trim(saveenefile), &
                 'CLUSTERINF', trim(clusterinfofile), &
                 'RESERVOIR',  trim(reservoirname), &
                 'REMDDIM',    trim(remd_dimension_file)
#  endif
#ifdef RISMSANDER
   if(len_trim(xvvfile) > 0)&
        write(6,9701) 'Xvv',trim(xvvfile)
   if(len_trim(guvfile) > 0)&
        write(6,9701) 'Guv',trim(Guvfile)
   if(len_trim(huvfile) > 0)&
        write(6,9701) 'Huv',trim(Huvfile)
   if(len_trim(cuvfile) > 0)&
        write(6,9701) 'Cuv',trim(Cuvfile)
   if(len_trim(uuvfile) > 0)&
        write(6,9701) 'Uuv',trim(Uuvfile)
   if(len_trim(asympfile) > 0)&
        write(6,9701) 'Asymptotics',trim(asympfile)
   if(len_trim(quvfile) > 0)&
        write(6,9701) 'Quv',trim(Quvfile)
   if(len_trim(chgDistfile) > 0)&
        write(6,9701) 'ChgDist',trim(chgDistfile)
   if(len_trim(exchemfile) > 0)&
        write(6,9701) 'ExChem',trim(exchemfile)
   if(len_trim(solvenefile) > 0)&
        write(6,9701) 'SolvEne',trim(solvenefile)
   if(len_trim(entropyfile) > 0)&
        write(6,9701) 'Entropy',trim(entropyfile)
   if(len_trim(exchemGFfile) > 0)&
        write(6,9701) 'ExChGF',trim(exchemGFfile)
   if(len_trim(solveneGFfile) > 0)&
        write(6,9701) 'SolvEneGF',trim(solveneGFfile)
   if(len_trim(entropyGFfile) > 0)&
        write(6,9701) '-TS_GF',trim(entropyGFfile)
   if(len_trim(exchemUCfile) > 0)&
        write(6,9701) 'ExChUC',trim(exchemUCfile)
   if(len_trim(solveneUCfile) > 0)&
        write(6,9701) 'SolvEneUC',trim(solveneUCfile)
   if(len_trim(entropyUCfile) > 0)&
        write(6,9701) '-TS_UC',trim(entropyUCfile)
   if(len_trim(potUVfile) > 0)&
        write(6,9701) 'PotUV',trim(potUVfile)
#endif /*RISMSANDER*/

   ! Echo the input file to the user:
   call echoin(5,6)
   !     ----- READ DATA CHARACTERIZING THE MD-RUN -----
   read(5,'(a80)') title
   !       ----read input in namelist format, first setting up defaults
#endif /* ifndef API */

   ierr = 0
   dtemp = RETIRED_INPUT_OPTION
   dxm   = RETIRED_INPUT_OPTION
   heat  = RETIRED_INPUT_OPTION
   timlim = RETIRED_INPUT_OPTION
   irest = 0
   ibelly = 0
   ipol = RETIRED_INPUT_OPTION
   iesp = 0
   ntx = 1
   ntxo = NO_INPUT_VALUE
   ig = 71277
   tempi = ZERO
   ntb = NO_INPUT_VALUE
   ntt = 0
   nchain = 1
   temp0 = 300.0d0
! PLUMED
   plumed = 0
   plumedfile = 'plumed.dat'
! END PLUMED
#ifdef LES
   ! alternate temp for LES copies, if negative then use single bath
   ! single bath not the same as 2 baths with same target T
   temp0les = -ONE
#endif
   rdt = 0
   ipimd =0
   itimass = 0   ! Default = no TI w.r.t. mass.
   ineb  =0

   tautp = ONE
   ntp = 0
   barostat = 1
   mcbarint = 100
   pres0 = ONE
   comp = 44.6d0
   taup = ONE
   npscal = 1
   nscm = 1000
   nstlim = 1
   t = ZERO
   dt = 0.001d0
   ntc = 1
   tol = 0.00001
   ntf = 1
   nsnb = 25
   cut =  NO_INPUT_VALUE_FLOAT
   dielc = ONE
   ntpr = 50
   ntwr = 500
   ntwx = 0
   ntwv = 0
   ntwf = 0
   ntwe = 0
   ipb = 0
   inp = 2

#ifdef RISMSANDER
   irism = 0
#endif /*RISMSANDER*/

   ntave = 0
   ioutfm = 0
   ntr = 0
   ntrx = 1
   ivcap = 0
   natcap = 0
   fcap = 1.5d0
   cutcap = 0.0d0
   xcap = 0.0d0
   ycap = 0.0d0
   zcap = 0.0d0
   forth = 1.5d0
   xlorth = -1.0d0
   ylorth = -1.0d0
   zlorth = -1.0d0
   xorth = 47114711.0d0
   yorth = 47114711.0d0
   zorth = 47114711.0d0
   numexchg = 0
   repcrd   = 1
   lj1264 = 0

   profile_mpi = 0 !whether to write profile_mpi timing file - default = 0 (NO).

   ! number of waters to keep for hybrid model,
   ! numwatkeep: the number of closest
   ! waters to keep. close is defined as close to non-water.
   ! for simulations with ions, ions should be stripped too
   ! or at least ignored in the "closest" calculation. this
   ! is not currently done.

   ! if it stays at -1 then we keep all waters
   ! 0 would mean to strip them all

    numwatkeep=-1

   ! hybridgb: gb model to use with hybrid REMD.
   hybridgb=0
 
   ! carlos targeted MD, like ntr
   
   itgtmd=0
   tgtrmsd=0.
   tgtmdfrc=0.
   tgtfitmask=''
   tgtrmsmask=''

   pencut = 0.1d0
   taumet = 0.0001d0
   omega = 500.0d0
   ipnlty = 1
   scalm = 100.0d0
   iscale = 0
   noeskp = 1
   nmropt = 0
   jar = 0
   tausw = 0.1d0
   imin = 0
   isftrp = 0
   rwell = ONE
   maxcyc = 1
   ncyc = 10
   ntmin = 1
   dx0 = 0.01d0
   drms = 1.0d-4
   vlimit = 20.0d0
   mxsub = 1
   jfastw = 0
   watnam = '    '
   owtnm =  '    '
   hwtnm1 = '    '
   hwtnm2 = '    '
   ntwprt = 0
   igb = 0
   alpb = 0
   Arad = 15.0d0
   rgbmax = 25.d0
   saltcon = ZERO

   !  default offset depends on igb value, and users need to
   !  be able to modify it, so we need to set a dummy value. if it's still the
   !  dummy after we read the namelist, we set the default based on igb. if not,
   !  we leave it at what the user set.
   !  best solution would be to create a GB namelist.
   offset = -999999.d0 
   gbneckscale = -999999.d0 

   iyammp = 0
   imcdo = -1
   gbsa = 0
   vrand=1000
   surften = 0.005d0
   iwrap = 0
   nrespa = 1
   nrespai = 1
   irespa = 1
   gamma_ln = ZERO
   extdiel = 78.5d0
   intdiel = ONE
   gbgamma = ZERO
   gbbeta = ZERO
   gbalpha = ONE

   ! Parameters for isokinetic integrator (OIN)
   nkija  = 1     ! Number of Nose-Hoover chains per atom
   idistr = 0     ! Compute distribution functions (1=yes)

   !Hai Nguyen: set default parameters for igb = 8
   ! NOTE THAT NONE OF THESE ARE USED UNLESS IGB=8, SO USERS SHOULD NOT EVEN SET
   ! THEM
   gbalphaH = 0.788440d0
   gbbetaH = 0.798699d0
   gbgammaH = 0.437334d0
   gbalphaC = 0.733756d0
   gbbetaC = 0.506378d0
   gbgammaC = 0.205844d0
   gbalphaN = 0.503364d0
   gbbetaN = 0.316828d0
   gbgammaN = 0.192915d0
   gbalphaOS = 0.867814d0
   gbbetaOS = 0.876635d0
   gbgammaOS = 0.387882d0
   gbalphaP = 1.0d0    !P parameters are not optimized yet
   gbbetaP = 0.8d0     !P parameters are not optimized yet
   gbgammaP = 4.85d0   !P parameters are not optimized yet
   !scaling parameters below will only be used for igb=8. 
   ! the actual code does not use these variables, it uses X(l96)
   ! if igb=8, we will use these to set the X(l96) array.
   Sh = 1.425952d0
   Sc = 1.058554d0
   Sn = 0.733599d0
   So = 1.061039d0
   Ss = -0.703469d0
   Sp = 0.5d0          !P parameters are not optimized for protein
   ! update gbneck2nu pars
   ! using offset and gb_neckscale parameters from GB8-protein for gbneck2nu 
   ! Scaling factors
   ! name of variables are different from sander's igb8
   ! we use below names in pmemd
   ! gbneck2nu
   screen_hnu = 1.696538d0
   screen_cnu = 1.268902d0
   screen_nnu = 1.4259728d0
   screen_onu = 0.1840098d0
   screen_pnu = 1.5450597d0
   !alpha, beta, gamma for each atome element
   gb_alpha_hnu = 0.537050d0
   gb_beta_hnu = 0.362861d0
   gb_gamma_hnu = 0.116704d0
   gb_alpha_cnu = 0.331670d0
   gb_beta_cnu = 0.196842d0
   gb_gamma_cnu = 0.093422d0
   gb_alpha_nnu = 0.686311d0
   gb_beta_nnu = 0.463189d0
   gb_gamma_nnu = 0.138722d0
   gb_alpha_osnu = 0.606344d0
   gb_beta_osnu = 0.463006d0
   gb_gamma_osnu = 0.142262d0
   gb_alpha_pnu = 0.418365d0
   gb_beta_pnu = 0.290054d0
   gb_gamma_pnu = 0.1064245d0
   ! End gbneck2nu 
         
   iconstreff = 0
   cut_inner = EIGHT
   icfe = 0
   clambda = ZERO
   klambda = 1
   ievb = 0
   rbornstat = 0
   idecomp = 0
   ! added a flag to control output of BDC/SDC synonymous with MMPBSA.py's
   ! version of the same variable.
   dec_verbose = 3 
   lastrst = 1
   lastist = 1
   restraintmask=''
   restraint_wt = ZERO
   bellymask=''
   noshakemask=''
   iwrap_mask=''  ! GMS: mask to wrap around if iwrap == 2
   crgmask=''

   icnstph = 0
   solvph = SEVEN
   ntcnstph = 10
   ntrelax = 500 ! how long to let waters relax
   mccycles = 1  ! How many cycles of Monte Carlo steps to run
   skmin = 50 !used by neb calculation
   skmax = 100 !used by neb calculation
   vv = 0 !velocity verlet -- off if vv/=1
   vfac = 0 !velocity verlet scaling factor, 0 by default
   tmode = 1 !default tangent mode for NEB calculation

   ifqnt = NO_INPUT_VALUE

   ifcr = 0 ! no charge relocation
   cropt = 0 ! 1-4 EEL is calculated with the original charges
   crcut = 3.0
   crskin = 2.0 
   crin = ''
   crprintcharges = 0

   ips = 0    ! no isotropic periodic sum
   raips=-1.0d0   ! automatically determined
   mipsx=-1   ! number of grids in x direction, <0 for automatically determined
   mipsy=-1   ! number of grids in y direction, <0 for automatically determined
   mipsz=-1   ! number of grids in z direction, <0 for automatically determined
   mipso=4    ! default 4th order b-spline
   gridips=2   ! grid size. used to determine grid number if not defined
   dvbips=1.0d-8   ! Volume change tolerance. aips will be done when change more than dvbips

   iamd = 0 ! No accelerated MD used
   iamdlag = 0 !frequency of boosting in steps
   EthreshD = 0.d0
   alphaD = 0.d0
   EthreshP = 0.d0
   alphaP = 0.d0
   w_amd = 0 ! windowed amd
   EthreshD_w = 0.d0
   alphaD_w = 0.d0
   EthreshP_w = 0.d0
   alphaP_w = 0.d0
 
   scaledMD = 0 ! No scaled MD used
   scaledMD_lambda = 0.d0

   iemap=0     ! no emap constraint
   gammamap=1     ! default friction constant for map motion, 1/ps
   isgld = 0   ! no self-guiding
   isgsta=1    ! Begining index of SGLD range
   isgend=0    ! Ending index of SGLD range
   fixcom=-1    ! fix center of mass in SGLD simulation
   tsgavg=0.2d0    !  Local averaging time of SGLD simulation
   sgft=-1.0d3      !  Guiding factor of SGLD simulation
   sgff=-1.0d3      !  Guiding factor of SGLD simulation
   sgfd=-1.0d3      !  Guiding factor of SGLD simulation
   tempsg=0.0d0    !  Guiding temperature of SGLD simulation
   treflf=0.0d0    !  Reference low frequency temperature of SGLD simulation
   tsgavp=2.0d0    !  Convergency time of SGLD simulation

   !     Check to see if "cntrl" namelist has been defined.
   mdin_cntrl=.false.
   mdin_qmmm = .false.
   mdin_ewald=.false.
   mdin_pb=.false.
#ifdef APBS
   mdin_apbs = .false.
#endif /* APBS */
   mdin_lmod=.false.
   mdin_amoeba=.false.
   mdin_sebomd=.false.
   iamoeba = 0
#ifdef MPI /* SOFT CORE */
   scalpha=0.5
   scbeta=12.0
   sceeorder=2
   ifsc=0
   logdvdl=0
   dvdl_norest=0
   dynlmb=0.0
   ifmbar=0
   bar_intervall=100
   bar_l_min=0.1
   bar_l_max=0.9
   bar_l_incr=0.1
   tishake = 0
   emil_sc = 0
#endif
   ilrt = 0
   lrt_interval = 50
   lrtmask=''
#ifdef DSSP
   idssp = 0
#endif
   emil_do_calc = 0
   
!  Constant Surface Tension
   csurften = 0      !constant surface tension off (valid options are 0,1,2,3)
   gamma_ten = 0.0d0 !0.0 dyne/cm - default used in charmm. Ignored for csurften=0
   ninterface = 2   !Number of interfaces in the surface tension (Must be greater than 2)

#ifdef API
   igb = input_options%igb
   alpb = input_options%alpb
   gbsa = input_options%gbsa
   lj1264 = input_options%lj1264
   ipb = input_options%ipb
   inp = input_options%inp
   vdwmeth = input_options%vdwmeth
   ew_type = input_options%ew_type
   extdiel = input_options%extdiel
   intdiel = input_options%intdiel
   rgbmax = input_options%rgbmax
   saltcon = input_options%saltcon
   cut = input_options%cut
   dielc = input_options%dielc
   ifqnt = input_options%ifqnt
   jfastw = input_options%jfastw
   ntf = input_options%ntf
   ntc = input_options%ntc
#  ifdef LES
   rdt = input_options%rdt
#  endif
#else /* NOT the API */
   call nmlsrc('cntrl',5,ifind)
   if (ifind /= 0) mdin_cntrl=.true.

   call nmlsrc('ewald',5,ifind)
   if (ifind /= 0) mdin_ewald=.true.

   call nmlsrc('pb',5,ifind)
   if (ifind /= 0) mdin_pb=.true.
   
   call nmlsrc('qmmm', 5, ifind)
   if (ifind /= 0) mdin_qmmm = .true.

#  ifdef APBS
   call nmlsrc('apbs',5,ifind)
   if (ifind /= 0) mdin_apbs=.true.
#  endif /* APBS */

   call nmlsrc('lmod',5,ifind)
   if (ifind /= 0) mdin_lmod=.true.

   call nmlsrc('amoeba',5,ifind)
   if (ifind /= 0) mdin_amoeba=.true.

   call nmlsrc('sebomd',5,ifind)
   if (ifind /= 0) mdin_sebomd=.true.

   call nmlsrc('xray',5,ifind)
   xray_active = (ifind /= 0)

   rewind 5
   if ( mdin_cntrl ) then
      read(5,nml=cntrl,err=999)
   else
      write(6, '(1x,a,/)') 'Could not find cntrl namelist'
      FATAL_ERROR
   end if
#endif /* ifdef API */

   if ( igb == 10 .and. ipb == 0 ) ipb = 2
   if ( igb == 0  .and. ipb /= 0 ) igb = 10

   if (plumed.eq.1) then
     write(6, '(1x,a,/)') 'PLUMED is on'
     write(6, '(1x,a,a,/)') 'PLUMEDfile is ',plumedfile
   endif
   
   if (ifqnt == NO_INPUT_VALUE) then
      ifqnt = 0 ! default value
      if (mdin_qmmm) then
         write(6, '(1x,a,/)') &
            '| WARNING qmmm namelist found, but ifqnt was not set! QMMM NOT &
            &active.'
      end if
   end if

   ! Now that we've read the input file, set up the defaults for variables
   ! whose values depend on other input values (ntb, cut)
   if (ntb == NO_INPUT_VALUE) then
      if (ntp > 0) then
         ntb = 2
      else if (igb > 0) then
         ntb = 0
      else
         ntb = 1
      end if
   end if

   if (ntxo == NO_INPUT_VALUE) then
#ifdef MPI
      if (rem < 0) then
         ntxo = 2
      else
         ntxo = 1
      end if
#else
      ntxo = 1
#endif
   end if

   if (cut == NO_INPUT_VALUE_FLOAT) then
      if (igb == 0) then
         cut = EIGHT
      else
         cut = 9999.d0
      end if
   end if

#ifdef RISMSANDER
   !force igb=6 to get vacuum electrostatics.  This must be done ASAP to ensure SANDER's
   !electrostatics are initialized properly
   rismprm%irism=irism
   if(irism/=0) then
#   ifndef API
      write(6,'(a)') "|3D-RISM Forcing igb=6"
#   endif
      igb=6
   end if
#endif /*RISMSANDER*/   


   if (ifqnt>0) then
      qmmm_nml%ifqnt = .true.
      if (saltcon /= 0.0d0) then
         qm_gb%saltcon_on = .true.
      else
         qm_gb%saltcon_on = .false.
      end if
      if (alpb == 1) then
         qm_gb%alpb_on = .true.
      else
         qm_gb%alpb_on = .false.
      end if
      if (igb == 10 .or. ipb /= 0) then
         write(6, '(1x,a,/)') 'QMMM is not compatible with Poisson Boltzmann (igb=10 or ipb/=0).'
         FATAL_ERROR
      end if
   else
      qmmm_nml%ifqnt = .false.
   end if

   if ( mdin_lmod ) then
      rewind 5
      call read_lmod_namelist()
   end if

   !--------------------------------------------------------------------
   !     --- vars have been read ---
   !--------------------------------------------------------------------
   
#ifndef API
   write(6,9309)
#endif
   
   ! emit warnings for retired cntrl namelist variables

   if ( dtemp /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/,a)') 'Warning: dtemp has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
   end if
   if ( dxm /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/,a)') 'Warning: dxm has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
            ! '  The step length will be unlimited.'
   end if
   if ( heat /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/,a)') 'Warning: heat has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
   end if
   
   if ( timlim /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/,a)') 'Warning: timlim has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
   end if

! Constant surface tension valid options
   if (csurften > 0) then
      if (csurften < 0 .or. csurften > 3) then
         write(6,'(/2x,a)') &
         'Invalid csurften value. csurften must be between 0 and 3'
         FATAL_ERROR
      end if
      if (ntb /= 2) then
         write(6,'(/2x,a)') &
         'ntb invalid. ntb must be 2 for constant surface tension.'
         FATAL_ERROR
      end if
      if (ntp < 2) then
         write(6,'(/2x,a)') &
         'ntp invalid. ntp must be 2 or 3 for constant surface tension.'
         FATAL_ERROR
      end if
      if (ninterface < 2) then
         write(6,'(/2x,a)') &
         'ninterface must be greater than 2 for constant surface tension.'
         FATAL_ERROR
      end if
   
      if (iamoeba > 0 ) then
         write(6,'(/2x,a)') &
         'Constant Surface Tension is not compatible with Amoeba Runs.'
         FATAL_ERROR
      end if
   
      if (ipimd > 0 ) then
         write(6,'(/2x,a)') &
         'Constant Surface Tension is not compatible with PIMD Runs.'
         FATAL_ERROR
      end if

   end if

! MC Barostat valid options. Some of these may work, but disable them until they
! are fully tested.
   
   if (ntp > 0 .and. barostat == 2) then
      inerr = 0
      if (ievb /= 0) then
         write(6, '(/2x,a)') 'AMOEBA is not compatible with the MC Barostat'
         inerr = 1
      end if
      if (ipimd /= 0) then
         write(6, '(/2x,a)') 'PIMD is not compatible with the MC Barostat'
         inerr = 1
      end if
      if (icfe /= 0) then
         write(6, '(/2x,a)') 'TI is not compatible with the MC Barostat'
         inerr = 1
      end if
#ifdef LES
      write(6, '(/2x,a)') 'LES is not compatible with the MC Barostat'
      inerr = 1
#endif
      ! Any others? Hopefully most or all of the above can be made compatible.
      if (inerr == 1) then
         FATAL_ERROR
      end if
   end if

#ifndef API
   call printflags()
#endif

   !--------------------------------------------------------------------
   ! If user has requested ewald electrostatics, read some more input
   !--------------------------------------------------------------------
   
#ifdef API
   if( igb == 0 .and. ipb == 0 ) call load_ewald_info(ntp)
#else
   if( igb == 0 .and. ipb == 0 ) call load_ewald_info(inpcrd,ntp)
#endif

   !--------------------------------------------------------------------
   ! parameters for IPS and for SGLD:
   ! ips=1  3D IPS for electrostatic and Lennard-Jones potentials
   ! ips=2  3D IPS for electrostatic potential only
   ! ips=3  3D IPS for  Lennard-Jones potential only
   ! ips=4  3D IPS/DFFT for electrostatic and Lennard-Jones potentials
   ! ips=5  3D IPS/DFFT for electrostatic potential only
   ! ips=6  3D IPS/DFFT for Lennard-Jones potential only
   !--------------------------------------------------------------------
   
   teips=.false.
   tvips=.false.
   teaips=.false.
   tvaips=.false.
   if((ips-4)*(ips-6) == 0 )tvaips =.true.
   if ( (ips-4)*(ips-5) == 0 )teaips =.true.
   if( tvaips.OR.( (ips -1)*(ips-3) == 0 ))tvips =.true.
   if( teaips.OR.((ips -1)*(ips-2) == 0 ))teips =.true.
   if( teips ) then
      use_pme = 0
      eedmeth = 6
   end if
   if( tvips ) then
      vdwmeth = 2
      if(use_pme/=0.and.tvaips)then
        mipsx=nfft1   ! number of grids in x direction, <0 for automatically determined
        mipsy=nfft2   ! number of grids in y direction, <0 for automatically determined
        mipsz=nfft3   ! number of grids in z direction, <0 for automatically determined
        mipso=order    ! default 6th order b-spline
      endif
   end if
   temap=iemap>0
   ishake = 0
   if (ntc > 1) ishake = 1
   
   !--------------------------------------------------------------------
   ! Set up some parameters for AMD simulations:
   ! AMD initialization
   ! iamd=0 no boost is used, 1 boost on the total energy, 
   ! 2 boost on the dohedrals, 3 boost on dihedrals and total energy
   !--------------------------------------------------------------------
   if(iamd.gt.0)then
      if(iamd.eq.1)then !only total potential energy will be boosted
         EthreshD=0.d0
         alphaD=0.d0
      else if(iamd.eq.2)then !only dihedral energy will be boosted
         EthreshP=0.d0
         alphaP=0.d0
      endif
      if(w_amd.gt.0)then
         if(iamd.eq.1)then !only total potential energy will be boosted
            EthreshD_w=0.d0
            alphaD_w=0.d0
         else if(iamd.eq.2)then !only dihedral energy will be boosted
            EthreshP_w=0.d0
            alphaP_w=0.d0
         endif
#ifndef API
         write(6,'(a,i3)')'| Using Windowed Accelerated MD (wAMD) &
                          &LOWERING BARRIERS to enhance sampling w_amd =', w_amd
         write(6,'(a,2f22.12)')'| AMD boost to total energy: EthreshP,alphaP',&
                           EthreshP, alphaP
         write(6,'(a,2f22.12)')'| AMD boost to dihedrals: EthreshD,alphaD',&
                           EthreshD,alphaD
         write(6,'(a,2f22.12)')'| AMD extra parameters boost to total energy: &
                           &EthreshP_w,alphaP_w', EthreshP_w, alphaP_w
         write(6,'(a,2f22.12)')'| AMD extra parameters boost to dihedrals: &
                           &EthreshD_w,alphaD_w', EthreshD_w, alphaD_w
      else
         write(6,'(a,i3)')'| Using Accelerated MD (AMD) RASING VALLEYS to &
                           &enhance sampling iamd =',iamd
         write(6,'(a,2f22.12)')'| AMD boost to total energy: EthreshP,alphaP', &
                           EthreshP, alphaP
         write(6,'(a,2f22.12)')'| AMD boost to dihedrals: EthreshD,alphaD', &
                           EthreshD, alphaD
#endif
      endif
   endif
 

   !--------------------------------------------------------------------
   ! Set up some parameters for scaledMD simulations:
   ! scaledMD initialization
   ! scaledMD=0 no scaling is used, 1 scale the potential energy 
   !--------------------------------------------------------------------
#ifndef API
   if(scaledMD.gt.0)then
      write(6,'(a,i3)')'| Using Scaled MD to enhance sampling scaledMD =',&
                        scaledMD
      write(6,'(a,f22.12)')'| scaledMD scaling factor lambda: ',scaledMD_lambda
   endif
#endif
   
   !--------------------------------------------------------------------
   ! Set up some parameters for GB simulations:
   !--------------------------------------------------------------------
   !Hai Nguyen: update offset = 0.09d0 for igb /= 8
   !I add this step because I want to use different offset value as default value
   !for igb = 8
   
   if ( igb == 8 ) then
     if (offset == -999999.d0) then
        offset = 0.195141d0  !set to default for igb=8
     end if 
     if (gbneckscale == -999999.d0) then
        gbneckscale = 0.826836d0  
     end if 
   else 
      ! not igb=8, use old defaults
      if (offset == -999999.d0) then 
         offset = 0.09d0
      end if 
      if (gbneckscale == -999999.d0) then
         gbneckscale = 0.361825d0
      end if 
   endif
   
   if( igb == 2 .or. hybridgb == 2 ) then
      !       --- use our best guesses for Onufriev/Case GB  (GB^OBC I)
      
      gbgamma = 2.90912499999d0  ! (the "99999" to force roundoff on print)
      gbbeta = ZERO
      gbalpha = 0.8d0
   end if

   if( igb == 5 .or. hybridgb == 5 ) then
      
      !       --- use our second best guesses for Onufriev/Case GB (GB^OBC II)
      
      gbgamma = 4.850d0
      gbbeta = 0.8d0
      gbalpha = ONE
   end if
   
   if( igb == 7 ) then
      
      !       --- use parameters for Mongan et al. CFA GBNECK
      
      gbgamma = 2.50798245d0
      gbbeta = 1.90792938d0
      gbalpha = 1.09511284d0
   end if
   
   !--------------------------------------------------------------------
   ! If user has requested PB electrostatics, read some more input
   !--------------------------------------------------------------------

   if ( igb == 10 .or. ipb /= 0 ) then
#ifdef MPI
      write(6,'(a)') "PBSA currently doesn't work with MPI inside SANDER."
      FATAL_ERROR
#endif /*MPI*/
      call pb_read
   end if

#ifdef APBS
   if ( mdin_apbs ) then
      call apbs_read
   end if
#endif /* APBS */

#ifndef API
   call xray_read_mdin(mdin_lun=5)

   call sebomd_namelist_default
   if (mdin_sebomd) then
     rewind 5
     call read_sebomd_namelist
   endif

   if( iamoeba == 1 ) then
      if( mdin_amoeba ) then
         call AMOEBA_read_mdin(5)
      else
        write(6,*) ' iamoeba is set but the &amoeba namelist was not found'
        FATAL_ERROR
      end if
   end if
#endif /* API */

   ! -------------------------------------------------------------------
   ! If the user has requested NMR restraints, do a cursory read of the
   ! restraints file(s) now to determine the amount of memory necessary
   ! for these restraints:
   ! -------------------------------------------------------------------
   
   if (jar == 1 ) nmropt = 1
   intreq = 0
   irlreq = 0
   if (nmropt > 0) then
      mxgrp = 0
      itotst = 1
      
      ! Set ITOTST to 0 if IMIN equals 1 (i.e. if minimization, not dynamics)
      ! This will cause any "time-averaged" requests to be over-ridden.
      
      if (imin == 1) then 
         itotst = 0
      end if 
      !         CALL AMOPEN(31,NMR,'O','F','R')
      call restlx(5,itotst,mxgrp,dt,6,ierr)
      !         CLOSE(31)
   end if

   ! -------------------------------------------------------------------
   ! If EMIL was requested, make sure it was compiled in, and validate.
   ! -------------------------------------------------------------------
  if( emil_do_calc .gt. 0 ) then
#ifdef EMIL 
    if( ntc .ne. 1 ) then
       write (6, '(a,a)') error_hdr, 'emil_do_calc == 1,'
       write (6, '(a)') '      and ntc != 1.'
       write (6, '(a)') '      Current thinking is that SHAKE and '
       write (6, '(a)') '      EMIL do not mix well, consider setting ntc = 1, ntf = 1,'
       write (6, '(a)') '      and dt = 0.001.'
       FATAL_ERROR
    end if
#else
    write (6, '(a,a)') error_hdr, 'emil_do_calc = 1,'
    write (6, '(a)') '      but AMBER was compiled with EMIL switched out.'
    write (6, '(a)') '      Run $AMBERHOME/configure --help for more info.'
    FATAL_ERROR
#endif
  end if


   ! Set the definition of the water molecule. The default definition is in
   ! WATDEF(4).
   
   read(watdef(1),'(A4)') iwtnm
   read(watdef(2),'(A4)') iowtnm
   read(watdef(3),'(A4)') ihwtnm(1)
   read(watdef(4),'(A4)') ihwtnm(2)
   if (watnam /= '    ') read(watnam,'(A4)') iwtnm
   if (owtnm /= '    ') read(owtnm, '(A4)') iowtnm
   if (hwtnm1 /= '    ') read(hwtnm1,'(A4)') ihwtnm(1)
   if (hwtnm2 /= '    ') read(hwtnm2,'(A4)') ihwtnm(2)

#if !defined(DISABLE_NCSU) && defined(MPI)
   call ncsu_on_mdread1()
#endif

   return

#ifndef API
   999 continue   ! bad cntrl read
   write(6,*) 'error in reading namelist cntrl'
   FATAL_ERROR

   ! --- input file polar opts read err trapping:
   
   9308 format(/10x,55('-'),/10x, &
         'Amber 14 SANDER                              2014', &
         /10x,55('-')/)
   9309 format(/80('-')/'   1.  RESOURCE   USE: ',/80('-')/)
   9700 format(/,'File Assignments:',/,15('|',a6,': ',a,/))
#  ifdef RISMSANDER
   9701 format('|',a6,': ',a)
#  endif /* RISMSANDER */
#  ifdef MPI
   9702 format(7('|',a10,': ',a,/))
#  endif /* MPI */
#endif /* API */

#undef FATAL_ERROR
