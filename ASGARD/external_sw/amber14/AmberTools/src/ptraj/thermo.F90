
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
!+ [Enter a one-line description of subroutine thermo here]
subroutine thermo(natoms,nvecs,ilevel,c,amass,freq,vtemp,evibn, &
      cvibn,svibn,t,patm)
   implicit double precision(a-h,o-z)
   
   !     given the structure of a molecule and its normal mode vibrational
   !     frequencies this routine uses standard statistical mechanical
   !     formulas for an ideal gas (in the canonical ensemble, see,
   !     for example, d. a. mcquarrie, "statistical thermodynamics",
   !     harper & row, new york, 1973, chapters 5, 6, and 8) to compute
   !     the entropy, heat capacity, and internal energy.
   
   !     the si system of units is used internally.  conversion to units
   !     more familiar to most chemists is made for output.
   
   logical linear
   double precision jpcal
   
   !     amass:   atomic weights, in amu.
   !     pmom:    principal moments of inertia, in amu-bohr**2 and
   !              in ascending order.
   !     freq:    vibrational frequencies, in cm**-1 and in ascending
   !              order
   !     c    :   coordinates in Angstroms
   !     vtemp:   vibrational temperatures, in kelvin.
   !     evibn:   contribution to e from the vibration n.
   !     cvibn:   contribution to cv from the vibration n.
   !     svibn:   contribution to s from the vibration n.
   !     t:       temperature
   !     patm:    pressure, in atmospheres
   
   dimension amass(*),freq(*),c(*)
   dimension vtemp(*),evibn(*),cvibn(*),svibn(*)
   dimension pmom(10)
   
   data zero,pt2,half,one,onept5/0.0d0,0.2d0,0.5d0,1.0d0,1.5d0/
   data two,twopt5,four,eight,akilo/2.0d0,2.5d0,4.0d0,8.0d0,1000.d0/
   data thresh/900.d0/
   data iout /6/
   data pstd   /1.01325d+05/
   
   1000 format(/20x,19(1h*),/20x,'- Thermochemistry -',/20x,19(1h*),/ /)
   1010 format(1x,'molecular mass (principal isotopes) ',f11.5,' amu')
   1020 format(/1x,'temperature ',f9.3,' kelvin', &
         /1x,'pressure    ',f9.5,' atm')
   !1030 format(/1x,'Warning-- assumptions made about the electronic ',
   !    +           'partition function',
   !    +       /1x,'          are not valid for multiplets!'/)
   1040 format(/1x,'internal energy:   ',f10.3,' joule/mol',9x, &
         f10.3,' kcal/mol' &
         /1x,'entropy:           ',f10.3,' joule/k-mol',7x, &
         f10.3,' cal/k-mol' &
         /1x,'heat capacity cv:  ',f10.3,' joule/k-mol',7x, &
         f10.3,' cal/k-mol')
   1050 format(/1x,'rotational symmetry number ',f3.0)
   1060 format(/1x,'Warning-- assumption of classical behavior for ', &
         'rotation', &
         /1x,'          may cause significant error'/)
   1070 format(/1x,'rotational temperatures (kelvin) ',3f12.5)
   1080 format(/1x,'rotational temperature (kelvin) ',f12.5)
   1090 format(/1x,'zero point vibrational energy ', &
         f12.1,' (joules/mol) ', &
         /1x,30x,f12.5,' (kcal/mol)', &
         /1x,30x,f12.7,' (hartree/particle)')
   1100 format(/1x,'Warning-- ',i3,' vibrations have low frequencies', &
         ' and may represent hindered ', &
         /1x,'        internal rotations.  The contributions ', &
         'printed below assume that these ', &
         /1x,'        really are vibrations.')
   1110 format(/1x,'vibrational temperatures: ',5f9.2)
   1120 format(1x,9x,'(kelvin)',9x,5f9.2)
   1125 format(1x,26x,5f9.2)
   1130 format(/ /1x, 10x, 'freq.',9x,'E',9x, 9x,'Cv',8x, 9x,'S')
   1140 format(1x, 15x, 5x,'joules/mol',4x, 1x,'joules/mol-kelvin',1x, &
         2x,'joules/mol-kelvin')
   1150 format(1x, 'Total',10x, 3(4x,f11.3,4x))
   1160 format(1x, 'translational',2x, 3(4x,f11.3,4x))
   1170 format(1x, 'rotational',5x, 3(4x,f11.3,4x))
   1180 format(1x, 'vibrational',4x, 3(4x,f11.3,4x))
   1190 format(1x, i5,f10.3, 3(4x,f11.3,4x))
   1200 format(1x, 9x,'cm**-1', 6x,'kcal/mol',5x, 3x,'cal/mol-kelvin',2x, &
         2x,'cal/mol-kelvin', &
         /8('----------'))
   1210 format(1x,i3,' imaginary frequencies ignored')
   1220 format(/1x,'principal moments of inertia (nuclei only) in ', &
         'amu-A**2:', &
         /1x,5x,3f12.2)
   
   !     tokg:    kilograms per amu.
   !     boltz:   boltzman constant, in joules per kelvin.
   !     planck:  planck constant, in joule-seconds.
   !     avog:    avogadro constant, in mol**(-1).
   !     jpcal:   joules per calorie.
   !     tomet:   metres per Angstrom.
   !     hartre:  joules per hartree.
   
   tokg   = 1.660531d-27
   boltz  = 1.380622d-23
   planck = 6.626196d-34
   avog   = 6.022169d+23
   jpcal  = 4.18674d+00
   tomet  = 1.0d-10
   hartre = 4.35981d-18
   
   !     compute the gas constant, pi, pi**2, and e.
   !     compute the conversion factors cal per joule and kcal per joule.
   
   gas  = avog * boltz
   pi   = four * datan(one)
   pipi = pi * pi
   e    = dexp(one)
   tocal  = one / jpcal
   tokcal = tocal / akilo
   
   !     print the temperature and pressure.
   
   p = pstd * patm
   write(iout,1000)
   write(iout,1020) t,patm
   rt = gas * t
   
   !     compute and print the molecular mass in amu, then convert to
   !     kilograms.
   
   weight = zero
   do 20 iat=1,natoms
      weight = weight + amass(iat)
   20 continue
   write(iout,1010) weight
   weight = weight * tokg
   
   !     trap non-unit multiplicities.
   
   !     if (multip .ne. 1) write(iout,1030)
   
   !     compute contributions due to translation:
   !        etran-- internal energy
   !        ctran-- constant v heat capacity
   !        stran-- entropy
   
   dum1 = boltz * t
   dum2 = (two*pi) ** onept5
   arg  = dum1 ** onept5  / planck
   arg  = (arg/p) * (dum1/planck)
   arg  = arg * dum2 * (weight/planck)
   arg  = arg * sqrt(weight) * e**twopt5
   stran = gas * log(arg)
   etran = onept5 * rt
   ctran = onept5 * gas
   
   !     Compute contributions due to electronic motion:
   !        It is assumed that the first electronic excitation energy
   !        is much greater than kt and that the ground state has a
   !        degeneracy of one.  Under these conditions the electronic
   !        partition function can be considered to be unity.  The
   !        ground electronic state is taken to be the zero of
   !        electronic energy.
   
   !  40 continue
   
   !     for monatomics print and return.
   
   if (natoms <= 1) then
      s  = stran * tocal
      e  = etran * tokcal
      cv = ctran * tocal
      write(iout,1040) etran,e,stran,s,ctran,cv
      return
   end if
   
   !     compute contributions due to rotation.
   
   !     Compute the principal moments of inertia, get the rotational
   !     symmetry number, see if the molecule is linear, and compute
   !     the rotational temperatures.  Note the imbedded conversion
   !     of the moments to SI units.
   
   call mofi(natoms,c,amass,pmom)
   write(iout,1220) (pmom(i),i=1,3)
   linear = .false.
   call symnum(natoms,amass,sn,linear)
   write(iout,1050) sn
   con = planck / (boltz*eight*pipi)
   con = (con / tokg)  *  (planck / (tomet*tomet))
   if (linear) then
      rtemp = con / pmom(3)
      if (rtemp < pt2) write(iout,1060)
      write(iout,1080) rtemp
   else
      rtemp1 = con / pmom(1)
      rtemp2 = con / pmom(2)
      rtemp3 = con / pmom(3)
      if (rtemp1 < pt2) write(iout,1060)
      write(iout,1070) rtemp1,rtemp2,rtemp3
   end if
   
   !         erot-- rotational contribution to internal energy.
   !         crot-- rotational contribution to cv.
   !         srot-- rotational contribution to entropy.
   
   if (linear) then
      erot = rt
      crot = gas
      arg  = (t/rtemp) * (e/sn)
      srot = gas * log(arg)
   else
      erot = onept5 * rt
      crot = onept5 * gas
      arg  = sqrt(pi*e*e*e) / sn
      dum  = (t/rtemp1) * (t/rtemp2) * (t/rtemp3)
      arg  = arg * sqrt(dum)
      srot = gas * log(arg)
   end if
   
   !     compute contributions due to vibration.
   
   !     compute vibrational temperatures and zero point vibrational
   !     energy.  only real frequencies are included in the analysis.
   
   !     ndof = 3*natoms - 6 - nimag
   !     if (nimag .ne. 0) write(iout,1210) nimag
   !     if (linear) ndof = ndof + 1
   ndof = nvecs
   
   !       (---iff is the first frequency to include in thermo:)
   
   if (ilevel /= 0) then
      iff = 0
   else if (linear) then
      iff = 5
   else
      iff = 6
   end if
   con = planck / boltz
   ezpe = zero
   do 160 i=1,ndof
      vtemp(i) = freq(i+iff) * con * 3.0d10
      ezpe     = ezpe + freq(i+iff) * 3.0d10
   160 continue
   ezpe = half * planck * ezpe
   ezj  = ezpe * avog
   ezkc = ezpe * tokcal * avog
   ezau = ezpe / hartre
   write(iout,1090) ezj, ezkc, ezau
   
   !     compute the number of vibrations for which more than 5% of an
   !     assembly of molecules would exist in vibrational excited states.
   !     special printing for these modes is done to allow the user to
   !     easily take internal rotations into account.  the criterion
   !     corresponds roughly to a low frequency of 1.9(10**13) hz, or
   !     625 cm**(-1), or a vibrational temperature of 900 k.
   
   lofreq = 0
   do 180 i=1,ndof
      if (vtemp(i) < thresh) lofreq = lofreq + 1
   180 continue
   if (lofreq /= 0) write(iout,1100) lofreq
   
   !     compute:
   !        evib-- the vibrational component of the internal energy.
   !        cvib-- the vibrational component of the heat capacity.
   !        svib-- the vibrational component of the entropy.
   
   evib = zero
   cvib = zero
   svib = zero
   do 220 i=1,ndof
      
      !       compute some common factors.
      
      tovt  = vtemp(i) / t
      etovt = dexp(tovt)
      em1   = etovt - one
      
      !       compute contributions due to the i'th vibration.
      
      econt = tovt  *  (half + one/em1)
      ccont = etovt *  (tovt/em1)**2
      argd = one - one/etovt
      if (argd > 1.0d-7) then
         scont = tovt/em1 - log(argd)
      else
         scont = 0.0
         write (6,*) 'warning: setting vibrational entropy to zero ', &
               'for mode ',i,' with vtemp = ',vtemp(i)
      end if
      !       if (lofreq .ge. i) then
      evibn(i) = econt * rt
      cvibn(i) = ccont * gas
      svibn(i) = scont * gas
      !       end if
      evib = evib + econt
      cvib = cvib + ccont
      svib = svib + scont
   220 continue
   evib = evib * rt
   cvib = cvib * gas
   svib = svib * gas
   
   !     the units are now:
   !         e-- joules/mol
   !         c-- joules/mol-kelvin
   !         s-- joules/mol-kelvin
   
   etot = etran + erot + evib
   ctot = ctran + crot + cvib
   stot = stran + srot + svib
   
   !     print the sum of the hartree-fock energy and the thermal energy.
   
   !     call tread(501,gen,47,1,47,1,0)
   !     esum = gen(32) + etot/avog/hartre
   !     write(iout,1230) esum
   
   !     convert to the following and print
   !         e-- kcal/mol
   !         c-- cal/mol-kelvin
   !         s-- cal/mol-kelvin
   
   240 continue
   etran = etran * tokcal
   ctran = ctran * tocal
   stran = stran * tocal
   erot   = erot * tokcal
   crot   = crot * tocal
   srot   = srot * tocal
   evib   = evib * tokcal
   cvib   = cvib * tocal
   svib   = svib * tocal
   etot   = etran + erot + evib
   ctot   = ctran + crot + cvib
   stot   = stran + srot + svib
   do 280 i=1,ndof
      evibn(i) = evibn(i) * tokcal
      cvibn(i) = cvibn(i) * tocal
      svibn(i) = svibn(i) * tocal
   280 continue
   
   write(iout,1130)
   write(iout,1200)
   write(iout,1150) etot,ctot,stot
   write(iout,1160) etran,ctran,stran
   write(iout,1170) erot,crot,srot
   write(iout,1180) evib,cvib,svib
   do 318 i=1,iff
      write(iout,1190) i,freq(i)
   318 continue
   do 320 i=1,ndof
      write(iout,1190) i+iff,freq(i+iff),evibn(i),cvibn(i),svibn(i)
   320 continue
   
   return
end subroutine thermo 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine symnum here]
subroutine symnum(natoms,amass,sn,linear)
   implicit double precision(a-h,o-z)
   logical linear
   
   !     ----- routine to give the symmetry number. only for linear
   !           molecules. for others symmetry number is unity -----
   
   dimension amass(*)
   
   sn = 1.0d0
   if(natoms <= 2) then
      linear = .true.
      if(amass(1) == amass(2)) sn = 2.0d0
   end if
   return
end subroutine symnum 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mofi here]
subroutine mofi(natoms,c,amass,pmom)
   implicit double precision(a-h,o-z)
   
   !     compute the principal moments of inertia.
   !     units are amu-bohr**2
   
   character(len=1) uplo,jobz
   data uplo,jobz / 'U','V' /
   dimension c(*),amass(*),pmom(*)
   dimension com(3),t(9),e2(30),eigvec(9)
   
   data zero/0.0d0/
   
   !     stetement function to get ccom.
   
   ccom(ixyz,iat)=c(ixyz+3*(iat-1))-com(ixyz)
   
   
   !     compute the position of the center of mass and translate
   !     it to the origin.
   
   com(1) = zero
   com(2) = zero
   com(3) = zero
   
   totwt = zero
   do 20 iat=1,natoms
      iaind = 3*(iat-1)
      wt = amass(iat)
      totwt = totwt + wt
      com(1) = com(1) + wt*c(1+iaind)
      com(2) = com(2) + wt*c(2+iaind)
      com(3) = com(3) + wt*c(3+iaind)
   20 continue
   
   com(1) = com(1) / totwt
   com(2) = com(2) / totwt
   com(3) = com(3) / totwt
   
   !     compute the principal moments.
   
   do 60 i=1,9
      t(i) = zero
   60 continue
   
   do 80 iat=1,natoms
      wt = amass(iat)
      x  = ccom(1,iat)
      y  = ccom(2,iat)
      z  = ccom(3,iat)
      t(1) = t(1) + wt * (y*y+z*z)
      t(3) = t(3) + wt * (x*x+z*z)
      t(6) = t(6) + wt * (x*x+y*y)
      t(2) = t(2) - wt * x * y
      t(4) = t(4) - wt * x * z
      t(5) = t(5) - wt * y * z
   80 continue
   ier = 0
   call dspev(jobz,uplo,3,t,pmom,eigvec,3,e2,ier)
   return
end subroutine mofi 
