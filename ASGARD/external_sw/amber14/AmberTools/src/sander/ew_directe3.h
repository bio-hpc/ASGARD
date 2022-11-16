! epilogue: softcore (modified 12-6) LJ terms for V1, equals appearing atoms

do im_new = 1,icount
   j = cache_bckptr(im_new) ! atom# of atom j

   df =   cache_df(im_new)     ! electrostatic energy times 1/r^2 for the forces
   delx = cache_x(im_new)      ! delta x,y and z
   dely = cache_y(im_new)
   delz = cache_z(im_new)
   delr2 = cache_r2(im_new)    ! this contains r^2, note the difference from ew_directe2

   rfour=delr2*delr2
   r6=rfour*delr2

   ic = ico(iaci+iac(j))      ! locate the index in the vdW parameter array
   if ( ic < 0 ) then
      b0 = 1.0
      b1 = 0.0
   else
      b0 = sigma6(ic)
      b1 = foureps(ic)
   end if

   if ( nsc(i) == nsc(j) .and. emil_sc .eq. 0) then 
      !both atoms are softcore atoms and have normal 6-12 vdw
      denom = 1.0d0 / ( r6 * b0 )
      denom2 = denom * denom

      sc_ener(7) = sc_ener(7) + b1 * ( denom2 - denom ) ! Potential goes into the softcore energy array
      df = df + oneweight * b1 * ( 12.0d0 * denom2 - 6.0d0 * denom ) / delr2 ! scaled up by oneweight
   else 
      ! use the softcore potential
      denom = 1.0d0 / ( scalpha * ( 1.0d0 - clambda ) + r6 * b0 ) !b0 is 1/(sigma^6)
      denom2 = denom * denom
      denom3 = denom2 * denom

      evdw = evdw + b1 * ( denom2 - denom ) ! softcore potential is part of van der Waals energy

      ! -- ti decomp
      if(decpr .and. idecomp > 0) call decpair(3,i,j,b1*(denom2 - denom)/(nstlim/ntpr))

      sc_dvdl = sc_dvdl + b1 * ( 2.0d0 * scalpha * denom3 - scalpha * denom2 )
      ! -- ti decomp
      if(decpr .and. idecomp > 0) call decpair(3,i,j,weight1*b1*(2.0d0*scalpha*denom3 - scalpha*denom2)/(nstlim/ntpr))

      df = df + b1 * ( 12.0d0 * rfour * b0 * denom3 - 6.0d0 * rfour * b0 * denom2 )

      ! collect lambda-dependent contributions to BAR FEP energy differences
      if (ifmbar /= 0 .and. do_mbar) then
         do bar_i = 1, bar_states
            ! remove the current lambda cont
            bar_cont(bar_i) = bar_cont(bar_i) - b1 * ( denom2 - denom )
         end do
         do bar_i = 1, bar_states
            denom = 1.0d0 / ( scalpha * ( 1.0d0 - bar_lambda(bar_i) ) + r6 * b0 ) !b0 is 1/(sigma^6)
            denom2 = denom * denom
            bar_cont(bar_i) = bar_cont(bar_i) + b1 * ( denom2 - denom )
         end do
      end if

   end if

   dfx = delx*df
   dfy = dely*df
   dfz = delz*df

#ifndef noVIRIAL
   vxx = vxx - dfx*delx
   vxy = vxy - dfx*dely
   vxz = vxz - dfx*delz
   vyy = vyy - dfy*dely
   vyz = vyz - dfy*delz
   vzz = vzz - dfz*delz
#endif

   dumx = dumx + dfx
   dumy = dumy + dfy
   dumz = dumz + dfz
   force(1,j) = force(1,j) + dfx
   force(2,j) = force(2,j) + dfy
   force(3,j) = force(3,j) + dfz

end do  !  im_new = 1,icount
