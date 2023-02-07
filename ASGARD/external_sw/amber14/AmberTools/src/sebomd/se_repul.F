      subroutine se_repul(iai,norbsi,ai0,ai1,ai2,di1,di2,xi,yi,zi,
     .                 iaj,norbsj,aj0,aj1,aj2,dj1,dj2,xj,yj,zj,
     .                 rij,repij)

      implicit none
#include "sebomd.dim"
#include "se_constants.h"
!
!     routine to calculate repulsion integrals for a pair of atoms i,j.
!
!     input:
!
!     iai,iaj     = atomic numbers for atoms i and j.
!
!     norbsi,
!     norbsj,     = numbers of atomic orbitals centered on atoms i and j.
!
!     ai0,ai1,ai2,
!     aj0,aj1,aj2 = monopole, dipole, and quadrupole additive terms.
!
!     di1,di2,
!     dj1,dj2,    = dipole and quadrupole charge separations.
!
!     xi,yi,zi,
!     xj,yj,zj,   = cartesian coordinates (angstroms).
!
!     rij         = interatomic distance (angstroms).
!
!
!     returned:
!
!     repij = 10x10 array of repulsion integrals organized according to:
!
!             s_s  s_px  s_py  s_pz  px_px  px_py  px_pz  py_py  py_pz  pz_pz
!     s_s    ----------------------------------------------------------------
!     s_px   |                  .
!     s_py   |                  .
!     s_pz   |                  .
!     px_px  |                  .
!     px_py  |................(6,4) = (px(i)*py(i)|s(j)*pz(j))
!     px_pz  |
!     py_py  |
!     py_pz  |
!     pz_pz  |
!
      integer iai, norbsi
      integer iaj, norbsj
      double precision ai0, ai1, ai2, di1, di2, xi, yi, zi
      double precision aj0, aj1, aj2, dj1, dj2, xj, yj, zj
      double precision rij, repij(10,10)
!
!     local variables:
!
      double precision reploc(22),sqroot(100)
      double precision ee
      double precision eehalf
      double precision ee4th
      double precision ee8th
      double precision ee16th

      integer i, j
      double precision r
      double precision ra
      double precision rpsqr
      double precision rsqr
      double precision rp
      double precision rm
      double precision rma
      double precision rmm
      double precision rpa
      double precision rpp
      double precision rpm
      double precision rmp, rmj, rmi
      double precision rpisq, rmisq
      double precision rmsqr, rpjsq, rmjsq

      double precision a0
      double precision a00
      double precision a10
      double precision a01
      double precision a11
      double precision a20
      double precision a02
      double precision a22
      double precision a12
      double precision a21
      double precision aii2sq, ajj2sq

      double precision ri0, ri1, ri2, rpi
      double precision rj0, rj1, rj2, rpj

      double precision qq

      double precision di22
      double precision dj22
      double precision di22sq
      double precision dj22sq
      double precision dipj, dipjsq
      double precision dimj, dimjsq
      double precision dzqxx
      double precision dzqzz
      double precision di2a
      double precision dj2a
      double precision dpsqa
      double precision dmsqa

      double precision qdz
      double precision qxxq
      double precision qqxx
      double precision qqzz
      double precision qzzq
      double precision qxyqxy
      double precision qxzqxz
      double precision qxxqyy
      double precision qxxqzz
      double precision qzzqxx
      double precision qxxqxx
      double precision qzzqzz
      double precision qqqqxx
      double precision qqqqzz

      double precision dp
      double precision dm
      double precision qxxdz
      double precision qzzdz
      double precision qxzdx
      double precision dxqxz
      double precision dzq
      double precision dxdx
      double precision dzdz

!
!     ee = conversion factor from atomic units to ev.
!
      data ee     /27.210000d0/
      data eehalf /13.605000d0/
      data ee4th  / 6.802500d0/
      data ee8th  / 3.401250d0/
      data ee16th / 1.700625d0/
      save ee,eehalf,ee4th,ee8th,ee16th
!
      do 20 i=1,10
        do 10 j=1,10
          repij(i,j) = 0.0d0
 10     continue
 20   continue
!
!     convert rij to atomic units.
!
      a0 = bohr2ang
      r = rij/a0
      rsqr = r*r
!
!     assign semiempirical "additive" terms that determine one-center
!     limits for two-center integrals.
!
      ri0 = 0.5d0/ai0
      rj0 = 0.5d0/aj0
      if(norbsi.gt.1)then
        ri1 = 0.5d0/ai1
        ri2 = 0.5d0/ai2
      endif
      if(norbsj.gt.1)then
        rj1 = 0.5d0/aj1
        rj2 = 0.5d0/aj2
      endif
!
!     compute the nonzero multipole-multipole interactions for a local
!     frame with xi=xj, yi=yj, and zj-zi=r, i.e.,
!
!                           z
!                           |
!                           |
!                           |
!                   atom i  |_______ y
!                          /.
!                         / .
!                        x  .
!                           .
!                           .
!                           z
!                           |
!                           |
!                           |
!                   atom j  |_______ y
!                          /
!                         /
!                        x
!
!
!     [q(i),q(i)]:
!
      a00 = (ri0 + rj0)**2
      sqroot(1) = dsqrt(rsqr + a00)
      qq = ee/sqroot(1)
      if(norbsi.gt.1)then
!
!       [dz(i),q(j)]:
!
        a10 = (ri1 + rj0)**2
        sqroot(2) = dsqrt((r+di1)**2 + a10)
        sqroot(3) = dsqrt((r-di1)**2 + a10)
        dzq = eehalf/sqroot(2) - eehalf/sqroot(3)
!
!       [qzz(i),q(j)]:
!
        a20 = (ri2 + rj0)**2
        di22 = 2.0d0*di2
        sqroot(4) = dsqrt((r+di22)**2 + a20)
        sqroot(5) = dsqrt((r-di22)**2 + a20)
        sqroot(6) = dsqrt(rsqr + a20)
        qzzq = ee4th/sqroot(4) + ee4th/sqroot(5) - eehalf/sqroot(6)
!
!       [qxx(i),q(j)] or (qyy(i),q(j)]:
!
        di22sq = di22*di22
        sqroot(7) = dsqrt(rsqr + a20 + di22sq)
        sqroot(8) = sqroot(6)
        qxxq = eehalf/sqroot(7) - eehalf/sqroot(8)
!
        if(norbsj.gt.1)then
!
!         [qxz(i),dx(j)] or [qyz(i),dy(j)]:
!
          a21 = (ri2 + rj1)**2
          rp = (r + di2)**2
          rm = (r - di2)**2
          rpa = rp + a21
          rma = rm + a21
          dp = (di2 + dj1)**2
          dm = (di2 - dj1)**2
          sqroot(9) = dsqrt(rpa + dm)
          sqroot(10) = dsqrt(rma + dm)
          sqroot(11) = dsqrt(rma + dp)
          sqroot(12) = dsqrt(rpa + dp)
          qxzdx = ee4th/sqroot(9) - ee4th/sqroot(10)
     .          + ee4th/sqroot(11) - ee4th/sqroot(12)
!
!         [qzz(i),dz(j)]:
!
          rp = r + dj1
          rpsqr = rp*rp
          rm = r - dj1
          rmsqr = rm*rm
          rpp = (rp + di22)**2
          rmm = (rm - di22)**2
          rpm = (rp - di22)**2
          rmp = (rm + di22)**2
          sqroot(13) = dsqrt(rmm + a21)
          sqroot(14) = dsqrt(rpp + a21)
          sqroot(15) = dsqrt(rmp + a21)
          sqroot(16) = dsqrt(rpm + a21)
          sqroot(17) = dsqrt(rmsqr + a21)
          sqroot(18) = dsqrt(rpsqr + a21)
          qzzdz = ee8th/sqroot(13) - ee8th/sqroot(14) + ee8th/sqroot(15)
     .          - ee8th/sqroot(16) - ee4th/sqroot(17) + ee4th/sqroot(18)
!
!         [qxx(i),dz(j)] or [qyy(i),dz(j)]:
!
          di2a = di22sq + a21
          sqroot(19) = dsqrt(rmsqr + di2a)
          sqroot(20) = dsqrt(rpsqr + di2a)
          sqroot(21) = sqroot(17)
          sqroot(22) = sqroot(18)
          qxxdz = ee4th/sqroot(19) - ee4th/sqroot(20)
     .          - ee4th/sqroot(21) + ee4th/sqroot(22)
        endif
      endif
      if(norbsj.gt.1)then
!
!       [q(i),dz(j)]:
!
        if(iai.eq.iaj)then
          qdz = -dzq
        else
          a01 = (ri0 + rj1)**2
          sqroot(23) = dsqrt((r+dj1)**2 + a01)
          sqroot(24) = dsqrt((r-dj1)**2 + a01)
          qdz = -eehalf/sqroot(23) + eehalf/sqroot(24)
        endif
!
!       [q(i),qzz(j)]:
!
        dj22 = 2.0d0*dj2
        if(iai.eq.iaj)then
          qqzz = qzzq
        else
          a02 = (ri0 + rj2)**2
          sqroot(25) = dsqrt((r+dj22)**2 + a02)
          sqroot(26) = dsqrt((r-dj22)**2 + a02)
          sqroot(27) = dsqrt(rsqr + a02)
          qqzz = ee4th/sqroot(25) + ee4th/sqroot(26) - eehalf/sqroot(27)
        endif
!
!       [q(i),qxx(j)]:
!
        dj22sq = dj22*dj22
        if(iai.eq.iaj)then
          qqxx = qxxq
        else
          sqroot(28) = dsqrt(rsqr + a02 + dj22sq)
          sqroot(29) = sqroot(27)
          qqxx = eehalf/sqroot(28) - eehalf/sqroot(29)
        endif
        if(norbsi.gt.1)then
!
!         [dx(i),qxz(j)] or [dy(i),qyz(j)]:
!
          if(iai.eq.iaj)then
            dxqxz = -qxzdx
          else
            a12 = (ri1 + rj2)**2
            rp = (r + dj2)**2
            rm = (r - dj2)**2
            rpa = rp + a12
            rma = rm + a12
            dp = (dj2+di1)**2
            dm = (dj2-di1)**2
            sqroot(30) = dsqrt(rpa + dm)
            sqroot(31) = dsqrt(rma + dm)
            sqroot(32) = dsqrt(rma + dp)
            sqroot(33) = dsqrt(rpa + dp)
            dxqxz = -ee4th/sqroot(30) + ee4th/sqroot(31)
     .              -ee4th/sqroot(32) + ee4th/sqroot(33)
          endif
!
!         [dz(i),qzz(j)]:
!
          if(iai.eq.iaj)then
            dzqzz = -qzzdz
          else
            rp = r + di1
            rpsqr = rp*rp
            rm = r - di1
            rmsqr = rm*rm
            rpp = (rp + dj22)**2
            rmm = (rm - dj22)**2
            rpm = (rp - dj22)**2
            rmp = (rm + dj22)**2
            sqroot(34) = dsqrt(rmm + a12)
            sqroot(35) = dsqrt(rpp + a12)
            sqroot(36) = dsqrt(rmp + a12)
            sqroot(37) = dsqrt(rpm + a12)
            sqroot(38) = dsqrt(rmsqr + a12)
            sqroot(39) = dsqrt(rpsqr + a12)
            dzqzz = -ee8th/sqroot(34) + ee8th/sqroot(35)
     .              -ee8th/sqroot(36) + ee8th/sqroot(37)
     .              +ee4th/sqroot(38) - ee4th/sqroot(39)
          endif
!
!         [dz(i),qxx(j)] or [dz(i),qyy(j)]:
!
          if(iai.eq.iaj)then
            dzqxx = -qxxdz
          else
            dj2a = dj22sq + a12
            sqroot(40) = dsqrt(rmsqr + dj2a)
            sqroot(41) = dsqrt(rpsqr + dj2a)
            sqroot(42) = sqroot(38)
            sqroot(43) = sqroot(39)
            dzqxx = -ee4th/sqroot(40) + ee4th/sqroot(41)
     .              +ee4th/sqroot(42) - ee4th/sqroot(43)
          endif
!
!         [dx(i),dx(j)] or [dy(i),dy(j)]:
!
          a11 = (ri1 + rj1)**2
          ra = rsqr + a11
          dipj = di1 + dj1
          dimj = di1 - dj1
          sqroot(44) = dsqrt(ra + dimj*dimj)
          sqroot(45) = dsqrt(ra + dipj*dipj)
          dxdx = eehalf/sqroot(44) - eehalf/sqroot(45)
!
!         [dz(i),dz(j)]:
!
          sqroot(46) = dsqrt((r+dimj)**2 + a11)
          sqroot(47) = dsqrt((r+dipj)**2 + a11)
          sqroot(48) = dsqrt((r-dipj)**2 + a11)
          sqroot(49) = dsqrt((r-dimj)**2 + a11)
          dzdz = ee4th/sqroot(46) - ee4th/sqroot(47)
     .         - ee4th/sqroot(48) + ee4th/sqroot(49)


!
!         [qzz(i),qzz(j)]:
!
          a22 = (ri2 + rj2)**2
          ra = rsqr + a22
          rpi = r + di22
          rmi = r - di22
          rpj = r + dj22
          rmj = r - dj22
          sqroot(50) = dsqrt((rpi - dj22)**2 + a22)
          sqroot(51) = dsqrt((rpi + dj22)**2 + a22)
          sqroot(52) = dsqrt((rmi - dj22)**2 + a22)
          sqroot(53) = dsqrt((rmi + dj22)**2 + a22)
          sqroot(54) = dsqrt(rpi*rpi + a22)
          sqroot(55) = dsqrt(rmi*rmi + a22)
          sqroot(56) = dsqrt(rpj*rpj + a22)
          sqroot(57) = dsqrt(rmj*rmj + a22)
          sqroot(58) = dsqrt(ra)
          qzzqzz = ee16th/sqroot(50) + ee16th/sqroot(51)
     .           + ee16th/sqroot(52) + ee16th/sqroot(53)
     .           - ee8th/sqroot(54) - ee8th/sqroot(55)
     .           - ee8th/sqroot(56) - ee8th/sqroot(57)
     .           + ee4th/sqroot(58)
!
!         [qxx(i),qxx(j)] or [qyy(i),qyy(i)]:
!
          sqroot(59) = dsqrt(ra + (di22-dj22)**2)
          sqroot(60) = dsqrt(ra + (di22+dj22)**2)
          sqroot(61) = dsqrt(ra + di22sq)
          sqroot(62) = dsqrt(ra + dj22sq)
          sqroot(63) = sqroot(58)
          qxxqxx = ee8th/sqroot(59) + ee8th/sqroot(60)
     .           - ee4th/sqroot(61) - ee4th/sqroot(62)
     .           + ee4th/sqroot(63)
!
!         [qzz(i),qxx(j)] or [qzz(i),qyy(i)]:
!
          ajj2sq = a22 + dj22sq
          rpisq = rpi*rpi
          rmisq = rmi*rmi
          sqroot(64) = dsqrt(rmisq + ajj2sq)
          sqroot(65) = dsqrt(rpisq + ajj2sq)
          sqroot(66) = sqroot(55)
          sqroot(67) = sqroot(54)
          sqroot(68) = sqroot(62)
          sqroot(69) = sqroot(58)
          qzzqxx = ee8th/sqroot(64) + ee8th/sqroot(65)
     .           - ee8th/sqroot(66) - ee8th/sqroot(67)
     .           - ee4th/sqroot(68) + ee4th/sqroot(69)
!
!         [qxx(i),qzz(j)] or [qyy(i),qzz(j)]
!
          if(iai.eq.iaj)then
            qxxqzz = qzzqxx
          else
            aii2sq = a22 + di22sq
            rpjsq = rpj*rpj
            rmjsq = rmj*rmj
            sqroot(70) = dsqrt(rmjsq + aii2sq)
            sqroot(71) = dsqrt(rpjsq + aii2sq)
            sqroot(72) = sqroot(57)
            sqroot(73) = sqroot(56)
            sqroot(74) = sqroot(61)
            sqroot(75) = sqroot(58)
            qxxqzz = ee8th/sqroot(70) + ee8th/sqroot(71)
     .             - ee8th/sqroot(72) - ee8th/sqroot(73)
     .             - ee4th/sqroot(74) + ee4th/sqroot(75)
          endif
!
!         [qxx(i),qyy(j)] or [qyy(i),qxx(j)]:
!
          sqroot(76) = dsqrt(ra + di22sq + dj22sq)
          sqroot(77) = sqroot(61)
          sqroot(78) = sqroot(62)
          sqroot(79) = sqroot(58)
          qxxqyy = ee4th/sqroot(76) - ee4th/sqroot(77)
     .           - ee4th/sqroot(78) + ee4th/sqroot(79)
!
!         [qxz(i),qxz(j)] or [qyz(i),qyz(j)]:
!

          dipj = di2 + dj2
          dimj = di2 - dj2
          dipjsq = dipj*dipj
          dimjsq = dimj*dimj
          rpp = (r+dipj)**2
          rmp = (r-dipj)**2
          rpm = (r+dimj)**2
          rmm = (r-dimj)**2
          dpsqa = dipjsq + a22
          dmsqa = dimjsq + a22
          sqroot(80) = dsqrt(rpm + dmsqa)
          sqroot(81) = dsqrt(rpm + dpsqa)
          sqroot(82) = dsqrt(rpp + dmsqa)
          sqroot(83) = dsqrt(rpp + dpsqa)
          sqroot(84) = dsqrt(rmp + dmsqa)
          sqroot(85) = dsqrt(rmp + dpsqa)
          sqroot(86) = dsqrt(rmm + dmsqa)
          sqroot(87) = dsqrt(rmm + dpsqa)
          qxzqxz = ee8th/sqroot(80) - ee8th/sqroot(81)
     .         - ee8th/sqroot(82) + ee8th/sqroot(83)
     .         - ee8th/sqroot(84) + ee8th/sqroot(85)
     .         + ee8th/sqroot(86) - ee8th/sqroot(87)
!
!     [qxy(i),qxy(j)] or  [qxx-yy(i),qxx-yy(j)]:
!
          qxyqxy = 0.5d0*(qxxqxx - qxxqyy)

        endif
      endif
!
!     now use the local multipole-multipole interactions to construct
!     the 22 nonzero local repulsion integrals.
!
!     (s s | s s):
!
      reploc(1) = qq
      if(norbsi.gt.1)then
!
!       (s pz | s s):
!
        reploc(2) = dzq
!
!       (pz pz | s s):
!
        reploc(3) = qzzq + qq
!
!       (px px | s s) or (py py | s s):
!
        reploc(4) = qxxq + qq
      endif
      if(norbsj.gt.1)then
!
!       (s s | s pz):
!
        reploc(5) = qdz
        if(norbsi.gt.1)then
!
!         (s pz | s pz):
!
          reploc(6) = dzdz
!
!         (s px | s px) or (s py | s py):
!
          reploc(7) = dxdx
!
!         (pz pz | s pz):
!
          reploc(8) = qzzdz + qdz
!
!         (px px | s pz) or (py py | s pz):
!
          reploc(9) = qxxdz + qdz
!
!         (px pz | s px):
!
          reploc(10) = qxzdx
        endif
!
!       (s s | pz pz):
!
        reploc(11) = qqzz + qq
!
!       (s s | px px) or (s s | py py):
!
        reploc(12) = qqxx + qq
        if(norbsi.gt.1)then
!
!         (s pz | pz pz):
!
          reploc(13) = dzqzz + dzq
!
!         (s pz | px px) or (s pz | py py):
!
          reploc(14) = dzqxx + dzq
!
!         (s px | px pz) or (s py | py pz):
!
          reploc(15) = dxqxz
          qqqqzz = qq + qqzz
!
!         (pz pz | pz pz):
!
          reploc(16) = qqqqzz + qzzq + qzzqzz
!
!         (px px | pz pz) or (py py | pz pz):
!
          reploc(17) = qqqqzz + qxxq + qxxqzz
          qqqqxx = qq + qqxx
!
!         (pz pz | px px) or (pz pz | py py):
!
          reploc(18) = qqqqxx + qzzq + qzzqxx
!
!         (px px | px px) or (py py | py py):
!
          reploc(19) = qqqqxx + qxxq + qxxqxx
!
!         (px pz | px pz) or (py pz | py pz):
!
          reploc(20) = qxzqxz
!
!         (px px | py py) or (py py | px px):
!
          reploc(21) = qqqqxx + qxxq +  qxxqyy
!
!         (px py | px py):
!
          reploc(22) = qxyqxy
        endif
      endif



!
!     rotate local repulsion integrals to molecular frame.
!
      if(norbsi.eq.1.and.norbsj.eq.1)then
        repij(1,1) = reploc(1)
      else
        call se_rrep(norbsi,norbsj,xi,yi,zi,xj,yj,zj,rij,reploc,repij)
      endif
       return
      end
!--------------------------------------------------------------------------------
      subroutine se_rrep(norbsi,norbsj,xi,yi,zi,xj,yj,zj,rij,reploc,
     &   repij)
      implicit none
      integer norbsi
      integer norbsj
      double precision xi, yi, zi
      double precision xj, yj, zj
      double precision rij
      double precision reploc(22),repij(10,10)

      double precision z12, z12sqr, yinv

      double precision r11r11, r11r12, r11r13 ! 1-1
      double precision         r12r12, r12r13
      double precision                 r13r13

      double precision r11r21, r11r22, r11r23 ! 1-2
      double precision         r12r22, r12r23
      double precision                 r13r23

      double precision r11r31, r11r32, r11r33 ! 1-3
      double precision r12r31, r12r32, r12r33
      double precision                 r13r33

      double precision r21r11, r21r12, r21r13 ! 2-1
      double precision         r22r12, r22r13
      double precision                 r23r13

      double precision r21r21, r21r22, r21r23 ! 2-2
      double precision         r22r22, r22r23
      double precision                 r23r23

      double precision r21r31, r21r32, r21r33 ! 2-3
      double precision r22r31, r22r32, r22r33
      double precision                 r23r33

      double precision r31r11, r31r12, r31r13 ! 3-1
      double precision         r32r12, r32r13
      double precision                 r33r13

      double precision r31r21, r31r22, r31r23 ! 3-2
      double precision         r32r22, r32r23
      double precision                 r33r23

      double precision r31r31, r31r32, r31r33 ! 3-3
      double precision         r32r32, r32r33
      double precision                 r33r33

      double precision radd01
      double precision radd02
      double precision radd03
      double precision radd04
      double precision radd05
      double precision radd06
      double precision radd07
      double precision radd08
      double precision radd09
      double precision radd10
      double precision radd11
      double precision radd12
      double precision radd13
      double precision radd14
      double precision radd15
      double precision radd16
      double precision radd17
!
!     rotates the local frame repulsion integrals to the molecular frame.
!
      double precision r(3,3)
!
!     generate rotation matrix.
!
      r(3,1) = (xi-xj)/rij
      r(3,2) = (yi-yj)/rij
      r(3,3) = (zi-zj)/rij
      if(abs(r(3,3)).gt.0.99999999d0)then
        r(3,1) = 0.0d0
        r(3,2) = 0.0d0
        r(3,3) = sign(1.0d0,r(3,3))
        r(1,1) = 1.0d0
        r(1,2) = 0.0d0
        r(1,3) = 0.0d0
        r(2,1) = 0.0d0
        r(2,2) = r(3,3)
        r(2,3) = 0.0d0
      else
        z12sqr = r(3,1)**2 + r(3,2)**2
        z12 = dsqrt(z12sqr)
        r(1,1) = -r(3,2)/z12
        r(1,2) = r(3,1)/z12
        r(1,3) = 0.0d0
        if(abs(r(3,3)).lt.1.0d-8)then
          r(2,1) = 0.0d0
          r(2,2) = 0.0d0
          r(2,3) = 1.0d0
        else
          yinv = 1.0d0/(z12*dsqrt(1.0d0 + z12sqr/r(3,3)**2))
          r(2,1) = r(3,1)*yinv
          r(2,2) = r(3,2)*yinv
          r(2,3) = -(z12sqr/r(3,3))*yinv
        endif
      endif
      r11r11 = r(1,1)*r(1,1)
      r11r12 = r(1,1)*r(1,2)
      r11r13 = r(1,1)*r(1,3)
      r12r12 = r(1,2)*r(1,2)
      r12r13 = r(1,2)*r(1,3)
      r13r13 = r(1,3)*r(1,3)
      r21r21 = r(2,1)*r(2,1)
      r21r22 = r(2,1)*r(2,2)
      r21r23 = r(2,1)*r(2,3)
      r22r22 = r(2,2)*r(2,2)
      r22r23 = r(2,2)*r(2,3)
      r23r23 = r(2,3)*r(2,3)
      r31r31 = r(3,1)*r(3,1)
      r31r32 = r(3,1)*r(3,2)
      r31r33 = r(3,1)*r(3,3)
      r32r32 = r(3,2)*r(3,2)
      r32r33 = r(3,2)*r(3,3)
      r33r33 = r(3,3)*r(3,3)
      if(norbsi.gt.1.or.norbsj.gt.1)then
        radd01 = r11r11+r21r21
        radd02 = r11r12+r21r22
        radd03 = r11r13+r21r23
        radd04 = r12r12+r22r22
        radd05 = r12r13+r22r23
        radd06 = r13r13+r23r23
      endif
      if(norbsi.gt.1.and.norbsj.gt.1)then
        r11r21 = r(1,1)*r(2,1)*2.0d0
        r11r22 = r(1,1)*r(2,2)
        r11r23 = r(1,1)*r(2,3)
        r11r31 = r(1,1)*r(3,1)*2.0d0
        r11r32 = r(1,1)*r(3,2)
        r11r33 = r(1,1)*r(3,3)
        r12r22 = r(1,2)*r(2,2)*2.0d0
        r12r23 = r(1,2)*r(2,3)
        r12r31 = r(1,2)*r(3,1)
        r12r32 = r(1,2)*r(3,2)*2.0d0
        r12r33 = r(1,2)*r(3,3)
        r13r23 = r(1,3)*r(2,3)*2.0d0
        r13r33 = r(1,3)*r(3,3)*2.0d0
        r21r11 = r11r21
        r21r12 = r(2,1)*r(1,2)
        r21r13 = r(2,1)*r(1,3)
        r21r21 = r(2,1)*r(2,1)
        r21r31 = r(2,1)*r(3,1)*2.0d0
        r21r32 = r(2,1)*r(3,2)
        r21r33 = r(2,1)*r(3,3)
        r22r12 = r12r22
        r22r13 = r(2,2)*r(1,3)
        r22r31 = r(2,2)*r(3,1)
        r22r32 = r(2,2)*r(3,2)*2.0d0
        r22r33 = r(2,2)*r(3,3)
        r23r13 = r13r23
        r23r33 = r(2,3)*r(3,3)*2.0d0
        r31r11 = r11r31
        r31r12 = r(3,1)*r(1,2)
        r31r13 = r(3,1)*r(1,3)
        r31r21 = r21r31
        r31r22 = r22r31
        r31r23 = r(3,1)*r(2,3)
        r32r12 = r12r32
        r32r13 = r(3,2)*r(1,3)
        r32r22 = r22r32
        r32r23 = r(3,2)*r(2,3)
        r33r13 = r13r33
        r33r23 = r23r33
        radd07 = r11r32+r12r31
        radd08 = r21r32+r22r31
        radd09 = r11r33+r31r13
        radd10 = r21r33+r31r23
        radd11 = r11r32+r31r12
        radd12 = r21r32+r31r22
        radd13 = r12r33+r32r13
        radd14 = r22r33+r32r23
        radd15 = r11r23+r21r13
        radd16 = r12r23+r22r13
        radd17 = r11r22+r21r12
      endif
!
!     (s s | s s):
!
      repij(1,1) = reploc(1)
      if(norbsj.gt.1)then
!
!       (s s | s px):
!
        repij(1,2) = reploc(5)*r(3,1)
!
!       (s s | s py):
!
        repij(1,3) = reploc(5)*r(3,2)
!
!       (s s | s pz):
!
        repij(1,4) = reploc(5)*r(3,3)
!
!       (s s | px px):
!
        repij(1,5) = reploc(11)*r31r31
     .             + reploc(12)*radd01
!
!       (s s | px py):
!
        repij(1,6) = reploc(11)*r31r32
     .             + reploc(12)*radd02
!
!       (s s | px pz):
!
        repij(1,7) = reploc(11)*r31r33
     .             + reploc(12)*radd03
!
!       (s s | py py):
!
        repij(1,8) = reploc(11)*r32r32
     .             + reploc(12)*radd04
!
!       (s s | py pz):
!
        repij(1,9) = reploc(11)*r32r33
     .             + reploc(12)*radd05
!
!       (s s | pz pz):
!
        repij(1,10) = reploc(11)*r33r33
     .              + reploc(12)*radd06
      endif
      if(norbsi.gt.1)then
!
!       (s px | s s):
!
        repij(2,1) = reploc(2)*r(3,1)
!
!       (s py | s s):
!
        repij(3,1) = reploc(2)*r(3,2)
!
!       (s pz | s s):
!
        repij(4,1) = reploc(2)*r(3,3)
!
!       (px px | s s):
!
        repij(5,1) = reploc(3)*r31r31
     .             + reploc(4)*radd01
!
!       (px py | s s):
!
        repij(6,1) = reploc(3)*r31r32
     .             + reploc(4)*radd02
!
!       (px pz | s s):
!
        repij(7,1) = reploc(3)*r31r33
     .             + reploc(4)*radd03
!
!       (py py | s s):
!
        repij(8,1) = reploc(3)*r32r32
     .             + reploc(4)*radd04
!
!       (py pz | s s);
!
        repij(9,1) = reploc(3)*r32r33
     .             + reploc(4)*radd05
!
!       (pz pz | s s):
!
        repij(10,1) = reploc(3)*r33r33
     .              + reploc(4)*radd06
      endif
      if(norbsi.gt.1.and.norbsj.gt.1)then
!
!       (s px | s px):
!
        repij(2,2) = reploc(6)*r31r31
     .             + reploc(7)*radd01
!
!       (s px | s py):
!
        repij(2,3) = reploc(6)*r31r32
     .             + reploc(7)*radd02
!
!       (s py | s px):
!
        repij(3,2) = repij(2,3)
!
!       (s px | s pz):
!
        repij(2,4) = reploc(6)*r31r33
     .             + reploc(7)*radd03
!
!       (s pz | s px):
!
        repij(4,2) = repij(2,4)
!
!       (s px | px px):
!
        repij(2,5) = reploc(13)*r(3,1)*r31r31
     .             + reploc(14)*r(3,1)*radd01
     .             + reploc(15)*(r(1,1)*r11r31
     .                         + r(2,1)*r21r31)
!
!       (px px | s px):
!
        repij(5,2) = reploc(8)*r(3,1)*r31r31
     .             + reploc(9)*r(3,1)*radd01
     .             + reploc(10)*(r(1,1)*r11r31
     .                         + r(2,1)*r21r31)
!
!       (s px | px py):
!
        repij(2,6) = reploc(13)*r(3,1)*r31r32
     .             + reploc(14)*r(3,1)*radd02
     .             + reploc(15)*(r(1,1)*radd07
     .                         + r(2,1)*radd08)
!
!       (px py | s px):
!
        repij(6,2) = reploc(8)*r(3,1)*r31r32
     .             + reploc(9)*r(3,1)*radd02
     .             + reploc(10)*(r(1,1)*radd07
     .                         + r(2,1)*radd08)
!
!       (s px | px pz):
!
        repij(2,7) = reploc(13)*r(3,1)*r31r33
     .             + reploc(14)*r(3,1)*radd03
     .             + reploc(15)*(r(1,1)*radd09
     .                         + r(2,1)*radd10)
!
!       (px pz | s px):
!
        repij(7,2) = reploc(8)*r(3,1)*r31r33
     .             + reploc(9)*(r(3,1)*r11r13
     .                        + r(3,1)*r21r23)
     .             + reploc(10)*(r(1,1)*radd09
     .                         + r(2,1)*radd10)
!
!       (s px | py py):
!
        repij(2,8) = reploc(13)*r(3,1)*r32r32
     .             + reploc(14)*r(3,1)*radd04
     .             + reploc(15)*(r(1,1)*r12r32
     .                         + r(2,1)*r22r32)
!
!       (py py | s px):
!
        repij(8,2) = reploc(8)*r(3,1)*r32r32
     .             + reploc(9)*r(3,1)*radd04
     .             + reploc(10)*(r(1,1)*r12r32
     .                         + r(2,1)*r22r32)
!
!       (s px | py pz):
!
!
        repij(2,9) = reploc(13)*r(3,1)*r32r33
     .             + reploc(14)*r(3,1)*radd05
     .             + reploc(15)*(r(1,1)*radd13
     .                         + r(2,1)*radd14)
!
!       (py pz | s px):
!
        repij(9,2) = reploc(8)*r(3,1)*r32r33
     .             + reploc(9)*r(3,1)*radd05
     .             + reploc(10)*(r(1,1)*radd13
     .                         + r(2,1)*radd14)
!
!       (s px | pz pz):
!
        repij(2,10) = reploc(13)*r(3,1)*r33r33
     .              + reploc(14)*r(3,1)*radd06
     .              + reploc(15)*(r(1,1)*r13r33
     .                          + r(2,1)*r23r33)
!
!       (pz pz | s px):
!
        repij(10,2) = reploc(8)*r(3,1)*r33r33
     .              + reploc(9)*r(3,1)*radd06
     .              + reploc(10)*(r(1,1)*r13r33
     .                          + r(2,1)*r23r33)
!
!       (s py | s py):
!
        repij(3,3) = reploc(6)*r32r32
     .             + reploc(7)*radd04
!
!       (s py | s pz):
!
        repij(3,4) = reploc(6)*r32r33
     .             + reploc(7)*radd05
!
!       (s pz | s py):
!
        repij(4,3) = repij(3,4)
!
!       (s py | px px):
!
        repij(3,5) = reploc(13)*r(3,2)*r31r31
     .             + reploc(14)*r(3,2)*radd01
     .             + reploc(15)*(r(1,2)*r11r31
     .                         + r(2,2)*r21r31)
!
!       (px px | s py):
!
        repij(5,3) = reploc(8)*r(3,2)*r31r31
     .             + reploc(9)*r(3,2)*radd01
     .             + reploc(10)*(r(1,2)*r11r31
     .                         + r(2,2)*r21r31)
!
!       (s py | px py):
!
        repij(3,6) = reploc(13)*r(3,2)*r31r32
     .             + reploc(14)*r(3,2)*radd02
     .             + reploc(15)*(r(1,2)*radd11
     .                         + r(2,2)*radd12)
!
!       (px py | s py):
!
        repij(6,3) = reploc(8)*r(3,2)*r31r32
     .             + reploc(9)*r(3,2)*radd02
     .             + reploc(10)*(r(1,2)*radd11
     .                         + r(2,2)*radd12)
!
!       (s py | px pz):
!
        repij(3,7) = reploc(13)*r(3,2)*r31r33
     .             + reploc(14)*r(3,2)*radd03
     .             + reploc(15)*(r(1,2)*radd09
     .                         + r(2,2)*radd10)
!
!       (px pz | s py):
!
        repij(7,3) = reploc(8)*r(3,2)*r31r33
     .             + reploc(9)*r(3,2)*radd03
     .             + reploc(10)*(r(1,2)*radd09
     .                         + r(2,2)*radd10)
!
!       (s py | py py):
!
        repij(3,8) = reploc(13)*r(3,2)*r32r32
     .             + reploc(14)*r(3,2)*radd04
     .             + reploc(15)*(r(1,2)*r12r32
     .                         + r(2,2)*r22r32)
!
!       (py py | s py):
!
        repij(8,3) = reploc(8)*r(3,2)*r32r32
     .             + reploc(9)*r(3,2)*radd04
     .             + reploc(10)*(r(1,2)*r12r32
     .                         + r(2,2)*r22r32)
!
!       (s py | py pz):
!
        repij(3,9) = reploc(13)*r(3,2)*r32r33
     .             + reploc(14)*r(3,2)*radd05
     .             + reploc(15)*(r(1,2)*radd13
     .                         + r(2,2)*radd14)
!
!       (py pz | s py):
!
        repij(9,3) = reploc(8)*r(3,2)*r32r33
     .             + reploc(9)*r(3,2)*radd05
     .             + reploc(10)*(r(1,2)*radd13
     .                         + r(2,2)*radd14)
!
!       (s py | pz pz):
!
        repij(3,10) = reploc(13)*r(3,2)*r33r33
     .              + reploc(14)*r(3,2)*radd06
     .              + reploc(15)*(r(1,2)*r13r33
     .                          + r(2,2)*r23r33)
!
!       (pz pz | s py):
!
        repij(10,3) = reploc(8)*r(3,2)*r33r33
     .              + reploc(9)*r(3,2)*radd06
     .              + reploc(10)*(r(1,2)*r13r33
     .                          + r(2,2)*r23r33)
!
!       (s pz | s pz):
!
        repij(4,4) = reploc(6)*r33r33
     .             + reploc(7)*radd06
!
!       (s pz | px px):
!
        repij(4,5) = reploc(13)*r(3,3)*r31r31
     .             + reploc(14)*r(3,3)*radd01
     .             + reploc(15)*(r(1,3)*r11r31
     .                         + r(2,3)*r21r31)
!
!       (px px | s pz):
!
        repij(5,4) = reploc(8)*r(3,3)*r31r31
     .             + reploc(9)*r(3,3)*radd01
     .             + reploc(10)*(r(1,3)*r11r31
     .                         + r(2,3)*r21r31)
!
!       (s pz | px py):
!
        repij(4,6) = reploc(13)*r(3,3)*r31r32
     .             + reploc(14)*r(3,3)*radd02
     .             + reploc(15)*(r(1,3)*radd11
     .                         + r(2,3)*radd12)
!
!       (px py | s pz):
!
        repij(6,4) = reploc(8)*r(3,3)*r31r32
     .             + reploc(9)*r(3,3)*radd02
     .             + reploc(10)*(r(1,3)*radd11
     .                         + r(2,3)*radd12)
!
!       (s pz | px pz):
!
        repij(4,7) = reploc(13)*r(3,3)*r31r33
     .             + reploc(14)*r(3,3)*radd03
     .             + reploc(15)*(r(1,3)*radd09
     .                         + r(2,3)*radd10)
!
!       (px pz | s pz):
!
        repij(7,4) = reploc(8)*r(3,3)*r31r33
     .             + reploc(9)*r(3,3)*radd03
     .             + reploc(10)*(r(1,3)*radd09
     .                         + r(2,3)*radd10)
!
!       (s pz | py py):
!
        repij(4,8) = reploc(13)*r(3,3)*r32r32
     .             + reploc(14)*r(3,3)*radd04
     .             + reploc(15)*(r(1,3)*r12r32
     .                         + r(2,3)*r22r32)
!
!       (py py | s pz):
!
        repij(8,4) = reploc(8)*r(3,3)*r32r32
     .             + reploc(9)*r(3,3)*radd04
     .             + reploc(10)*(r(1,3)*r12r32
     .                         + r(2,3)*r22r32)
!
!       (s pz | py pz):
!
        repij(4,9) = reploc(13)*r(3,3)*r32r33
     .             + reploc(14)*r(3,3)*radd05
     .             + reploc(15)*(r(1,3)*radd13
     .                         + r(2,3)*radd14)
!
!       (py pz | s pz):
!
        repij(9,4) = reploc(8)*r(3,3)*r32r33
     .             + reploc(9)*r(3,3)*radd05
     .             + reploc(10)*(r(1,3)*radd13
     .                         + r(2,3)*radd14)
!
!       (s pz | pz pz):
!
        repij(4,10) = reploc(13)*r(3,3)*r33r33
     .              + reploc(14)*r(3,3)*radd06
     .              + reploc(15)*(r(1,3)*r13r33
     .                          + r(2,3)*r23r33)
!
!       (pz pz | s pz):
!
        repij(10,4) = reploc(8)*r(3,3)*r33r33
     .              + reploc(9)*r(3,3)*radd06
     .              + reploc(10)*(r(1,3)*r13r33
     .                          + r(2,3)*r23r33)
!
!       (px px | px px):
!
        repij(5,5) = reploc(16)*r31r31*r31r31
     .             + reploc(17)*radd01*r31r31
     .             + reploc(18)*radd01*r31r31
     .             + reploc(19)*(r11r11*r11r11
     .                         + r21r21*r21r21)
     .             + reploc(20)*(r11r31*r11r31
     .                         + r21r31*r21r31)
     .             + reploc(21)*r11r11*r21r21*2.0d0
     .             + reploc(22)*r11r21*r11r21
!
!       (px px | px py):
!
        repij(5,6) = reploc(16)*r31r31*r31r32
     .             + reploc(17)*radd01*r31r32
     .             + reploc(18)*radd02*r31r31
     .             + reploc(19)*(r11r11*r11r12
     .                         + r21r21*r21r22)
     .             + reploc(20)*(r11r31*r11r32
     .                         + r31r11*r31r12
     .                         + r21r31*r21r32
     .                         + r31r21*r31r22)
     .             + reploc(21)*(r11r11*r21r22
     .                         + r21r21*r11r12)
     .             + reploc(22)*(r11r21*r11r22
     .                         + r21r11*r21r12)
!
!       (px py | px px):
!
        repij(6,5) = reploc(16)*r31r32*r31r31
     .             + reploc(17)*radd02*r31r31
     .             + reploc(18)*radd01*r31r32
     .             + reploc(19)*(r11r12*r11r11
     .                         + r21r22*r21r21)
     .             + reploc(20)*(radd11*r11r31
     .                         + radd12*r21r31)
     .             + reploc(21)*(r11r12*r21r21
     .                         + r21r22*r11r11)
     .             + reploc(22)*radd17*r11r21
!
!       (px px | px pz):
!
        repij(5,7) = reploc(16)*r31r31*r31r33
     .             + reploc(17)*radd01*r31r33
     .             + reploc(18)*radd03*r31r31
     .             + reploc(19)*(r11r11*r11r13
     .                         + r21r21*r21r23)
     .             + reploc(20)*(r11r31*r11r33
     .                         + r31r11*r31r13
     .                         + r21r31*r21r33
     .                         + r31r21*r31r23)
     .             + reploc(21)*(r11r11*r21r23
     .                         + r21r21*r11r13)
     .             + reploc(22)*(r11r21*r11r23
     .                         + r21r11*r21r13)
!
!       (px pz | px px):
!
        repij(7,5) = reploc(16)*r31r33*r31r31
     .             + reploc(17)*radd03*r31r31
     .             + reploc(18)*radd01*r31r33
     .             + reploc(19)*(r11r13*r11r11
     .                         + r21r23*r21r21)
     .             + reploc(20)*(radd09*r11r31
     .                          +radd10*r21r31)
     .             + reploc(21)*(r11r13*r21r21
     .                         + r21r23*r11r11)
     .             + reploc(22)*radd15*r11r21
!
!       (px px | py py):
!
        repij(5,8) = reploc(16)*r31r31*r32r32
     .             + reploc(17)*radd01*r32r32
     .             + reploc(18)*radd04*r31r31
     .             + reploc(19)*(r11r11*r12r12
     .                         + r21r21*r22r22)
     .             + reploc(20)*(r11r31*r12r32
     .                         + r21r31*r22r32)
     .             + reploc(21)*(r11r11*r22r22
     .                         + r21r21*r12r12)
     .             + reploc(22)*r11r21*r12r22
!
!     (py py | px px):
!
        repij(8,5) = reploc(16)*r32r32*r31r31
     .             + reploc(17)*radd04*r31r31
     .             + reploc(18)*radd01*r32r32
     .             + reploc(19)*(r12r12*r11r11
     .                         + r22r22*r21r21)
     .             + reploc(20)*(r12r32*r11r31
     .                         + r22r32*r21r31)
     .             + reploc(21)*(r12r12*r21r21
     .                         + r22r22*r11r11)
     .             + reploc(22)*r12r22*r11r21
!
!       (px px | py pz):
!
        repij(5,9) = reploc(16)*r31r31*r32r33
     .             + reploc(17)*radd01*r32r33
     .             + reploc(18)*radd05*r31r31
     .             + reploc(19)*(r11r11*r12r13
     .                         + r21r21*r22r23)
     .             + reploc(20)*(r11r31*r12r33
     .                         + r31r11*r32r13
     .                         + r21r31*r22r33
     .                         + r31r21*r32r23)
     .             + reploc(21)*(r11r11*r22r23
     .                         + r21r21*r12r13)
     .             + reploc(22)*(r11r21*r12r23
     .                         + r21r11*r22r13)
!
!       (py pz | px px):
!
        repij(9,5) = reploc(16)*r32r33*r31r31
     .             + reploc(17)*radd05*r31r31
     .             + reploc(18)*radd01*r32r33
     .             + reploc(19)*(r12r13*r11r11
     .                         + r22r23*r21r21)
     .             + reploc(20)*(radd13*r11r31
     .                         + radd14*r21r31)
     .             + reploc(21)*(r12r13*r21r21
     .                         + r22r23*r11r11)
     .             + reploc(22)*radd16*r11r21
!
!       (px px | pz pz):
!
        repij(5,10) = reploc(16)*r31r31*r33r33
     .              + reploc(17)*radd01*r33r33
     .              + reploc(18)*radd06*r31r31
     .              + reploc(19)*(r11r11*r13r13
     .                          + r21r21*r23r23)
     .              + reploc(20)*(r11r31*r13r33
     .                          + r21r31*r23r33)
     .              + reploc(21)*(r11r11*r23r23
     .                          + r21r21*r13r13)
     .              + reploc(22)*r11r21*r13r23
!
!       (pz pz | px px):
!
        repij(10,5) = reploc(16)*r33r33*r31r31
     .              + reploc(17)*radd06*r31r31
     .              + reploc(18)*radd01*r33r33
     .              + reploc(19)*(r13r13*r11r11
     .                          + r23r23*r21r21)
     .              + reploc(20)*(r13r33*r11r31
     .                          + r23r33*r21r31)
     .              + reploc(21)*(r13r13*r21r21
     .                          + r23r23*r11r11)
     .              + reploc(22)*r13r23*r11r21
!
!       (px py | px py):
!
        repij(6,6) = reploc(16)*r31r32*r31r32
     .             + reploc(17)*radd02*r31r32
     .             + reploc(18)*radd02*r31r32
     .             + reploc(19)*(r11r12*r11r12
     .                         + r21r22*r21r22)
     .             + reploc(20)*(radd11**2
     .                         + radd12**2)
     .             + reploc(21)*2.0d0*r11r12*r21r22
     .             + reploc(22)*radd17**2
!
!       (px py | px pz):
!
        repij(6,7) = reploc(16)*r31r32*r31r33
     .             + reploc(17)*radd02*r31r33
     .             + reploc(18)*r31r32*radd03
     .             + reploc(19)*(r11r12*r11r13
     .                         + r21r22*r21r23)
     .             + reploc(20)*(radd11*radd09
     .                         + radd12*radd10)
     .             + reploc(21)*(r11r12*r21r23
     .                         + r21r22*r11r13)
     .             + reploc(22)*radd17*(r11r23+r21r13)
!
!       (px pz | px py):
!
        repij(7,6) = repij(6,7)
!
!       (px py | py py):
!
        repij(6,8) = reploc(16)*r31r32*r32r32
     .             + reploc(17)*radd02*r32r32
     .             + reploc(18)*radd04*r31r32
     .             + reploc(19)*(r11r12*r12r12
     .                         + r21r22*r22r22)
     .             + reploc(20)*(radd11*r12r32
     .                         + radd12*r22r32)
     .             + reploc(21)*(r11r12*r22r22
     .                         + r21r22*r12r12)
     .             + reploc(22)*radd17*r12r22
!
!       (py py | px py):
!
        repij(8,6) = reploc(16)*r32r32*r31r32
     .             + reploc(17)*radd04*r31r32
     .             + reploc(18)*radd02*r32r32
     .             + reploc(19)*(r12r12*r11r12
     .                         + r22r22*r21r22)
     .             + reploc(20)*(r12r32*r11r32
     .                         + r32r12*r31r12
     .                         + r22r32*r21r32
     .                         + r32r22*r31r22)
     .             + reploc(21)*(r12r12*r21r22
     .                         + r22r22*r11r12)
     .             + reploc(22)*(r12r22*r11r22
     .                         + r22r12*r21r12)
!
!       (px py | py pz):
!
        repij(6,9) = reploc(16)*r31r32*r32r33
     .             + reploc(17)*radd02*r32r33
     .             + reploc(18)*r31r32*radd05
     .             + reploc(19)*(r11r12*r12r13
     .                         + r21r22*r22r23)
     .             + reploc(20)*(radd11*radd13
     .                         + radd12*radd14)
     .             + reploc(21)*(r11r12*r22r23
     .                         + r21r22*r12r13)
     .             + reploc(22)*radd17*(r12r23+r22r13)
!
!       (py pz | px py):
!
        repij(9,6) = repij(6,9)
!
!       (px py | pz pz):
!
        repij(6,10) = reploc(16)*r31r32*r33r33
     .              + reploc(17)*radd02*r33r33
     .              + reploc(18)*radd06*r31r32
     .              + reploc(19)*(r11r12*r13r13
     .                          + r21r22*r23r23)
     .              + reploc(20)*(radd11*r13r33
     .                          + radd12*r23r33)
     .              + reploc(21)*(r11r12*r23r23
     .                          + r21r22*r13r13)
     .              + reploc(22)*radd17*r13r23
!
!       (pz pz | px py):
!
        repij(10,6) = reploc(16)*r33r33*r31r32
     .              + reploc(17)*radd06*r31r32
     .              + reploc(18)*radd02*r33r33
     .              + reploc(19)*(r13r13*r11r12
     .                          + r23r23*r21r22)
     .              + reploc(20)*(r13r33*r11r32
     .                          + r33r13*r31r12
     .                          + r23r33*r21r32
     .                          + r33r23*r31r22)
     .              + reploc(21)*(r13r13*r21r22
     .                          + r23r23*r11r12)
     .              + reploc(22)*(r13r23*r11r22
     .                          + r23r13*r21r12)
!
!       (px pz | px pz):
!
        repij(7,7) = reploc(16)*r31r33*r31r33
     .             + reploc(17)*radd03*r31r33
     .             + reploc(18)*radd03*r31r33
     .             + reploc(19)*(r11r13*r11r13
     .                         + r21r23*r21r23)
     .             + reploc(20)*(radd09**2
     .                         + radd10**2)
     .             + reploc(21)*r11r13*r21r23*2.0d0
     .             + reploc(22)*radd15**2
!
!       (px pz | py py):
!
        repij(7,8) = reploc(16)*r31r33*r32r32
     .             + reploc(17)*radd03*r32r32
     .             + reploc(18)*radd04*r31r33
     .             + reploc(19)*(r11r13*r12r12
     .                         + r21r23*r22r22)
     .             + reploc(20)*(radd09*r12r32
     .                         + radd10*r22r32)
     .             + reploc(21)*(r11r13*r22r22
     .                         + r21r23*r12r12)
     .             + reploc(22)*radd15*r12r22
!
!       (py py | px pz):
!
        repij(8,7) = reploc(16)*r32r32*r31r33
     .             + reploc(17)*radd04*r31r33
     .             + reploc(18)*radd03*r32r32
     .             + reploc(19)*(r12r12*r11r13
     .                         + r22r22*r21r23)
     .             + reploc(20)*(r12r32*r11r33
     .                         + r32r12*r31r13
     .                         + r22r32*r21r33
     .                         + r32r22*r31r23)
     .             + reploc(21)*(r12r12*r21r23
     .                         + r22r22*r11r13)
     .             + reploc(22)*(r12r22*r11r23
     .                         + r22r12*r21r13)
!
!       (px pz | py pz):
!
        repij(7,9) = reploc(16)*r31r33*r32r33
     .             + reploc(17)*radd03*r32r33
     .             + reploc(18)*r31r33*radd05
     .             + reploc(19)*(r11r13*r12r13
     .                         + r21r23*r22r23)
     .             + reploc(20)*(radd09*radd13
     .                         + radd10*r22r33
     .                         + radd10*r32r23)
     .             + reploc(21)*(r11r13*r22r23
     .                         + r21r23*r12r13)
     .             + reploc(22)*radd15*radd16
!
!       (py pz | px pz):
!
        repij(9,7) = repij(7,9)
!
!       (px pz | pz pz):
!
        repij(7,10) = reploc(16)*r31r33*r33r33
     .              + reploc(17)*radd03*r33r33
     .              + reploc(18)*radd06*r31r33
     .              + reploc(19)*(r11r13*r13r13
     .                          + r21r23*r23r23)
     .              + reploc(20)*(radd09*r13r33
     .                          + radd10*r23r33)
     .              + reploc(21)*(r11r13*r23r23
     .                          + r21r23*r13r13)
     .              + reploc(22)*radd15*r13r23
!
!       (pz pz | px pz):
!
        repij(10,7) = reploc(16)*r33r33*r31r33
     .              + reploc(17)*radd06*r31r33
     .              + reploc(18)*radd03*r33r33
     .              + reploc(19)*(r13r13*r11r13
     .                          + r23r23*r21r23)
     .              + reploc(20)*(r13r33*r11r33
     .                          + r33r13*r31r13
     .                          + r23r33*r21r33
     .                          + r33r23*r31r23)
     .              + reploc(21)*(r13r13*r21r23
     .                          + r23r23*r11r13)
     .              + reploc(22)*(r13r23*r11r23
     .                          + r23r13*r21r13)
!
!       (py py | py py):
!
        repij(8,8) = reploc(16)*r32r32*r32r32
     .             + reploc(17)*radd04*r32r32
     .             + reploc(18)*radd04*r32r32
     .             + reploc(19)*(r12r12*r12r12
     .                         + r22r22*r22r22)
     .             + reploc(20)*(r12r32*r12r32
     .                         + r22r32*r22r32)
     .             + reploc(21)*r12r12*r22r22*2.0d0
     .             + reploc(22)*r12r22*r12r22
!
!       (py py | py pz):
!
        repij(8,9) = reploc(16)*r32r32*r32r33
     .             + reploc(17)*radd04*r32r33
     .             + reploc(18)*radd05*r32r32
     .             + reploc(19)*(r12r12*r12r13
     .                         + r22r22*r22r23)
     .             + reploc(20)*(r12r32*r12r33
     .                         + r32r12*r32r13
     .                         + r22r32*r22r33
     .                         + r32r22*r32r23)
     .             + reploc(21)*(r12r12*r22r23
     .                         + r22r22*r12r13)
     .             + reploc(22)*(r12r22*r12r23
     .                         + r22r12*r22r13)
!
!       (py pz | py py):
!
        repij(9,8) = reploc(16)*r32r33*r32r32
     .             + reploc(17)*radd05*r32r32
     .             + reploc(18)*radd04*r32r33
     .             + reploc(19)*(r12r13*r12r12
     .                         + r22r23*r22r22)
     .             + reploc(20)*(radd13*r12r32
     .                         + radd14*r22r32)
     .             + reploc(21)*(r12r13*r22r22
     .                         + r22r23*r12r12)
     .             + reploc(22)*radd16*r12r22
!
!       (py py | pz pz):
!
        repij(8,10) = reploc(16)*r32r32*r33r33
     .              + reploc(17)*radd04*r33r33
     .              + reploc(18)*radd06*r32r32
     .              + reploc(19)*(r12r12*r13r13
     .                          + r22r22*r23r23)
     .              + reploc(20)*(r12r32*r13r33
     .                          + r22r32*r23r33)
     .              + reploc(21)*(r12r12*r23r23
     .                          + r22r22*r13r13)
     .              + reploc(22)*r12r22*r13r23
!
!       (pz pz | py py):
!
        repij(10,8) = reploc(16)*r33r33*r32r32
     .              + reploc(17)*radd06*r32r32
     .              + reploc(18)*radd04*r33r33
     .              + reploc(19)*(r13r13*r12r12
     .                          + r23r23*r22r22)
     .              + reploc(20)*(r13r33*r12r32
     .                          + r23r33*r22r32)
     .              + reploc(21)*(r13r13*r22r22
     .                          + r23r23*r12r12)
     .              + reploc(22)*r13r23*r12r22
!
!       (py pz | py pz):
!
        repij(9,9) = reploc(16)*r32r33*r32r33
     .             + reploc(17)*radd05*r32r33
     .             + reploc(18)*radd05*r32r33
     .             + reploc(19)*(r12r13*r12r13
     .                         + r22r23*r22r23)
     .             + reploc(20)*(radd13**2
     .                         + radd14**2)
     .             + reploc(21)*(r12r13*r22r23)*2.0d0
     .             + reploc(22)*radd16**2
!
!       (py pz | pz pz):
!
        repij(9,10) = reploc(16)*r32r33*r33r33
     .              + reploc(17)*radd05*r33r33
     .              + reploc(18)*radd06*r32r33
     .              + reploc(19)*(r12r13*r13r13
     .                          + r22r23*r23r23)
     .              + reploc(20)*(radd13*r13r33
     .                          + radd14*r23r33)
     .              + reploc(21)*(r12r13*r23r23
     .                          + r22r23*r13r13)
     .              + reploc(22)*radd16*r13r23
!
!       (pz pz | py pz):
!
        repij(10,9) = reploc(16)*r33r33*r32r33
     .              + reploc(17)*radd06*r32r33
     .              + reploc(18)*radd05*r33r33
     .              + reploc(19)*(r13r13*r12r13
     .                          + r23r23*r22r23)
     .              + reploc(20)*(r13r33*r12r33
     .                          + r33r13*r32r13
     .                          + r23r33*r22r33
     .                          + r33r23*r32r23)
     .              + reploc(21)*(r13r13*r22r23
     .                          + r23r23*r12r13)
     .              + reploc(22)*(r13r23*r12r23
     .                          + r23r13*r22r13)
!
!       (pz pz | pz pz):
!
        repij(10,10) = reploc(16)*r33r33*r33r33
     .               + reploc(17)*radd06*r33r33
     .               + reploc(18)*radd06*r33r33
     .               + reploc(19)*(r13r13*r13r13
     .                           + r23r23*r23r23)
     .               + reploc(20)*(r13r33*r13r33
     .                           + r23r33*r23r33)
     .               + reploc(21)*r13r13*r23r23*2.0d0
     .               + reploc(22)*r13r23*r13r23
      endif
      return
      end
