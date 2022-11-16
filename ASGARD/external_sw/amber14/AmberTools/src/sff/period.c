c
c         -- truncated octahedral box
c
                  dwx = fw(1,jn)*boxi(1)
                  dwy = fw(2,jn)*boxi(2)
                  dwz = fw(3,jn)*boxi(3)
                  dux = xwij(1,jn)*boxi(1)
                  duy = xwij(2,jn)*boxi(2)
                  duz = xwij(3,jn)*boxi(3)
                  dux = dux - dnint( dwx )
                  duy = duy - dnint( dwy )
                  duz = duz - dnint( dwz )
                  dwx = dwx - dnint( dwx )
                  dwy = dwy - dnint( dwy )
                  dwz = dwz - dnint( dwz )
                  corr= 0.5 * aint( fothi* (abs( dwx ) +
     .                                     abs( dwy ) +
     .                                     abs( dwz )))
                  dux = dux - sign( corr,dwx )
                  duy = duy - sign( corr,dwy )
                  duz = duz - sign( corr,dwz )
                  xwij(1,jn) = dux * box(1)
                  xwij(2,jn) = duy * box(2)
                  xwij(3,jn) = duz * box(3)
              ! r2(JN) = ONE / (XWIJ(1,JN)**2+XWIJ(2,JN)**2+XWIJ(3,JN)**2)
