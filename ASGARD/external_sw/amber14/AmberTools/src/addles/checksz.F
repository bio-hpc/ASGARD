        subroutine checksz( error )
c
c checks current array sizes against limits.
c
        
c error is true if a limit is exceeded; false if none are exceeded
        logical error

#include "SIZE.h"
#include "TOP.h"
#include "MISC.h"
c
c
        write (6,*) 'Checking topology sizes against compiled limits'
c
c check against maximum sizes in SIZE.BLOCK
c
        error=.false.
c
        if (natom.gt.maxnatom) then
                write (6,*) 'exceeded MAXNATOM in SIZE.BLOCK.'
                write (6,*) 'MAXNATOM = ',maxnatom
                write (6,*) 'NATOM = ',natom
                error=.true.
        endif
        if (ntypes.gt.maxntypes) then
                write (6,*) 'exceeded MAXNTYPES in SIZE.BLOCK.'
                write (6,*) 'MAXNTYPES = ',maxntypes
                write (6,*) 'ntypes = ',ntypes
                error=.true.
        endif
        if (nres.gt.maxnres) then
                write (6,*) 'exceeded MAXNRES in SIZE.BLOCK.'
                write (6,*) 'MAXNRES = ',maxnres
                write (6,*) 'nres = ',nres
                error=.true.
        endif
c
c bond, angle dihedral types
c
        if (numbnd.gt.maxbndt) then
                write (6,*) 'exceeded MAXBNDT in SIZE.BLOCK.'
                write (6,*) 'MAXBNDT = ',maxbndt
                write (6,*) 'numbnd = ',numbnd
                error=.true.
        endif
        if (numang.gt.maxangt) then
                write (6,*) 'exceeded MAXANGT in SIZE.BLOCK.'
                write (6,*) 'MAXANGT = ',maxangt
                write (6,*) 'numang = ',numang
                error=.true.
        endif
        if (nptra.gt.maxdiht) then
                write (6,*) 'exceeded MAXDIHT in SIZE.BLOCK.'
                write (6,*) 'MAXDIHT = ',maxdiht
                write (6,*) 'nptra = ',nptra
                error=.true.
        endif
        if (nptra.gt.maxdiht) then
                write (6,*) 'exceeded MAXDIHT in SIZE.BLOCK.'
                write (6,*) 'MAXDIHT = ',maxdiht
                write (6,*) 'nptra = ',nptra
                error=.true.
        endif
c
c actual bonds, angles, dihedrals
c
        if (nbper.gt.maxbnd) then
                write (6,*) 'exceeded MAXBND in SIZE.BLOCK.'
                write (6,*) 'MAXBND = ',maxbnd
                write (6,*) 'nbper = ',nbper
                error=.true.
        endif
        if (nbona.gt.maxbnd) then
                write (6,*) 'exceeded MAXBND in SIZE.BLOCK.'
                write (6,*) 'MAXBND = ',maxbnd
                write (6,*) 'nbona = ',nbona
                error=.true.
        endif
        if (nbonh.gt.maxbnd) then
                write (6,*) 'exceeded MAXBND in SIZE.BLOCK.'
                write (6,*) 'MAXBND = ',maxbnd
                write (6,*) 'nbonh = ',nbonh
                error=.true.
        endif
        if (ngper.gt.maxang) then
                write (6,*) 'exceeded MAXANG in SIZE.BLOCK.'
                write (6,*) 'MAXANG = ',maxang
                write (6,*) 'ngper = ',ngper
                error=.true.
        endif
        if (ntheta.gt.maxang) then
                write (6,*) 'exceeded MAXANG in SIZE.BLOCK.'
                write (6,*) 'MAXANG = ',maxang
                write (6,*) 'ntheta = ',ntheta
                error=.true.
        endif
        if (ntheth.gt.maxang) then
                write (6,*) 'exceeded MAXANG in SIZE.BLOCK.'
                write (6,*) 'MAXANG = ',maxang
                write (6,*) 'ntheth = ',ntheth
                error=.true.
        endif
        if (ndper.gt.maxdih) then
                write (6,*) 'exceeded MAXDIH in SIZE.BLOCK.'
                write (6,*) 'MAXDIH = ',maxdih
                write (6,*) 'ndper = ',ndper
                error=.true.
        endif
        if (nphia.gt.maxdih) then
                write (6,*) 'exceeded MAXDIH in SIZE.BLOCK.'
                write (6,*) 'MAXDIH = ',maxdih
                write (6,*) 'nphia = ',nphia
                error=.true.
        endif
        if (nphih.gt.maxdih) then
                write (6,*) 'exceeded MAXDIH in SIZE.BLOCK.'
                write (6,*) 'MAXDIH = ',maxdih
                write (6,*) 'nphih = ',nphih
                error=.true.
        endif
        if (natyp.gt.maxnatyp) then
                write (6,*) 'exceeded MAXNATYp in SIZE.BLOCK.'
                write (6,*) 'MAXNATYp = ',maxnatyp
                write (6,*) 'natyp = ',natyp
                error=.true.
        endif
        if (next.gt.maxnext) then
                write (6,*) 'exceeded MAXNEXT in SIZE.BLOCK.'
                write (6,*) 'MAXNEXT = ',maxnext
                write (6,*) 'next = ',next
                error=.true.
        endif
        if (nphb.gt.maxnphb) then
                write (6,*) 'exceeded MAXNPHB in SIZE.BLOCK.'
                write (6,*) 'MAXNPHB = ',maxnphb
                write (6,*) 'nphb = ',nphb
                error=.true.
        endif
        if (nspm.gt.maxnspm) then
                write (6,*) 'exceeded MAXNSPM in SIZE.BLOCK.'
                write (6,*) 'MAXNSPM = ',maxnspm
                write (6,*) 'nspm = ',nspm
                error=.true.
        endif
c
        if (error) then
                write (6,*) 'Change SIZE.BLOCK and recompile.'
        endif

        return
c
        end
