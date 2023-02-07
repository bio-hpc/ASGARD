!<compile=optimized>
#include "../include/dprec.fh"

module les_data

public

!  parameters for LES:

integer maxles,maxlestyp,maxlesadj
parameter (maxles=500000)
parameter (maxlestyp=100)
parameter (maxlesadj=3000000)

integer BC_LESR,BC_LESI
parameter( BC_LESI=1+maxles*3+maxlesadj*2+1)
parameter( BC_LESR=maxlestyp*maxlestyp+1)

_REAL_ lesfac(maxlestyp*maxlestyp),lfac

! for separate LES and non-LES temperature coupling

_REAL_ ekinles0,temp0les,rndfles,sdfacles
_REAL_ scaltles,tempsules,ekeles,rsdles
_REAL_ ekmhles,ekphles

common/lesr/lesfac,temp0les

integer ileslst(maxlesadj),jleslst(maxlesadj),nlesadj
integer lestyp(maxles),nlesty,lestmp,cnum(maxles),subsp(maxles)

! this one is 1+maxles*3+maxlesadj*2+1

common/lesi/nlesty,lestyp,cnum,subsp,ileslst,jleslst,nlesadj

! some PME variables
! these are communicated in places other than parallel.f! but the sizes should
! not change

_REAL_ eeles,les_vir(3,3)
common/ewlescomm/eeles,les_vir

! some partial REM variables
_REAL_ elesa,elesb,elesd,elesp
common/eprem/elesa,elesb,elesd,elesp

end module les_data
