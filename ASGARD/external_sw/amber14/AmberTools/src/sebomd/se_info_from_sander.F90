subroutine se_info_from_sander(netcharge, &
    ipolyn, fullscf, screen, wrtscr, prtsub, &
    pme, pmeqm, mullewald, cmewald, chewald, &
    stand, clust, resclust, smallclust, &
    debug, &
    mndo, am1, pm3, rm1, am1d, pm3pddg, &
    nncore, dcbuff1, dcbuff2, &
    diag_routine, dpmaxcrt)

  use sebomd_module, only: sebomd_obj

  implicit none

  integer netcharge

  logical ipolyn
  logical fullscf
  logical screen
  logical wrtscr
  logical prtsub

  logical pme
  logical pmeqm
  logical mullewald
  logical cmewald
  logical chewald

  logical mndo
  logical am1
  logical pm3
  logical rm1
  logical am1d
  logical pm3pddg

  logical debug

  logical stand
  logical clust
  logical resclust
  logical smallclust

  integer nncore
  double precision dcbuff1
  double precision dcbuff2

  integer diag_routine
  double precision dpmaxcrt

  dcbuff1 = sebomd_obj%dbuff1
  dcbuff2 = sebomd_obj%dbuff2

  netcharge = sebomd_obj%charge

  ipolyn = (sebomd_obj%ipolyn.eq.1)
  debug = (sebomd_obj%debugmsg.eq.1)

  diag_routine = sebomd_obj%diag_routine
  dpmaxcrt = sebomd_obj%dpmax

  if (sebomd_obj%hamiltonian.eq.'MNDO') then
    mndo = .true.
  else if (sebomd_obj%hamiltonian.eq.'AM1') then
    am1 = .true.
  else if (sebomd_obj%hamiltonian.eq.'PM3') then
    pm3 = .true.
  else if (sebomd_obj%hamiltonian.eq.'RM1') then
    rm1 = .true.
  else if (sebomd_obj%hamiltonian.eq.'AM1D') then
    am1d = .true.
  else if (sebomd_obj%hamiltonian.eq.'PM3PDDG') then
    pm3pddg = .true.
  else ! default
    pm3 = .true.
  endif

  stand = .false.
  clust = .false.
  resclust = .false.
  smallclust = .false.
  if (sebomd_obj%method.eq.0) then
    stand = .true.
  else if (sebomd_obj%method.eq.1) then
    clust = .true.
  else if (sebomd_obj%method.eq.2) then
    resclust = .true.
  else if (sebomd_obj%method.eq.3) then
    resclust = .true.
    smallclust = .true.
  else ! default
    stand = .true.
  endif

  fullscf = (sebomd_obj%fullscf.eq.1)

  screen = .false.
  wrtscr = .false.
  prtsub = .false.
  if (sebomd_obj%screen.eq.1) then
    screen = .true.
    wrtscr = .true.
  else if (sebomd_obj%screen.eq.2) then
    screen = .true.
    wrtscr = .true.
    prtsub = .true.
  endif

  pme = .false.
  pmeqm = .false.
  mullewald = .false.
  cmewald = .false.
  chewald = .false.

  if (sebomd_obj%longrange == 1) then
     pme = .true.
  elseif (sebomd_obj%longrange == 2) then
     pmeqm = .true.
     mullewald = .true.
  elseif (sebomd_obj%longrange == 3) then
     pmeqm = .true.
     cmewald = .true.
  elseif (sebomd_obj%longrange == 4) then
     pmeqm = .true.
     chewald = .true.
     sebomd_obj%chtype = sebomd_obj%chewald
  else
     continue
  endif

  nncore = sebomd_obj%ncore  

end subroutine se_info_from_sander
