/***=======================================================================***/
/*** Like the nonbonded loop, there are different versions of the bond,    ***/
/*** angle, and dihedral routines.  The different variants can calculate   ***/
/*** forces alone, energies alone, forces and energies, or forces,         ***/
/*** energies, and virial contributions.                                   ***/
/***=======================================================================***/
#if NEEDFORCE == 1
  #if NEEDENERGY == 1
    #if NEEDVIRIAL == 1 
      #define ATTNPFRC AttenuatePairVir
      #define BONDCALC BondVir
      #define ANGLCALC AnglVir
      #define DIHECALC DiheVir
    #else
      #define ATTNPFRC AttenuatePairFrcNrg
      #define BONDCALC BondFrcNrg
      #define ANGLCALC AnglFrcNrg
      #define DIHECALC DiheFrcNrg
    #endif
  #else
    #define ATTNPFRC AttenuatePairFrc
    #define BONDCALC BondFrc
    #define ANGLCALC AnglFrc
    #define DIHECALC DiheFrc
  #endif
#else
  #define ATTNPFRC AttenuatePairNrg
  #define BONDCALC BondNrg
  #define ANGLCALC AnglNrg
  #define DIHECALC DiheNrg
#endif

/***=======================================================================***/
/*** AttenuatePairForce: routine to modify the pair force between atoms in ***/
/***                     order to account for exclusion of basic van-der   ***/
/***                     Waals or electrostatic interactions.  After any   ***/
/***                     necessary corrections are added, the forces,      ***/
/***                     energies, and virials are accumulated.            ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***=======================================================================***/
#if NEEDFORCE == 1
#if NEEDENERGY == 1
void ATTNPFRC(prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, int atmA, int atmB,
	      double dx, double dy, double dz, double fmag, double* afrc,
	      double* bfrc, Energy *sysUV, double elec14fac, double lj14fac)
#else
void ATTNPFRC(prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, int atmA, int atmB,
	      double dx, double dy, double dz, double fmag, double* afrc,
	      double* bfrc, double elec14fac, double lj14fac)
#endif
#else
void ATTNPFRC(prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, int atmA, int atmB,
	      double dx, double dy, double dz, Energy *sysUV, double elec14fac,
	      double lj14fac)
#endif
{
  int irc, irh;
  double r, r2, invr, invr2, qq;
#if NEEDFORCE == 1
  CSpln *CdSD, *HdSD;
#endif
#if NEEDENERGY == 1
  CSpln *CSD, *HSD;
#endif

  /*** (Re)compute r and related quantities ***/
  r2 = dx*dx + dy*dy + dz*dz;
  r = sqrt(r2);
  invr = 1.0/r;
  invr2 = invr*invr;

  /*** Correct 1:2 electrostatics.  All electrostatic interactions are ***/
  /*** counted in the nonbonded routine (direct and reciprocal space), ***/
  /*** so exclusions or 1:4 screened interactions must be addressed.   ***/
  /*** However, electrostatic direct-space interactions are computed   ***/
  /*** with a coarse lookup table, which is inaccurate at short range. ***/
  /*** If an electrostatic interaction was calculated at short range,  ***/
  /*** we subtract off the contributions from the coarse lookup table  ***/
  /*** and replace them with calculations from a more accurate fine    ***/
  /*** lookup table before subtracting the analytic qq/r interaction.  ***/
  qq = tp->Charges[atmA] * tp->Charges[atmB];
  if (r2 < MINNB2) {
    irc = r2*Cfrc->ivdr;
#if NEEDFORCE == 1
    CdSD = Cfrc->dSD;
    fmag -= qq*(((CdSD[irc].A*r2 + CdSD[irc].B)*r2 + CdSD[irc].C)*r2 +
		CdSD[irc].D);
#endif
#if NEEDENERGY == 1
    CSD = Cfrc->SD;
    sysUV->delec -= qq*(((CSD[irc].A*r2 + CSD[irc].B)*r2 + CSD[irc].C)*r2 +
			CSD[irc].D);
#endif
    irh = r2*Hfrc->ivdr;
#if NEEDFORCE == 1
    HdSD = Hfrc->dSD;
    fmag += qq*(((HdSD[irh].A*r2 + HdSD[irh].B)*r2 + HdSD[irh].C)*r2 +
		HdSD[irh].D);
#endif
#if NEEDENERGY == 1
    HSD = Hfrc->SD;
    sysUV->delec += qq*(((HSD[irh].A*r2 + HSD[irh].B)*r2 + HSD[irh].C)*r2 +
			HSD[irh].D);
#endif
  }
#if NEEDFORCE == 1
  fmag += elec14fac*BIOQ*qq*invr*invr2;
#endif
#if NEEDENERGY == 1
  sysUV->delec -= elec14fac*BIOQ*qq*invr;
#endif

  /*** Correct 1:2 van-der Waals interactions ***/
  if (r2 > MINNB2) {

    int ljA, ljB;
    double invr4, invr6;

    ljA = tp->LJIdx[atmA];
    ljB = tp->LJIdx[atmB];
    if (ljA >= 0 && ljB >= 0) {
      invr4 = invr2*invr2;
      invr6 = invr4*invr2;
#if NEEDFORCE == 1
      fmag -= lj14fac*invr4*invr4*(tp->LJftab.map[ljA][2*ljB]*invr6 +
				   tp->LJftab.map[ljA][2*ljB+1]);
#endif
#if NEEDENERGY == 1
      sysUV->vdw12 -= lj14fac*invr6*invr6*tp->LJutab.map[ljA][2*ljB];
      sysUV->vdw6 -= lj14fac*invr6*tp->LJutab.map[ljA][2*ljB+1];
#endif
    }
  }

#if NEEDFORCE == 1
  /*** Accumulate the force ***/
  afrc[0] += fmag*dx;
  afrc[1] += fmag*dy;
  afrc[2] += fmag*dz;
  bfrc[0] -= fmag*dx;
  bfrc[1] -= fmag*dy;
  bfrc[2] -= fmag*dz;

#if NEEDVIRIAL == 1
  /*** Accumulate the stress tensor ***/
  sysUV->Vir[0] += dx*fmag*dx;
  sysUV->Vir[1] += dx*fmag*dy;
  sysUV->Vir[2] += dx*fmag*dz;
  sysUV->Vir[3] += dy*fmag*dx;
  sysUV->Vir[4] += dy*fmag*dy;
  sysUV->Vir[5] += dy*fmag*dz;
  sysUV->Vir[6] += dz*fmag*dx;
  sysUV->Vir[7] += dz*fmag*dy;
  sysUV->Vir[8] += dz*fmag*dz;
#endif
#endif
}

/***=======================================================================***/
/*** BondFrc: computes the force on two atoms due to a bonded interaction. ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   [ab]ptr:    the coordinates for atoms A or B                        ***/
/***   [ab]frc:    the forces for atoms A or B                             ***/
/***   bcom:       the bond command                                        ***/
/***   tp:         the topology                                            ***/
/***   Cfrc:       the "coarse" electrostatic force/energy spline table    ***/
/***   Hfrc:       the "fine" electrostatic force/energy spline table      ***/
/***   sysUV:      the system energy and virial (state information)        ***/
/***=======================================================================***/
#if NEEDFORCE == 1
#if NEEDENERGY == 1
void BONDCALC(double *aptr, double *bptr, double *afrc, double *bfrc,
	      bondcomm *bcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
	      Energy *sysUV)
#else
void BONDCALC(double *aptr, double *bptr, double *afrc, double *bfrc,
              bondcomm *bcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc)
#endif
#else
void BONDCALC(double *aptr, double *bptr, bondcomm *bcom, prmtop *tp,
              FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV)
#endif
{
  double r, r2, dx, dy, dz, dl;
#if NEEDFORCE == 1
  double fmag;
#endif

  /*** Compute displacement ***/
  dx = bptr[0] - aptr[0];
  dy = bptr[1] - aptr[1];
  dz = bptr[2] - aptr[2];

  /*** Accumulate the bond force and energy ***/
  r2 = dx*dx + dy*dy + dz*dz;
  r = sqrt(r2);
  dl = tp->BParam[bcom->t].l0 - r;
#if NEEDFORCE == 1
  fmag = -2.0*tp->BParam[bcom->t].K*dl/r;
#endif
#if NEEDENERGY == 1
  sysUV->bond += tp->BParam[bcom->t].K*dl*dl;
  sysUV->BondUdc[3*bcom->t] += r*r;
  sysUV->BondUdc[3*bcom->t+1] += r*tp->BParam[bcom->t].l0;
  sysUV->BondUdc[3*bcom->t+2] += pow(tp->BParam[bcom->t].l0, 2.0);
#endif

#if NEEDFORCE == 1
#if NEEDVIRIAL == 0
#if NEEDENERGY == 1
  AttenuatePairFrcNrg(tp, Cfrc, Hfrc, bcom->a, bcom->b, dx, dy, dz, fmag,
		      afrc, bfrc, sysUV, 1.0, 1.0);
#else
  AttenuatePairFrc(tp, Cfrc, Hfrc, bcom->a, bcom->b, dx, dy, dz, fmag,
		   afrc, bfrc, 1.0, 1.0);
#endif
#else
  AttenuatePairVir(tp, Cfrc, Hfrc, bcom->a, bcom->b, dx, dy, dz, fmag,
                   afrc, bfrc, sysUV, 1.0, 1.0);
#endif
#else
  AttenuatePairNrg(tp, Cfrc, Hfrc, bcom->a, bcom->b, dx, dy, dz, sysUV, 1.0,
		   1.0);
#endif
}

/***=======================================================================***/
/*** AnglFrc: computes the force on two atoms due to an angle interaction. ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   [abc]ptr:   the coordinates for atoms A, B, or C                    ***/
/***   [abc]frc:   the forces for atoms A, B, or C                         ***/
/***   acom:       the angle command                                       ***/
/***   tp:         the topology                                            ***/
/***   Cfrc:       the "coarse" electrostatic force/energy spline table    ***/
/***   Hfrc:       the "fine" electrostatic force/energy spline table      ***/
/***   sysUV:      the system energy and virial (state information)        ***/
/***=======================================================================***/
#if NEEDFORCE == 1
#if NEEDENERGY == 1
void ANGLCALC(double *aptr, double *bptr, double *cptr, double *afrc,
	      double *bfrc, double *cfrc, anglcomm *acom, prmtop *tp,
	      FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV)
#else
void ANGLCALC(double *aptr, double *bptr, double *cptr, double *afrc,
	      double *bfrc, double *cfrc, anglcomm *acom, prmtop *tp,
	      FrcTab *Cfrc, FrcTab *Hfrc)
#endif
#else
void ANGLCALC(double *aptr, double *bptr, double *cptr, anglcomm *acom,
              prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV)
#endif
{
  int i;
  double costheta, theta, mgba, mgbc;
#if NEEDFORCE == 1
  double dA, mbabc, sqba, sqbc, adf, cdf;
#endif
  double invbabc, dtheta;
  double ac[3], ba[3], bc[3];

  /*** Compute displacements ***/
  for (i = 0; i < 3; i++) {
    ba[i] = aptr[i] - bptr[i];
    bc[i] = cptr[i] - bptr[i];
    ac[i] = cptr[i] - aptr[i];
  }

  /*** On to the angle force computation ***/
  mgba = ba[0]*ba[0] + ba[1]*ba[1] + ba[2]*ba[2];
  mgbc = bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2];
  invbabc = 1.0/sqrt(mgba*mgbc);
  costheta = (ba[0]*bc[0] + ba[1]*bc[1] + ba[2]*bc[2]) * invbabc;
  costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
  theta = acos(costheta);
  dtheta = theta - tp->AParam[acom->t].th0;
#if NEEDFORCE == 1
  dA = -2.0*tp->AParam[acom->t].K*dtheta /
    sqrt(1.0 - costheta*costheta);
  sqba = dA/mgba;
  sqbc = dA/mgbc;
  mbabc = dA * invbabc;

  /*** Accumulate the angle forces and stress tensor ***/
  for (i = 0; i < 3; i++) {
    adf = bc[i]*mbabc - costheta*ba[i]*sqba;
    cdf = ba[i]*mbabc - costheta*bc[i]*sqbc;
    afrc[i] -= adf;
    bfrc[i] += adf + cdf;
    cfrc[i] -= cdf;
#if NEEDVIRIAL == 1
    sysUV->Vir[0+i] += ba[0]*adf + bc[0]*cdf;
    sysUV->Vir[3+i] += ba[1]*adf + bc[1]*cdf;
    sysUV->Vir[6+i] += ba[2]*adf + bc[2]*cdf;
#endif
  }
#endif

#if NEEDENERGY == 1
  /*** Accumulate the angle energy ***/
  sysUV->angl += tp->AParam[acom->t].K*dtheta*dtheta;
#endif

  /*** Bail out if this angle's A:C interactions do not need to be ***/
  /*** excluded (i.e. they were already excluded by some bond)     ***/
  if (acom->excl == 0) {
    return; 
  }

#if NEEDFORCE == 1
#if NEEDVIRIAL == 0
#if NEEDENERGY == 1
  AttenuatePairFrcNrg(tp, Cfrc, Hfrc, acom->a, acom->c, ac[0], ac[1], ac[2],
		      0.0, afrc, cfrc, sysUV, 1.0, 1.0);
#else
  AttenuatePairFrc(tp, Cfrc, Hfrc, acom->a, acom->c, ac[0], ac[1], ac[2], 0.0,
                   afrc, cfrc, 1.0, 1.0);
#endif
#else
  AttenuatePairVir(tp, Cfrc, Hfrc, acom->a, acom->c, ac[0], ac[1], ac[2], 0.0,
                   afrc, cfrc, sysUV, 1.0, 1.0);
#endif
#else
  AttenuatePairNrg(tp, Cfrc, Hfrc, acom->a, acom->c, ac[0], ac[1], ac[2],
		   sysUV, 1.0, 1.0);
#endif

}

/***=======================================================================***/
/*** DiheFrc: computes the force on four atoms due to a dihedral motion.   ***/
/***          This also applies to improper dihedrals.  Standard and       ***/
/***          improper dihedrals are treated using the following notation: ***/
/***                                                                       ***/
/***             Standard Dihedral          Improper Dihedral              ***/
/***                                                                       ***/
/***                       D                           D                   ***/
/***                      /                           /                    ***/
/***                  B--C                        B--C                     ***/
/***                 /                                \                    ***/
/***                A                                  A                   ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   [abcd]ptr:   the coordinates for atoms A, B, C, or D                ***/
/***   [abcd]frc:   the forces for atoms A, B, C, or D                     ***/
/***   hcom:        the dihedral command, giving details on all Fourier    ***/
/***                series terms and a flag to tell whether this is an     ***/
/***                improper or a standard dihedral                        ***/
/***   hdef:        the array of dihedral definitions                      ***/
/***=======================================================================***/
#if NEEDFORCE == 1
#if NEEDENERGY == 1
void DIHECALC(double *aptr, double *bptr, double *cptr, double *dptr,
	      double *afrc, double *bfrc, double *cfrc, double *dfrc,
	      dihecomm *hcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
	      Energy *sysUV)
#else
void DIHECALC(double *aptr, double *bptr, double *cptr, double *dptr,
	      double *afrc, double *bfrc, double *cfrc, double *dfrc,
	      dihecomm *hcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc)
#endif
#else
void DIHECALC(double *aptr, double *bptr, double *cptr, double *dptr,
	      dihecomm *hcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
	      Energy *sysUV)
#endif
{
  int i, idx;
  double theta, costheta;
#if NEEDFORCE == 1
  double fr, fa, fb1, fc1, fb2, fc2, fd, isinb2, isinc2;
  double mgab, mgbc, mgcd, invab, invbc, invcd, invabc, invbcd, cosb, cosc;
#endif
  double sangle;
  double ab[3], bc[3], cd[3], crabbc[3], crbccd[3], scr[3];
  dihedef *hdef;

  /*** Compute displacements ***/
  for (i = 0; i < 3; i++) {
    ab[i] = bptr[i] - aptr[i];
    bc[i] = cptr[i] - bptr[i];
    cd[i] = dptr[i] - cptr[i];
  }
  CrossP(ab, bc, crabbc);
  CrossP(bc, cd, crbccd);
  costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
  costheta /=
    sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
	 (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  CrossP(crabbc, crbccd, scr);
  costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
  if (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0) {
    theta = acos(costheta);
  }
  else {
    theta = -acos(costheta);
  }

  /*** Compute the magnitude of the force and accumulate the energy ***/
#if NEEDFORCE == 1
  fr = 0.0;
#endif
  hdef = tp->HParam;
  for (i = 0; i < hcom->nt; i++) {
    idx = hcom->t[i];
    sangle = hdef[idx].N*theta - hdef[idx].Phi;
#if NEEDFORCE == 1
    fr += hdef[idx].K * hdef[idx].N * sin(sangle);
#endif
#if NEEDENERGY == 1
    sysUV->dihe += hdef[idx].K * (1.0 + cos(sangle));
#endif
  }

#if NEEDFORCE == 1
  /*** Other pre-computations ***/
  mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
  invab = 1.0/mgab;
  mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
  invbc = 1.0/mgbc;
  mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
  invcd = 1.0/mgcd;
  cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2])*invab*invbc;
  isinb2 = (cosb*cosb < 0.9999) ? 1.0/(1.0 - cosb*cosb) : 0.0;
  cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2])*invbc*invcd;
  isinc2 = (cosc*cosc < 0.9999) ? 1.0/(1.0 - cosc*cosc) : 0.0;
  isinb2 *= fr;
  isinc2 *= fr;
  invabc = invab*invbc;
  invbcd = invbc*invcd;
  for (i = 0; i < 3; i++) {
    crabbc[i] *= invabc;
    crbccd[i] *= invbcd;
  }

  /*** Transform the dihedral forces to cartesian coordinates ***/
  fa = -invab * isinb2;
  fb1 = (mgbc - mgab*cosb) * invabc * isinb2;
  fb2 = cosc * invbc * isinc2;
  fc1 = (mgbc - mgcd*cosc) * invbcd * isinc2;
  fc2 = cosb * invbc * isinb2;
  fd = -invcd * isinc2;
  for (i = 0; i < 3; i++) {
    afrc[i] += crabbc[i] * fa;
    bfrc[i] += fb1 * crabbc[i] - fb2 * crbccd[i];
    cfrc[i] += -fc1 * crbccd[i] + fc2 * crabbc[i];
    dfrc[i] += -fd * crbccd[i];
  }
#endif

  /*** Evaluate 1-4 interactions ***/
  if (hcom->eval14 == 1) {
#if NEEDFORCE == 1
#if NEEDVIRIAL == 0
#if NEEDENERGY == 1
  AttenuatePairFrcNrg(tp, Cfrc, Hfrc, hcom->a, hcom->d, dptr[0] - aptr[0],
		      dptr[1] - aptr[1], dptr[2] - aptr[2], 0.0, afrc, dfrc,
		      sysUV, hcom->scee, hcom->scnb);
#else
  AttenuatePairFrc(tp, Cfrc, Hfrc, hcom->a, hcom->d, dptr[0] - aptr[0],
		   dptr[1] - aptr[1], dptr[2] - aptr[2], 0.0, afrc, dfrc,
		   hcom->scee, hcom->scnb);
#endif
#else
  AttenuatePairVir(tp, Cfrc, Hfrc, hcom->a, hcom->d, dptr[0] - aptr[0],
		   dptr[1] - aptr[1], dptr[2] - aptr[2], 0.0, afrc, dfrc,
		   sysUV, hcom->scee, hcom->scnb);
#endif
#else
  AttenuatePairNrg(tp, Cfrc, Hfrc, hcom->a, hcom->d, dptr[0] - aptr[0],
		   dptr[1] - aptr[1], dptr[2] - aptr[2], sysUV, hcom->scee,
		   hcom->scnb);
#endif
  }
}

#undef BONDCALC
#undef ANGLCALC
#undef DIHECALC
#undef ATTNPFRC
