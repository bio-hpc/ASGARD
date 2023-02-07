#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mdgxVector.h"
#include "Matrix.h"
#include "ParamFit.h"
#include "Topology.h"
#include "CrdManip.h"
#include "Parse.h"
#include "ChargeFit.h"
#include "VirtualSites.h"
#include "Manual.h"
#include "Nonbonded.h"
#include "Trajectory.h"
#include "Macros.h"
#include "ParamRead.h"
#include "ParamFit.h"

#include "CrdManipDS.h"

/***=======================================================================***/
/*** SetParmOutputFlags: make decisions regarding whether to print the     ***/
/***                     parameters to output.  By default, all parameters ***/
/***                     which were present in the original parm.dat and   ***/
/***                     frcmod files are printed to output, but if there  ***/
/***                     are cloned parameters then they will be printed   ***/
/***                     only if they are found in one of the systems      ***/
/***                     under consideration.  By turning up the report    ***/
/***                     level, all cloned parameters can appear in the    ***/
/***                     output, and by turning down the report level only ***/
/***                     parameters needed for the systems involved in the ***/
/***                     fitting exercise will make it into the fit.       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***=======================================================================***/
static void SetParmOutputFlags(prmset *mp)
{
  int i, j, k, atmid, replvl;
  prmtop *tp;

  replvl = mp->reportall;

  /*** Atoms ***/
  for (i = 0; i < mp->natom; i++) {
    if (replvl == 2 || (replvl == 1 && mp->atoms[i].dup == 0)) {
      mp->atoms[i].inreport = 1;
    }
    else {
      mp->atoms[i].inreport = 0;
    }
  }
  for (i = 0; i < mp->nunisys; i++) {
    for (j = 0; j < mp->nconf; j++) {
      if (mp->conf[j].GroupNum == i) {
        tp = mp->conf[j].tp;
        for (k = 0; k < tp->natom; k++) {
          atmid = CrossRefAtomType(mp, &tp->AtomTypes[4*k]);
          mp->atoms[atmid].inreport = 1;
        }
	break;
      }
    }
  }

  /*** Bonds, angles, dihedrals ***/
  for (i = 0; i < mp->nbond; i++) {
    if (replvl == 2 || (replvl == 1 && mp->bonds[i].dup == 0)) {
      mp->bonds[i].inreport = 1;
    }
    else {
      mp->bonds[i].inreport = 0;
    }
  }
  for (i = 0; i < mp->nangl; i++) {
    if (replvl == 2 || (replvl == 1 && mp->angls[i].dup == 0)) {
      mp->angls[i].inreport = 1;
    }
    else {
      mp->angls[i].inreport = 0;
    }
  }
  for (i = 0; i < mp->ntor; i++) {
    if (replvl == 2 || (replvl == 1 && mp->torsions[i].dup == 0)) {
      mp->torsions[i].inreport = 1;
    }
    else {
      mp->torsions[i].inreport = 0;
    }
  }
  for (i = 0; i < mp->nhb1012; i++) {
    if (replvl == 2 || (replvl == 1 && mp->hb1012[i].dup == 0)) {
      mp->hb1012[i].inreport = 1;
    }
    else {
      mp->hb1012[i].inreport = 0;
    }
  }

  /*** Determine whether these bonds, angles, or dihedrals ***/
  /*** are present in ANY system used in this fit.         ***/
  for (i = 0; i < mp->nunisys; i++) {
    for (j = 0; j < mp->nconf; j++) {
      if (mp->conf[j].GroupNum == i) {
        for (k = 0; k < mp->conf[j].bmap.nbond; k++) {
	  mp->bonds[mp->conf[j].bmap.id[k].key].inreport = 1;
        }
        for (k = 0; k < mp->conf[j].amap.nangl; k++) {
	  mp->angls[mp->conf[j].amap.id[k].key].inreport = 1;
        }
        for (k = 0; k < mp->conf[j].hmap.ntterm; k++) {
	  mp->torsions[mp->conf[j].hmap.id[k].key].inreport = 1;
	}
	break;
      }
    }
  }
}

/***=======================================================================***/
/*** PrintParmAtoms: print the masses involved in a frcmod file for this   ***/
/***                   parameter fit.                                      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   outp:     the output frcmod file (pointer)                          ***/
/***=======================================================================***/
static void PrintParmAtoms(prmset *mp, FILE *outp)
{
  int i;

  for (i = 0; i < mp->natom; i++) {
    if (mp->atoms[i].inreport == 1) {
      fprintf(outp, "%.2s               %10.4lf  %10.4lf  %s\n",
	      mp->atoms[i].atype, mp->atoms[i].mass, mp->atoms[i].apol,
	      mp->atoms[i].comment);
    }
  }
  fprintf(outp, "\n");
}

/***=======================================================================***/
/*** PrintParmBond: this procedure is much like PrintParmAtoms above.      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   order:    the order of bonded terms to seek out (2 = bonds,         ***/
/***             3 = angles, 4 = dihedrals)                                ***/
/***   outp:     the output file                                           ***/
/***   x:        the vector of fitted parameters                           ***/
/***=======================================================================***/
static void PrintParmBond(prmset *mp, int order, double* x, FILE *outp)
{
  int i, ilim, slen;
  double stiffness;
  char tmpline[128];

  /*** Print bonds, angles, and dihedrals to the new parameter file ***/
  if (order == 2) ilim = mp->nbond;
  else if (order == 3) ilim = mp->nangl;
  else if (order == 4) ilim = mp->ntor;
  else if (order == 5) ilim = mp->ntor;
  else if (order == 9) ilim = mp->nhb1012;
  for (i = 0; i < ilim; i++) {
    if (order == 2 && mp->bonds[i].inreport == 1) {
      if (mp->bonds[i].fitcol >= 0) {
	stiffness = x[mp->bonds[i].fitcol];
	if (strncmp(mp->bonds[i].comment, "Branched from ", 14) == 0) {
	  slen = strlen(mp->bonds[i].comment);
	  mp->bonds[i].comment[slen] = ',';
	  mp->bonds[i].comment[slen+1] = ' ';
	  mp->bonds[i].comment[slen+2] = '\0';
	  slen += 2;
	}
	else {
	  slen = 0;
	}
        strcpy(mp->bonds[i].comment + slen, mp->icomm);
      }
      else {
	stiffness = mp->bonds[i].K;
      }
      sprintf(tmpline, "%.2s-%.2s            %10.4lf  %10.4lf  %s\n",
	      mp->bonds[i].atype, mp->bonds[i].btype, stiffness,
	      mp->bonds[i].l0, mp->bonds[i].comment);
      fprintf(outp, "%s", tmpline);
    }
    if (order == 3 && mp->angls[i].inreport == 1) {
      if (mp->angls[i].fitcol >= 0) {
	stiffness = x[mp->angls[i].fitcol];
	if (strncmp(mp->angls[i].comment, "Branched from ", 14) == 0) {
	  slen = strlen(mp->angls[i].comment);
	  mp->angls[i].comment[slen] = ',';
	  mp->angls[i].comment[slen+1] = ' ';
	  mp->angls[i].comment[slen+2] = '\0';
	  slen += 2;
	}
	else {
	  slen = 0;
	}
        strcpy(mp->angls[i].comment + slen, mp->icomm);
      }
      else {
	stiffness = mp->angls[i].K;
      }
      sprintf(tmpline, "%.2s-%.2s-%.2s         %10.4lf  %10.2lf  %s\n",
              mp->angls[i].atype, mp->angls[i].btype, mp->angls[i].ctype,
              stiffness, (180.0/PI)*mp->angls[i].th0, mp->angls[i].comment);
      fprintf(outp, "%s", tmpline);
    }
    if (order == 4 && mp->torsions[i].impr == 0 &&
	mp->torsions[i].inreport == 1) {
      if (mp->torsions[i].fitcol >= 0) {
	stiffness = x[mp->torsions[i].fitcol];
	if (strncmp(mp->torsions[i].comment, "Branched from ", 14) == 0) {
	  slen = strlen(mp->torsions[i].comment);
	  mp->torsions[i].comment[slen] = ',';
	  mp->torsions[i].comment[slen+1] = ' ';
	  mp->torsions[i].comment[slen+2] = '\0';
	  slen += 2;
	}
	else {
	  slen = 0;
	}
        strcpy(mp->torsions[i].comment + slen, mp->icomm);
      }
      else {
	stiffness = mp->torsions[i].K;
      }
      sprintf(tmpline, "%.2s-%.2s-%.2s-%.2s   1  %10.5lf %6.1lf %4.1lf  "
	      "%s\n", mp->torsions[i].atype, mp->torsions[i].btype,
              mp->torsions[i].ctype, mp->torsions[i].dtype, stiffness,
	      mp->torsions[i].phase*180.0/PI,
	      mp->torsions[i].singlet * mp->torsions[i].pn,
	      mp->torsions[i].comment);
      fprintf(outp, "%s", tmpline);
    }
    if (order == 5 && mp->torsions[i].impr == 1 &&
	mp->torsions[i].inreport == 1) {
      if (mp->torsions[i].fitcol >= 0) {
	stiffness = x[mp->torsions[i].fitcol];
        strcpy(mp->torsions[i].comment, mp->icomm);
      }
      else {
	stiffness = mp->torsions[i].K;
      }
      sprintf(tmpline, "%.2s-%.2s-%.2s-%.2s      %10.5lf %6.1lf %4.1lf  "
	      "%s\n", mp->torsions[i].atype, mp->torsions[i].btype,
              mp->torsions[i].ctype, mp->torsions[i].dtype, stiffness,
	      mp->torsions[i].phase*180.0/PI,
              mp->torsions[i].singlet * mp->torsions[i].pn,
	      mp->torsions[i].comment);
      fprintf(outp, "%s", tmpline);
    }
    if (order == 9 && mp->hb1012[i].inreport == 1) {

      /*** Currently there is no support for ***/
      /*** fitting H-bond 10-12 potentials.  ***/
      sprintf(tmpline, "  %.2s  %.2s         %10.4lf  %10.4lf  %s\n",
	      mp->hb1012[i].atype, mp->hb1012[i].btype, mp->hb1012[i].Aterm,
	      mp->hb1012[i].Bterm, mp->hb1012[i].comment);
      fprintf(outp, "%s", tmpline);
    }
  }
  fprintf(outp, "\n");
}

/***=======================================================================***/
/*** MapSamplingIndices: while the fitcol indices in each bond, angle, or  ***/
/***                     torsion map each adjustable term into the fitting ***/
/***                     matrix, the samprow indices needed to map each    ***/
/***                     term into the distinct sampling matrices.         ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   order:    the order of bonded terms to seek out (2 = bonds,         ***/
/***             3 = angles, 4 = dihedrals)                                ***/
/***=======================================================================***/
static void MapSamplingIndices(prmset *mp, int order)
{
  int i, j, ilim;

  /*** Map the indices of the sampling matrices ***/
  j = 0;
  ilim = (order == 2) ? mp->nbond : (order == 3) ? mp->nangl : mp->ntor;
  for (i = 0; i < ilim; i++) {
    if (order == 2) {
      if (mp->bonds[i].fitcol >= 0) {
	mp->bonds[i].samprow = j;
	j++;
      }
      else {
	mp->bonds[i].samprow = -1;
      }
    }
    else if (order == 3) {
      if (mp->angls[i].fitcol >= 0) {
	mp->angls[i].samprow = j;
	j++;
      }
      else {
	mp->angls[i].samprow = -1;
      }
    }
    else if (order == 4) {
      if (mp->torsions[i].fitcol >= 0) {
	mp->torsions[i].samprow = j;
	j++;
      }
      else {
	mp->torsions[i].samprow = -1;
      }
    }
  }
}

/***=======================================================================***/
/*** AccumulateSamplingTable: accumulate the sampling histogram based on   ***/
/***                          all conformations.  This routine is          ***/
/***                          encapsulated to allow it to be used in       ***/
/***                          multiple cases.                              ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   S:        the sampling matrix                                       ***/
/***   order:    the order of bonded terms to seek out (2 = bonds,         ***/
/***             3 = angles, 4 = dihedrals)                                ***/
/***   sysid:    system ID number, -1 for accumulation over all systems    ***/
/***   track:    activates tracking of the exact values of the variable;   ***/
/***             only useable if there is only 1 row in S                  ***/
/***   record:   all values of the variable encountered                    ***/
/***=======================================================================***/
static void AccumulateSamplingTable(prmset *mp, imat *S, int order, int sysid,
				    int track, double* record)
{
  int i, j, dev, nrec;
  double bdisp, adisp, hdisp;
  bidx *tbond;
  aidx *tangl;
  hidx *tdihe;

  /*** Compile sampling histograms ***/
  MapSamplingIndices(mp, order);
  nrec = 0;
  for (i = 0; i < mp->nconf; i++) {
    if (sysid >= 0 && mp->conf[i].GroupNum != sysid) {
      continue;
    }
    if (order == 2) {
      for (j = 0; j < mp->conf[i].bmap.nbond; j++) {
        tbond = &mp->conf[i].bmap.id[j];
        if (mp->bonds[tbond->key].fitcol < 0) {
          continue;
        }
        bdisp = mp->conf[i].bmap.val[j];
        dev = (bdisp / 0.027777) + 9.00001;
        if (dev >= 0 && dev < 18) {
          S->map[mp->bonds[tbond->key].samprow][dev] += 1;
        }
	if (track == 1) {
	  record[nrec] = bdisp;
	  nrec++;
	}
      }
    }
    else if (order == 3) {
      for (j = 0; j < mp->conf[i].amap.nangl; j++) {
        tangl = &mp->conf[i].amap.id[j];
        if (mp->angls[tangl->key].fitcol < 0) {
          continue;
        }
        adisp = mp->conf[i].amap.val[j];
        dev = (adisp / 0.023271) + 9.00001;
        if (dev >= 0 && dev < 18) {
          S->map[mp->angls[tangl->key].samprow][dev] += 1;
        }
	if (track == 1) {
	  record[nrec] = adisp;
	  nrec++;
	}
      }
    }
    else {
      for (j = 0; j < mp->conf[i].hmap.ntterm; j++) {
        tdihe = &mp->conf[i].hmap.id[j];
        if (mp->torsions[tdihe->key].fitcol < 0) {
          continue;
        }
        hdisp = mp->conf[i].hmap.val[j];
        dev = hdisp * 9.0 / PI + 9.00001;
        if (dev >= 0 && dev < 18) {
          S->map[mp->torsions[tdihe->key].samprow][dev] += 1;
        }
	if (track == 1) {
	  record[nrec] = hdisp;
	  nrec++;
	}
      }
    }
  }
}

/***=======================================================================***/
/*** CharacterHistogram: print a character-keyed histogram, values ranging ***/
/***                     from 0 to 10.                                     ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   S:      integer matrix containing the histogram data                ***/
/***   ncol:   the number of bins in the histogram                         ***/
/***   outp:   the output file                                             ***/
/***=======================================================================***/
static void CharacterHistogram(int* S, int ncol, FILE *outp)
{
  int i;
  char occmap[16];

  sprintf(occmap, " .:=eoUO0@");
  for (i = 0; i < ncol; i++) {
    if (S[i] < 10) {
      fprintf(outp, "%c", occmap[S[i]]);
    }
    else {
      fprintf(outp, "X");
    }
  }
}

/***=======================================================================***/
/*** MatchSamprow2Column: search through all adjustable variables of a     ***/
/***                      specified order and find the variable's index    ***/
/***                      into the fitting matrix columns based on its     ***/
/***                      index into the sampling matrix.                  ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:      fitting control data (contains notes about the fitting     ***/
/***            matrix)                                                    ***/
/***   sidx:    the index into the sampling table                          ***/
/***   order:   2 = bonds, 3 = angles, 4 = torsions                        ***/
/***=======================================================================***/
static int MatchSamprow2Column(prmset *mp, int sidx, int order)
{
  int i;

  if (order == 2) {
    for (i = 0; i < mp->nbond; i++) {
      if (mp->bonds[i].samprow == sidx && mp->bonds[i].fitcol != -1) {
	return mp->bonds[i].fitcol;
      }
    }
    printf("MatchSamprow2Column >> Unable to find match for bond "
	   "%d, %.4s-%.4s.\n", i, mp->bonds[i].atype, mp->bonds[i].btype);
    exit(1);
  }
  else if (order == 3) {
    for (i = 0; i < mp->nangl; i++) {
      if (mp->angls[i].samprow == sidx && mp->angls[i].fitcol != -1) {
	return mp->angls[i].fitcol;
      }
    }
    printf("MatchSamprow2Column >> Unable to find match for angle "
	   "%d, %.4s-%.4s-%.4s.\n", i, mp->angls[i].atype, mp->angls[i].btype,
	   mp->angls[i].ctype);
    exit(1);
  }
  else if (order == 4) {
    for (i = 0; i < mp->ntor; i++) {
      if (mp->torsions[i].samprow == sidx && mp->torsions[i].fitcol != -1) {
	return mp->torsions[i].fitcol;
      }
    }
    printf("MatchSamprow2Column >> Unable to find match for angle "
	   "%d, %.4s-%.4s-%.4s-%.4s.\n", i, mp->torsions[i].atype,
	   mp->torsions[i].btype, mp->torsions[i].ctype,
	   mp->torsions[i].dtype);
    exit(1);
  }

  return -1;
}

/***=======================================================================***/
/*** FindColumnTerm: find the term, be it a bond, angle, or dihedral,      ***/
/***                 corresponding to the column of interest.  It is       ***/
/***                 terribly inefficient to search through all bonds,     ***/
/***                 angles, and dihedrals just to find the one that feeds ***/
/***                 into a particular matrix column when all of that data ***/
/***                 could be recorded, but this routine is only called at ***/
/***                 output and therefore saves some complexity in the     ***/
/***                 data structures.                                      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:      fitting control data (contains notes about the fitting     ***/
/***            matrix)                                                    ***/
/***   cc:      the column to match a term against                         ***/
/***   outp:    the output file                                            ***/
/***   order:   returned value, identifying the column as associated with  ***/
/***            a bond (return 2), angle (return 3), or torsion (return 4) ***/
/***   rval:    flag to set the actions taken upon finding the term        ***/
/***            corresponding to the column cc                             ***/
/***=======================================================================***/
static double FindColumnTerm(prmset *mp, int cc, FILE *outp, int *order,
                             int rval)
{
  int i;

  *order = -1;

  /*** Leading spaces for output in correlated columns output ***/
  if (rval == 0) {
    fprintf(outp, "  ");
  }
  for (i = 0; i < mp->nbond; i++) {
    if (mp->bonds[i].fitcol == cc) {
      *order = 2;
      if (rval <= 0) {
        fprintf(outp, " BOND %.2s %.2s      ", mp->bonds[i].atype,
                mp->bonds[i].btype);
      }
      else if (rval == 1) {
        return  mp->bonds[i].K;
      }
    }
  }
  for (i = 0; i < mp->nangl; i++) {
    if (mp->angls[i].fitcol == cc) {
      *order = 3;
      if (rval <= 0) {
        fprintf(outp, " ANGL %.2s %.2s %.2s   ", mp->angls[i].atype,
                mp->angls[i].btype, mp->angls[i].ctype);
      }
      else if (rval == 1) {
        return mp->angls[i].K;
      }
    }
  }
  for (i = 0; i < mp->ntor; i++) {
    if (mp->torsions[i].fitcol == cc) {
      if (rval <= 0) {
        *order = 4;
        if (mp->torsions[i].impr == 0) {
          fprintf(outp, " DIHE ");
        }
        else {
          fprintf(outp, " IMPR ");
        }
        fprintf(outp, "%.2s %.2s %.2s %.2s", mp->torsions[i].atype,
                mp->torsions[i].btype, mp->torsions[i].ctype,
                mp->torsions[i].dtype);
      }
      else if (rval == 1) {
        return mp->torsions[i].K;
      }
    }
  }

  /*** Carriage return for output relating to term glossary ***/
  if (rval == -1) {
    fprintf(outp, "\n");
  }

  /*** By default return zero ***/
  return 0.0;
}

/***=======================================================================***/
/*** SamplingMatrixOutput: print the output of the sampling matrix, in     ***/
/***                       human-readable format.                          ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   order:    the order of bonded terms to seek out (2 = bonds,         ***/
/***             3 = angles, 4 = dihedrals)                                ***/
/***   outp:     the output file                                           ***/
/***   S:        the sampling matrix                                       ***/
/***   x:        the solution vector                                       ***/
/***=======================================================================***/
static void SamplingMatrixOutput(prmset *mp, int order, FILE *outp, imat *S,
                                 double* x)
{
  int i, j, nvar, ninst, ifitcol;
  double Korig;

  if (order == 2) {
    fprintf(outp, "  Term  Amber Types  Count -0.25   +0    0.25"
	    "   Fitted  Original\n");
  }
  else if (order == 3) {
    fprintf(outp, "  Term  Amber Types  Count -15     +0      15"
	    "   Fitted  Original\n");
  }
  else if (order == 4) {
    fprintf(outp, "  Term  Amber Types  Count -PI     +0      PI"
	    "   Fitted  Original\n");
  }
  fprintf(outp, "  ----  -----------  ----- "
          "------------------  -------- --------\n");
  nvar = 0;
  for (i = 0; i < S->row; i++) {
    if (order == 2) {
      while (mp->bonds[nvar].fitcol < 0) {
        nvar++;
      }
      ninst = mp->bonds[nvar].ninst;
      fprintf(outp, "  BOND  %.2s %.2s       ", mp->bonds[nvar].atype,
              mp->bonds[nvar].btype);
    }
    if (order == 3) {
      while (mp->angls[nvar].fitcol < 0) {
        nvar++;
      }
      ninst = mp->angls[nvar].ninst;
      fprintf(outp, "  ANGL  %.2s %.2s %.2s    ", mp->angls[nvar].atype,
              mp->angls[nvar].btype, mp->angls[nvar].ctype);
    }
    if (order == 4) {
      while (mp->torsions[nvar].fitcol < 0) {
        nvar++;
      }
      ninst = mp->torsions[nvar].ninst;
      if (mp->torsions[nvar].impr == 0) {
        fprintf(outp, "  DIHE  ");
      }
      else {
        fprintf(outp, "  IMPR  ");
      }
      fprintf(outp, "%.2s %.2s %.2s %.2s ", mp->torsions[nvar].atype,
              mp->torsions[nvar].btype, mp->torsions[nvar].ctype,
              mp->torsions[nvar].dtype);
    }
    nvar++;
    fprintf(outp, "%6d ", ninst);
    CharacterHistogram(S->map[i], S->col, outp);
    ifitcol = MatchSamprow2Column(mp, i, order);
    Korig = FindColumnTerm(mp, ifitcol, outp, &j, 1);
    fprintf(outp, "  %8.2lf %8.2lf\n", x[ifitcol], Korig);
  }
  fprintf(outp, "\n");
}

/***=======================================================================***/
/*** ParameterWarnings: runs a series of diagnostics on the parameter and  ***/
/***                    its possible effects on the resulting model.       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   n:        the column of the parameter to detail                     ***/
/***   outp:     the output file                                           ***/
/***=======================================================================***/
static void ParameterWarnings(prmset *mp, int n, FILE *outp)
{
  int order;
  double Korig;

  /*** Get the original stiffness and order of this parameter ***/
  Korig = FindColumnTerm(mp, n, outp, &order, 1);

  /*** This does something with Korig, keeps ***/
  /*** the compiler from complaining         ***/
  if (abs(Korig) > -1.0e5) {
    fprintf(outp, " ");
  }

  /*** Make a complete list of every time this parameter is ***/
  /*** encountered in the context of a particular system    ***/

  /*** Carriage return to end the line in the output file ***/
  fprintf(outp, "\n");
}

/***=======================================================================***/
/*** DetailParameter: details the instances in which a parameter occurs in ***/
/***                  the fitting set.                                     ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   n:        the column of the parameter to detail                     ***/
/***   outp:     the output file                                           ***/
/***=======================================================================***/
static void DetailParameter(prmset *mp, int n, FILE *outp)
{
  int h, i, j, bkey, akey, hkey, ninst, atma, atmb, atmc, atmd;
  int resa, resb, resc, resd, nspc, rescon, bcol, acol, hcol;
  double* dtmp;
  char *tpatoms, *tpres;
  imat* BondSamp;
  imat* AnglSamp;
  imat* DiheSamp;
  itrack* inst;

  /*** Determine the maximum possible number of instances ***/
  j = 0;
  for (i = 0; i < mp->nconf; i++) {
    j += mp->conf[i].bmap.nbond;
    j += mp->conf[i].amap.nangl;
    j += mp->conf[i].hmap.ntterm;
  }
  inst = (itrack*)malloc(j*sizeof(itrack));

  /*** Check for this term in exactly one example of each system ***/
  ninst = 0;
  for (h = 0; h < mp->nunisys; h++) {
    i = mp->FirstConf[h];
    tpatoms = mp->conf[i].tp->AtomNames;
    tpres = mp->conf[i].tp->ResNames;
    for (j = 0; j < mp->conf[i].bmap.nbond; j++) {
      bkey = mp->conf[i].bmap.id[j].key;
      if (mp->bonds[bkey].fitcol == n) {

        /*** Column n pertains to a bond in this molecule ***/
        inst[ninst].sysid = mp->conf[i].GroupNum;
        inst[ninst].sysno = i;
        inst[ninst].order = 2;
        inst[ninst].tnum = j;
        atma = mp->conf[i].bmap.id[j].a;
        atmb = mp->conf[i].bmap.id[j].b;
	resa = LocateResID(mp->conf[i].tp, atma, 0, mp->conf[i].tp->nres);
	resb = LocateResID(mp->conf[i].tp, atmb, 0, mp->conf[i].tp->nres);
        sprintf(inst[ninst].atom, "%.4s%.4s", &tpatoms[4*atma],
                &tpatoms[4*atmb]);
        sprintf(inst[ninst].res, "%.4s%.4s", &tpres[4*resa], &tpres[4*resb]);
        ninst++;
      }
    }
    for (j = 0; j < mp->conf[i].amap.nangl; j++) {
      akey = mp->conf[i].amap.id[j].key;
      if (mp->angls[akey].fitcol == n) {

        /*** Column n pertains to an angle in this molecule ***/
        inst[ninst].sysid = mp->conf[i].GroupNum;
        inst[ninst].sysno = i;
        inst[ninst].order = 3;
        inst[ninst].tnum = j;
        atma = mp->conf[i].amap.id[j].a;
        atmb = mp->conf[i].amap.id[j].b;
        atmc = mp->conf[i].amap.id[j].c;
	resa = LocateResID(mp->conf[i].tp, atma, 0, mp->conf[i].tp->nres);
	resb = LocateResID(mp->conf[i].tp, atmb, 0, mp->conf[i].tp->nres);
	resc = LocateResID(mp->conf[i].tp, atmc, 0, mp->conf[i].tp->nres);
        sprintf(inst[ninst].atom, "%.4s%.4s%.4s", &tpatoms[4*atma],
                &tpatoms[4*atmb], &tpatoms[4*atmc]);
        sprintf(inst[ninst].res, "%.4s%.4s%.4s", &tpres[4*resa],
                &tpres[4*resb], &tpres[4*resc]);
        ninst++;
      }
    }
    for (j = 0; j < mp->conf[i].hmap.ntterm; j++) {
      hkey = mp->conf[i].hmap.id[j].key;
      if (mp->torsions[hkey].fitcol == n) {

        /*** Column n pertains to a dihedral in this molecule ***/
        inst[ninst].sysid = mp->conf[i].GroupNum;
        inst[ninst].sysno = i;
        inst[ninst].order = 4;
        inst[ninst].tnum = j;
        atma = mp->conf[i].hmap.id[j].a;
        atmb = mp->conf[i].hmap.id[j].b;
        atmc = mp->conf[i].hmap.id[j].c;
        atmd = mp->conf[i].hmap.id[j].d;
        resa = LocateResID(mp->conf[i].tp, atma, 0, mp->conf[i].tp->nres);
        resb = LocateResID(mp->conf[i].tp, atmb, 0, mp->conf[i].tp->nres);
        resc = LocateResID(mp->conf[i].tp, atmc, 0, mp->conf[i].tp->nres);
        resd = LocateResID(mp->conf[i].tp, atmd, 0, mp->conf[i].tp->nres);
        sprintf(inst[ninst].atom, "%.4s%.4s%.4s%.4s", &tpatoms[4*atma],
                &tpatoms[4*atmb], &tpatoms[4*atmc], &tpatoms[4*atmd]);
        sprintf(inst[ninst].res, "%.4s%.4s%.4s%.4s", &tpres[4*resa],
                &tpres[4*resb], &tpres[4*resc], &tpres[4*resd]);
        ninst++;
      }
    }
  }

  /*** Create sampling tables for conformations  ***/
  /*** sampled by each individual system.  dtmp  ***/
  /*** is set just to avoid compiler complaints. ***/
  if (mp->nbvar > 0) BondSamp = (imat*)malloc(mp->nunisys*sizeof(imat));
  if (mp->navar > 0) AnglSamp = (imat*)malloc(mp->nunisys*sizeof(imat));
  if (mp->nhvar > 0) DiheSamp = (imat*)malloc(mp->nunisys*sizeof(imat));
  dtmp = mp->corrval;
  for (i = 0; i < mp->nunisys; i++) {
    if (mp->nbvar > 0) {
      BondSamp[i] = CreateImat(mp->nbvar, 18);
      AccumulateSamplingTable(mp, &BondSamp[i], 2, i, 0, dtmp);
    }
    if (mp->navar > 0) {
      AnglSamp[i] = CreateImat(mp->navar, 18);
      AccumulateSamplingTable(mp, &AnglSamp[i], 3, i, 0, dtmp);
    }
    if (mp->nhvar > 0) {
      DiheSamp[i] = CreateImat(mp->nhvar, 18);
      AccumulateSamplingTable(mp, &DiheSamp[i], 4, i, 0, dtmp);
    }
  }

  /*** Recount the instances ***/
  for (i = 0; i < ninst; i++) {
    fprintf(outp, "  ");
    fprintf(outp, "%.4s(", inst[i].res);
    rescon = 0;
    nspc = 7;
    for (j = 0; j < inst[i].order; j++) {
      if (str4cmp(&inst[i].res[4*j], &inst[i].res[4*rescon]) == 0) {
        fprintf(outp, "%.4s", &inst[i].atom[4*j]);
        nspc += 4;
      }
      else {
        fprintf(outp, ") %.4s(%.4s", &inst[i].res[4*j], &inst[i].atom[4*j]);
        rescon = j;
        nspc += 11;
      }
      if (j < inst[i].order-1 &&
          str4cmp(&inst[i].res[4*(j+1)], &inst[i].res[4*rescon]) == 0) {
        fprintf(outp, " ");
        nspc += 1;
      }
    }
    fprintf(outp, ")");
    nspc += 1;
    for (j = nspc; j < 41; j++) {
      fprintf(outp, " ");
    }
    fprintf(outp, "%5d  ", mp->GroupCount[inst[i].sysid]);
    if (inst[i].order == 2) {
      bkey = mp->conf[inst[i].sysno].bmap.id[inst[i].tnum].key;
      bcol = mp->bonds[bkey].samprow;
      CharacterHistogram(BondSamp[inst[i].sysid].map[bcol], 18, outp);
    }
    else if (inst[i].order == 3) {
      akey = mp->conf[inst[i].sysno].amap.id[inst[i].tnum].key;
      acol = mp->angls[akey].samprow;
      CharacterHistogram(AnglSamp[inst[i].sysid].map[acol], 18, outp);
    }
    else if (inst[i].order == 4) {
      hkey = mp->conf[inst[i].sysno].hmap.id[inst[i].tnum].key;
      hcol = mp->torsions[hkey].samprow;
      CharacterHistogram(DiheSamp[inst[i].sysid].map[hcol], 18, outp);
    }
    ParameterWarnings(mp, n, outp);
  }

  /*** Free allocated memory ***/
  for (i = 0; i < mp->nunisys; i++) {
    if (mp->nbvar > 0) DestroyImat(&BondSamp[i]);
    if (mp->navar > 0) DestroyImat(&AnglSamp[i]);
    if (mp->nhvar > 0) DestroyImat(&DiheSamp[i]);
  }
  if (mp->nbvar > 0) free(BondSamp);
  if (mp->navar > 0) free(AnglSamp);
  if (mp->nhvar > 0) free(DiheSamp);
  free(inst);
}

/***=======================================================================***/
/*** PrintSamplingTable: print a table to describe the sampling of each    ***/
/***                     fitted parameter across all conformations.        ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   order:    the order of bonded terms to seek out (2 = bonds,         ***/
/***             3 = angles, 4 = dihedrals)                                ***/
/***   x:        the solution vector (to print alongside sampling data)    ***/
/***   outp:     the output file                                           ***/
/***=======================================================================***/
static void PrintSamplingTable(prmset *mp, int order, double* x, FILE *outp)
{
  imat samp;

  /*** Allocate the histogram matrix ***/
  samp = (order == 2) ? CreateImat(mp->nbvar, 18) :
    (order == 3) ? CreateImat(mp->navar, 18) : CreateImat(mp->nhvar, 18);

  /*** Accumulate the histogram ***/
  AccumulateSamplingTable(mp, &samp, order, -1, 0, mp->corrval);

  /*** Print the histogram ***/
  if (order == 2 && mp->nbvar > 0) {
    fprintf(outp, " Bond sampling:\n"
	    "                            Bins in Angstroms\n");
    SamplingMatrixOutput(mp, 2, outp, &samp, x);
  }
  if (order == 3 && mp->navar > 0) {
    fprintf(outp, " Angle sampling:\n"
	    "                             Bins in Degrees \n");
    SamplingMatrixOutput(mp, 3, outp, &samp, x);
  }
  if (order == 4 && mp->nhvar > 0) {
    fprintf(outp, " Torsion term sampling:\n"
            "                             Bins in Radians \n");
    SamplingMatrixOutput(mp, 4, outp, &samp, x);
  }

  /*** Free allocated memory ***/
  DestroyImat(&samp);
}

/***=======================================================================***/
/*** PrintEAConstants: print the energy adjustment constants emerging from ***/
/***                   the fit.  Large values of these constants can       ***/
/***                   indicate that certain fitted parameters are being   ***/
/***                   allowed to contribute large baseline values to the  ***/
/***                   energy estimates.                                   ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   x:        the solution vector (to print alongside sampling data)    ***/
/***   outp:     the output file                                           ***/
/***=======================================================================***/
static void PrintEAConstants(prmset *mp, double* x, FILE *outp)
{
  int h, i, j, k, vstart;
  double maxwt;

  if (mp->nunisys == 1) {
    if (mp->wtfloor > 100.0) {
      fprintf(outp,
	      "      System      Adjustment\n"
	      " ---------------- ----------");
    }
    else {
      fprintf(outp,
	      "      System      Adjustment  Weight\n"
	      " ---------------- ---------- -------");
    }
  }
  else {
    if (mp->wtfloor > 100.0) {
      fprintf(outp,
	      "      System      Adjustment         System      "
	      "Adjustment\n"
	      " ---------------- ----------    ---------------- "
	      "----------");
    }
    else {
      fprintf(outp,
	      "      System      Adjustment  Weight         System      "
	      "Adjustment  Weight\n"
	      " ---------------- ---------- -------    ---------------- "
	      "---------- -------");
    }
  }
  h = 2;
  vstart = mp->nparm;
  for (i = 0; i < mp->nunisys; i++) {
    if (h == 2) {
      h = 0;
      fprintf(outp, "\n ");
    }
    for (j = 0; j < mp->nconf; j++) {
      if (mp->conf[j].GroupNum == i) {
	if (mp->wtfloor > 100.0) {
	  fprintf(outp, "%-16.16s %10.4lf    ", mp->conf[j].tp->source,
		  x[vstart]);
	}
	else {
	  maxwt = 0.0;
	  for (k = 0; k < mp->nconf; k++) {
	    if (mp->conf[k].GroupNum == i && mp->conf[k].wt > maxwt) {
	      maxwt = mp->conf[k].wt;
	    }
	  }
	  fprintf(outp, "%-16.16s %10.4lf %7.4lf    ", mp->conf[j].tp->source,
		  x[vstart], maxwt);
	}
	h++;
	break;
      }
    }
    vstart++;
  }
  if (h != 0) {
    fprintf(outp, "\n");
  }
  fprintf(outp, "\n");
  if (mp->fitscee == 1 || mp->fitscnb == 1) {
    PrintVADesc(0, " ", 1, " ", 1, "In addition, 1:4 scaling factors were "
		"also fitted, and should be applied in any simulations with "
		"the resulting parameters.\n", 77, 0, outp);
    if (mp->fitscee == 1) {
      fprintf(outp, " - SCEE = %10.6lf\n",
	      1.0/x[mp->nbvar+mp->navar+mp->nhvar]);
    }
    if (mp->fitscnb == 1) {
      fprintf(outp, " - SCNB = %10.6lf\n\n",
	      1.0/x[mp->nbvar+mp->navar+mp->nhvar+mp->fitscee]);
    }
  }
}

/***=======================================================================***/
/*** PrintMatrixAnalysis: print an analysis of the fitting matrix,         ***/
/***                      identifying highly correlated columns as well as ***/
/***                      columns with no fitting data in them.            ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:      fitting control data (contains notes about the fitting     ***/
/***            matrix)                                                    ***/
/***   x:       the solution vector                                        ***/
/***   outp:    the output file                                            ***/
/***=======================================================================***/
static void PrintMatrixAnalysis(prmset *mp, double* x, FILE *outp)
{
  int i, order;

  /*** Zero-data columns? ***/
  if (mp->nzerocol == 0) {
    fprintf(outp, " - No columns with zero data detected\n");
  }
  else {
    fprintf(outp, " - %d columns with zero data detected, corresponding to "
            "parameters:\n", mp->nzerocol);
    for (i = 0; i < mp->nzerocol; i++) {
      FindColumnTerm(mp, mp->zerocol[i], outp, &order, 0);
      fprintf(outp, "\n");
    }
  }

  /*** Correlated columns? ***/
  if (mp->ncorrcol == 0) {
    fprintf(outp, " - No correlated columns detected\n");
  }
  else {
    fprintf(outp, " - %d highly correlated column pairs detected, "
            "corresponding to parameters:\n\n", mp->ncorrcol);
    fprintf(outp, "     First term,        Second term,     Fitted   Fitted   "
            "   Pearson\n");
    fprintf(outp, "   Type    Atoms      Type    Atoms      value 1  value 2  "
            " Correlation\n");
    fprintf(outp, "   ---- -----------   ---- -----------   -------- -------- "
            " -----------\n");
    for (i = 0; i < mp->ncorrcol; i++) {
      FindColumnTerm(mp, mp->corrcol[2*i], outp, &order, 0);
      FindColumnTerm(mp, mp->corrcol[2*i+1], outp, &order, 0);
      fprintf(outp, "   %8.2lf %8.2lf  %11.8lf\n",
              x[mp->corrcol[2*i]], x[mp->corrcol[2*i+1]], mp->corrval[i]);
    }
  }
}

/***=======================================================================***/
/*** ParameterDescriptions: this routine will print information to put all ***/
/***                        fitted parameters in context.  It details the  ***/
/***                        residues and atoms that each parameter affects ***/
/***                        and gives counts on the number of instances of ***/
/***                        each occurrence.  This routine will also check ***/
/***                        for criteria that may indicate a parameter has ***/
/***                        been poorly fitted.  Output is tabulated in a  ***/
/***                        special section of the mdout file.             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   outp:     pointer to the mdout file                                 ***/
/***=======================================================================***/
static void ParameterDescription(prmset *mp, FILE *outp)
{
  int i, order;

  fprintf(outp, " [ Parameter ]                                  Sampling "
          "Histogram\n"
          " Residue (Atom Names)                    Count  MIN     +0     MAX "
          " Warnings\n"
          " --------------------------------------- -----  ------------------ "
          " --------\n");
  for (i = 0; i < mp->nparm; i++) {
    if (mp->verbose == 1) {
      fprintf(stderr, "\rmdgx >> Detailing parameter %4d of %4d.", i,
	      mp->nparm);
      fflush(stderr);
    }
    FindColumnTerm(mp, i, outp, &order, -1);
    DetailParameter(mp, i, outp);
  }
  if (mp->verbose == 1) {
    printf("\rmdgx >> Descriptions complete for all instances of %4d "
	   "parameters.\n", mp->nparm);
  }
}

/***=======================================================================***/
/*** CheckMinima: this function computes the minima of each fourier series ***/
/***              and determines whether any of them are significantly     ***/
/***              lower than the points at which the fourier series is     ***/
/***              sampled in the training data.                            ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   tj:       trajectory control information (for the input file name)  ***/
/***   mdout:    pointer to the mdout file                                 ***/
/***   x:        the solution vector from the fit                          ***/
/***=======================================================================***/
static void CheckMinima(prmset *mp, trajcon *tj, double* x)
{
  int i, j, ntt, nser, totsamp;
  int* examined;
  int* sterm;
  double hangle, sangle, tfunc;
  FILE *outp;
  time_t ct;
  torterm *tortmp;

  /*** Checks are done here to avoid cluttering ParamReport() ***/
  if (mp->sr[0] == '\0' || mp->nhvar == 0) {
    return;
  }

  /*** Array to flag torsions that have ***/
  /*** been included in fourier series  ***/
  examined = (int*)calloc(mp->ntor, sizeof(int));

  /*** Array to store the torsion terms in a series ***/
  sterm = (int*)malloc(mp->ntor*sizeof(int));

  /*** Prepare a file to plot each series ***/
  outp = fopen(mp->sr, "w");
  ct = time(&ct);
  fprintf(outp, "%%\n%% Printed by mdgx on %s", asctime(localtime(&ct)));
  fprintf(outp, "%% This MatLab-readable results file plots torsion "
	  "potentials fitted\n%% according to input found in %s.\n%%\n\n",
	  tj->inpname);

  /*** Loop over all fitted torsion parameters ***/
  nser = 0;
  for (i = 0; i < mp->ntor; i++) {
    if (examined[i] == 1 || mp->torsions[i].fitcol < 0) {
      continue;
    }

    /*** Mark this torsion term ***/
    examined[i] = 1;    
    sterm[0] = i;
    ntt = 1;
    totsamp = mp->torsions[i].ninst;

    /*** Determine whether this torsion is part of a series ***/
    for (j = i+1; j < mp->ntor; j++) {
      if (TypeCompare(mp->torsions[i].atype, mp->torsions[i].btype,
		      mp->torsions[i].ctype, mp->torsions[i].dtype,
		      mp->torsions[j].atype, mp->torsions[j].btype,
		      mp->torsions[j].ctype, mp->torsions[j].dtype,
		      4, 1) == 2) {
	examined[j] = 1;
	sterm[ntt] = j;
	ntt++;
	totsamp += mp->torsions[j].ninst;
      }
    }

    /*** Print the output file ***/
    nser++;
    fprintf(outp, "%% Series for atom types %.4s %.4s %.4s %.4s\n"
	    "%% %d samples, %d Fourier terms\n",
	    mp->torsions[i].atype, mp->torsions[i].btype,
	    mp->torsions[i].ctype, mp->torsions[i].dtype, totsamp, ntt);
    fprintf(outp, "FSeries%d = [\n", nser);
    for (hangle = 0.0; hangle < 2.0*PI; hangle += 0.001) {
      tfunc = 0.0;
      for (j = 0; j < ntt; j++) {
	tortmp = &mp->torsions[sterm[j]];
	sangle = tortmp->pn*hangle - tortmp->phase;
	tfunc += x[tortmp->fitcol] * (1.0 + cos(sangle));
      }
      fprintf(outp, "%9.4lf\n", tfunc);
    }
    fprintf(outp, "];\n");
  }

  /*** Close output file ***/
  fclose(outp);

  /*** Free allocated memory ***/
  free(examined);
  free(sterm);
}

/***=======================================================================***/
/*** CountFittedTerms: count the number of fitted terms in a system.       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   sysno:    the system number                                         ***/
/***   mp:       the master parameter and fitting set                      ***/
/***=======================================================================***/
static int CountFittedTerms(int sysno, prmset *mp)
{
  int i, exid;
  int* pterm;
  mmsys *exconf;

  /*** First, find one example of this system ***/
  exid = -1;
  for (i = 0; i < mp->nconf; i++) {
    if (mp->conf[i].GroupNum == sysno) {
      exid = i;
      break;
    }
  }
  if (exid == -1) {
    printf("CountFittedTerms >> Error.  System %d does not exist.\n", sysno);
    exit(1);
  }
  exconf = &mp->conf[exid];

  /*** Bonds, angles, and dihedrals ***/
  pterm = (int*)calloc(mp->nparm, sizeof(int));
  for (i = 0; i < exconf->bmap.nbond; i++) {
    if (mp->bonds[exconf->bmap.id[i].key].fitcol >= 0) {
      pterm[mp->bonds[exconf->bmap.id[i].key].fitcol] = 1;
    }
  }
  for (i = 0; i < exconf->amap.nangl; i++) {
    if (mp->angls[exconf->amap.id[i].key].fitcol >= 0) {
      pterm[mp->angls[exconf->amap.id[i].key].fitcol] = 1;
    }
  }
  for (i = 0; i < exconf->hmap.ntterm; i++) {
    if (mp->torsions[exconf->hmap.id[i].key].fitcol >= 0) {
      pterm[mp->torsions[exconf->hmap.id[i].key].fitcol] = 1;
    }
  }

  return ISum(pterm, mp->nparm);
}

/***=======================================================================***/
/*** PrintAccuracy: report the accuracy of this model in relation to the   ***/
/***                training set.                                          ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   A:        the fitting matrix, in its original form before the QR    ***/
/***             decomposition                                             ***/
/***   x:        the solution vector, containing all fitted parameters     ***/
/***   tj:       trajectory control information                            ***/
/***   outp:     the output file                                           ***/
/***=======================================================================***/
void PrintAccuracy(prmset *mp, dmat *A, double* x, trajcon *tj, FILE *outp)
{
  int i, j, k, m, ncfg, fsys, natom;
  int* nsysvar;
  int* confid;
  double fitnrg, eQQ, eLJ;
  double *atmp;
  double* enorm;
  double* eorig;
  double* efin;
  double* eelec;
  double* elj;
  double* ebond;
  double* eangl;
  double* edihe;
  double* efit;
  FILE *mscript;
  time_t ct;

  /*** If a MatLab readable display is requested, print it ***/
  if (mp->ao[0] != '\0') {
    mscript = FOpenSafe(mp->ao, tj->OverwriteOutput);
    ct = time(&ct);
    fprintf(mscript, "%% Model accuracy report for fitting ordered by\n"
	    "%% %s.\n", tj->inpname);
    fprintf(mscript, "%% Written on %s\n", asctime(localtime(&ct)));
  }
  fprintf(outp, 
	  "                      RMS Error         Correlation    Fitted   "
	  "Fitted   Model\n"
          "      System        Orig    Fitted    Orig    Fitted   Terms    "
	  "Energy   Diff.\n"
          " ----------------  -------  -------  -------  -------  ------  "
	  "-------  -------\n");

  /*** Count the number of fitted terms in each system ***/
  nsysvar = (int*)malloc(mp->nunisys*sizeof(int));
  for (i = 0; i < mp->nunisys; i++) {
    nsysvar[i] = CountFittedTerms(i, mp);
  }

  /*** Compute the accuracy of the model ***/
  confid = (int*)malloc(mp->nconf*sizeof(int));
  enorm = (double*)malloc(mp->nconf*sizeof(double));
  eorig = (double*)malloc(mp->nconf*sizeof(double));
  efin = (double*)malloc(mp->nconf*sizeof(double));
  eelec = (double*)malloc(mp->nconf*sizeof(double));
  elj = (double*)malloc(mp->nconf*sizeof(double));
  ebond = (double*)malloc(mp->nconf*sizeof(double));
  eangl = (double*)malloc(mp->nconf*sizeof(double));
  edihe = (double*)malloc(mp->nconf*sizeof(double));
  efit = (double*)calloc(mp->nconf, sizeof(double));
  for (i = 0; i < mp->nunisys; i++) {
    ncfg = 0;
    for (j = 0; j < mp->nconf; j++) {
      if (mp->conf[j].GroupNum != i) {
        continue;
      }
      if (mp->ao[0] != '\0' && ncfg == 0) {
	fprintf(mscript, "%% System %s\ndata%d = [\n"
		"%% Target  Original    Model     Error\n",
		mp->conf[j].tp->source, i);
      }
      fsys = j;
      fitnrg = 0.0;
      atmp = A->map[j];
      for (k = 0; k < A->col; k++) {
        fitnrg += atmp[k] * x[k];
      }
      confid[ncfg] = j;
      ebond[ncfg] = DSum(mp->conf[j].bmap.Ucontrib, mp->conf[j].bmap.nbond);
      eangl[ncfg] = DSum(mp->conf[j].amap.Ucontrib, mp->conf[j].amap.nangl);
      edihe[ncfg] = DSum(mp->conf[j].hmap.Ucontrib, mp->conf[j].hmap.ntterm);
      natom = mp->conf[j].tp->natom;
      eLJ = 0.0;
      eQQ = 0.0;
      for (k = 0; k < natom-1; k++) {
	for (m = k+1; m < natom; m++) {
	  eQQ += mp->conf[j].nbnrg.map[k][m] * mp->conf[j].excl.map[k][m];
	  eLJ += mp->conf[j].nbnrg.map[m][k] * mp->conf[j].excl.map[m][k];
	}
      }
      eelec[ncfg] = eQQ;
      elj[ncfg] = eLJ;
      efin[ncfg] = mp->conf[j].nonfitmm + fitnrg;
      enorm[ncfg] = mp->conf[j].enorm;
      eorig[ncfg] = mp->conf[j].eorig;
      efit[i] += fitnrg;
      ncfg++;
    }
    efit[i] /= ncfg;
    if (mp->ao[0] != '\0') {
      for (j = 0; j < ncfg; j++) {
	fprintf(mscript, "%9.4lf %9.4lf %9.4lf %9.4lf   %% %s\n", enorm[j],
		eorig[j], efin[j], enorm[j]-efin[j],
		mp->conf[confid[j]].crdsrc);
      }
      fprintf(mscript, "];\n\n");
      fprintf(mscript, "%% System %s\n", mp->conf[confid[0]].tp->source);
      fprintf(mscript, "components%d = [\n", i);
      fprintf(mscript, "%%    elec      LJ        bond      angl      dihe\n");
      for (j = 0; j < ncfg; j++) {
	fprintf(mscript, "%9.4lf %9.4lf %9.4lf %9.4lf %9.4lf    %% %s\n",
		eelec[j], elj[j], ebond[j], eangl[j], edihe[j],
		mp->conf[confid[j]].crdsrc);
      }
      fprintf(mscript, "];\nfigure; hold on;\n"
	      "plot(data%d(:,1), data%d(:,2), 'k.');\n"
	      "plot(data%d(:,1), data%d(:,3), 'r.');\n", i, i, i, i);
      fprintf(mscript, "legend('Original', 'Fitted');\n"
	      "xlabel('Target energy');\nylabel('Model energy');\n");
      fprintf(mscript, "title('Model results for %s');\n",
	      mp->conf[fsys].tp->source);
      fprintf(mscript, "n = input('Ready?');\n\n");
    }
    fprintf(outp, " %-16.16s  %7.4lf  %7.4lf  %7.4lf  %7.4lf  %6d  %7.2lf  "
	    "%7.3lf\n", mp->conf[fsys].tp->source, VecRMSD(eorig, enorm, ncfg),
            VecRMSD(efin, enorm, ncfg), Pearson(eorig, enorm, ncfg),
            Pearson(efin, enorm, ncfg), nsysvar[i], efit[i],
	    VecRMSD(efin, eorig, ncfg));
  }

  /*** Complete the MatLab output if requested ***/
  if (mp->ao[0] != '\0') {  
    fclose(mscript);
  }

  /*** Free allocated memory ***/
  free(enorm);
  free(eorig);
  free(efin);
  free(ebond);
  free(eangl);
  free(edihe);
  free(efit);
  free(elj);
  free(eelec);
  free(nsysvar);
}

/***=======================================================================***/
/*** SortEpacket: function for quicksort of energy packet contributions.   ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***=======================================================================***/
static int SortEpacket(const void *pcA, const void *pcB)
{
  double eA = fabs(((epacket*)pcA)[0].eave);
  double eB = fabs(((epacket*)pcB)[0].eave);

  if (eA < eB) {
    return 1;
  }
  else if (eA > eB) {
    return -1;
  }
  else {
    return 0;
  }
}

/***=======================================================================***/
/*** PrintEnergyContrib: print the energetic contributions of each of the  ***/
/***                     fitted terms, system-by-system.                   ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   A:        the fitting matrix, in its original form before the QR    ***/
/***             decomposition                                             ***/
/***   x:        the solution vector, containing all fitted parameters     ***/
/***   outp:     the output file                                           ***/
/***=======================================================================***/
static void PrintEnergyContrib(prmset *mp, dmat *A, double* x, FILE *outp)
{
  int i, j, k, order;
  double edr, invN;
  double *dtmp;
  epacket* epc;

  /*** Array to store the values of each parameter's contribution ***/
  epc = (epacket*)malloc(mp->nparm*sizeof(epacket));

  /*** Loop over all systems ***/
  for (i = 0; i < mp->nunisys; i++) {

    /*** Print subsection heading ***/
    fprintf(outp, " System: %s\n", mp->tpencyc[i].source);
    fprintf(outp, "                      Fitted      Orig\n");

    /*** Sum up the energies ***/
    for (j = 0; j < mp->nparm; j++) {
      epc[j].fitcol = j;
      epc[j].eave = 0.0;
      epc[j].estd = 0.0;
    }
    for (j = 0; j < mp->nconf; j++) {
      if (mp->conf[j].GroupNum == i) {
	dtmp = A->map[j];
	for (k = 0; k < mp->nparm; k++) {
	  edr = dtmp[k]*x[k];
	  epc[k].eave += edr;
	  epc[k].estd += edr*edr;
	}
      }
    }

    /*** Compute average and standard deviation, then print ***/
    invN = 1.0/mp->GroupCount[i];
    for (j = 0; j < mp->nparm; j++) {
      epc[j].eave *= invN;
      epc[j].estd = sqrt(invN*epc[j].estd - epc[j].eave*epc[j].eave);
    }
    qsort(epc, mp->nparm, sizeof(epacket), SortEpacket);
    for (j = 0; j < mp->nparm; j++) {
      if (fabs(epc[j].eave) > 1.0e-8) {
	FindColumnTerm(mp, epc[j].fitcol, outp, &order, 0);
	if (order >= 0) {
	  fprintf(outp, "  %9.4lf  %9.4lf\n", epc[j].eave, epc[j].estd);
	}
      }
    }
  }
  fprintf(outp, "\n\n");

  /*** Free allocated memory ***/
  free(epc);
}

/***=======================================================================***/
/*** PrintChangeLog: print the atom type change log.                       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   outp:     the output file                                           ***/
/***=======================================================================***/
static void PrintChangeLog(prmset *mp, FILE *outp)
{
  int i;

  fprintf(outp,
	  "      System       Atom No Res    Change   Oper'n  "
	  "Instances (ambmask format)\n"
	  " ----------------  ---- -- ----  --------  ------  "
	  "---------------------------\n");
  for (i = 0; i < mp->nchng; i++) {
    fprintf(outp, "%s", mp->ChangeLog.map[i]);
  }
}

/***=======================================================================***/
/*** ParamReport: report the best parameters and their fit to the target   ***/
/***              data.  New parameters are reported in a frcmod-format    ***/
/***              file, while statistics from the run are reported in a    ***/
/***              text file.                                               ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   A:        the fitting matrix, in its original form before the QR    ***/
/***             decomposition                                             ***/
/***   x:        the solution vector, containing all fitted parameters     ***/
/***   tj:       trajectory control information                            ***/
/***=======================================================================***/
void ParamReport(prmset *mp, dmat *A, double* x, trajcon *tj)
{
  int i, j;
  FILE *outp;
  time_t ct;

  /*** Input file header ***/
  ct = time(&ct);
  outp = FOpenSafe(tj->outbase, tj->OverwriteOutput);
  PrintSplash(outp);
  fprintf(outp, "Run on %s", asctime(localtime(&ct)));

  /*** Reprint the input file ***/
  HorizontalRule(outp, 0);
  fprintf(outp, "\nINPUT LINE TEXT:\n\n");
  PrintParagraph(tj->inpline, 79, outp);
  fprintf(outp, "\nINPUT FILE TEXT:\n\n");
  for (i = 0; i < tj->inptext.row; i++) {
    fprintf(outp, "%s", tj->inptext.map[i]);
  }
  HorizontalRule(outp, 1);

  /*** Print the energies according to the derived model ***/
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(1.)", 4, " ", 1, "Energies of each system, according to "
              "the fitted parameters.  Units on errors are kcal/mol.  The "
              "original model's energies are compared to target energies "
              "after applying an offset to make the average of each set of "
              "energies for any given system equal (see Section 3).\n", 74, 0,
	      outp);
  PrintAccuracy(mp, A, x, tj, outp);
  HorizontalRule(outp, 1);

  /*** Print the sampling in each fitted ***/
  /*** parameter, across all systems     ***/
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(2.)", 4, " ", 1, "Parameter sampling.  Bonds and angles "
              "are binned over 18 intervals based on deviations, positive or "
              "negative, from the ideal bond length.  Each interval "
              "represents 0.028 Angstrom or 0.023 radian deviation for bonds "
              "or angles, respectively.  Dihedrals and impropers are binned "
              "at 10-degree intervals, ranging from zero to 180 (the range of "
              "the arccosine function).\n", 74, 0, outp);
  fprintf(outp, " Sampling histogram key:  (zero)   . : = e o U O 0 @ X "
          "( > ten)\n\n");
  PrintSamplingTable(mp, 2, x, outp);
  PrintSamplingTable(mp, 3, x, outp);
  PrintSamplingTable(mp, 4, x, outp);
  HorizontalRule(outp, 1);

  /*** Print energy adjustment constants         ***/
  /*** (otherwise these will never get noticed). ***/
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(3.)", 4, " ", 1, "Energy adjustment constants.  These "
	      "constants are included to bring the overall energy of quantum "
	      "mechanical and molecular mechanical calculations into general "
	      "agreement so that the fitted parameters can address the energy "
	      "differences between multiple states of each system.  Gross "
	      "differences between quantum and molecular mechanics are first "
	      "eliminated by subtracting off the average energy of all states "
	      "for each system; these parameters are then included to \"mop "
	      "up.\"\n", 74, 0, outp);
  if (mp->wtfloor <= 100.0) {
    PrintVADesc(0, " ", 4, " ", 1, "Before fitting parameters, the importance "
		"of each system was measured as a function of the target "
		"energy, the lowest energy gaining the most importance and "
		"higher energies receiving less importance down to a minimum "
		"weight specified in the input.  This expanded table includes "
		"the mamimum weights applied to any one conformation of each "
		"system.\n", 74, 0, outp);
  }
  PrintEAConstants(mp, x, outp);
  HorizontalRule(outp, 1);

  /*** Print an analysis of the fitting matrix, ***/
  /*** checking for bad conditioning.           ***/
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(4.)", 4, " ", 1, "Fitting matrix analysis.  This is a "
	      "check against features of the fitting set that might create a "
	      "poorly conditioned matrix.\n", 74, 0, outp);
  PrintMatrixAnalysis(mp, x, outp);
  HorizontalRule(outp, 1);

  /*** Print the contributions of each parameter ***/
  /*** to the energy of each system              ***/
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(5.)", 4, " ", 1, "Energy contributions of each fitted "
	      "parameter.  All parameters that contribute to each system are "
	      "printed, this time with respect to the amount of energy that "
	      "they contribute to the molecular mechanical system.\n", 74, 0,
	      outp);
  PrintEnergyContrib(mp, A, x, outp);
  HorizontalRule(outp, 1);

  /*** Print the contexts in which each parameter ***/
  /*** is sampled; this expounds on section (2.). ***/
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(6.)", 4, " ", 1, "Context of each parameter.  Instances of "
              "each parameter found in the fitting conformations are "
              "enumerated below.\n", 74, 0, outp);
  ParameterDescription(mp, outp);
  HorizontalRule(outp, 1);

  /*** If changes have been made to any atom types, describe them. ***/
  if (mp->nchng > 0) {
    HorizontalRule(outp, 1);
    PrintVADesc(0, "(7.)", 4, " ", 1, "Changes mades to atom types in each "
		"system's topology file.  These changes were made to expand "
		"or reduce the parameters available for fitting.  The "
		"resulting parameter file must be paired with library files "
		"which reflect all of these changes, and a leap source file "
		"which includes all new atom types.\n", 74, 0, outp);
    PrintChangeLog(mp, outp);
    HorizontalRule(outp, 1);
  }
  fclose(outp);

  /*** FIX ME!!!  The CheckMinima function may not work so well. ***/
  CheckMinima(mp, tj, x);

  /*** Print the frcmod file for this parameter set ***/
  SetParmOutputFlags(mp);
  outp = FOpenSafe(tj->dumpname, tj->OverwriteOutput);
  fprintf(outp, "%s\n", mp->ititl);
  PrintParmAtoms(mp, outp);
  if (mp->Hydrophilics.row > 0) {
    for (i = 0; i < mp->Hydrophilics.row; i++) {
      fprintf(outp, "%.4s", mp->Hydrophilics.map[i]);
    }
    fprintf(outp, "\n");
  }
  PrintParmBond(mp, 2, x, outp);
  PrintParmBond(mp, 3, x, outp);
  PrintParmBond(mp, 4, x, outp);
  PrintParmBond(mp, 5, x, outp);
  PrintParmBond(mp, 9, x, outp);
  for (i = 0; i < mp->neqgroups; i++) {
    for (j = 0; j < mp->eqgroups[i].natom; j++) {
      fprintf(outp, "%.4s", &mp->eqgroups[i].types[4*j]);
    }
    fprintf(outp, "\n");
  }
  fprintf(outp, "\nMOD4      RE\n");
  for (i = 0; i < mp->natom; i++) {
    fprintf(outp, "  %.2s             %10.6lf  %10.6lf  %s\n",
	    mp->atoms[i].atype, mp->atoms[i].ljsig, mp->atoms[i].ljeps,
	    mp->atoms[i].comment);
  }
  fprintf(outp, "\nEND\n");
  fclose(outp);
}
