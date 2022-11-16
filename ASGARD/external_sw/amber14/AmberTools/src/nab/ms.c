/* ms.f -- translated by f2c (version 20030306).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "../f2c/f2c.h"

/* Common Block Declarations */

struct {
    real colpre[3], radpre, r2pre;
} mpck_;

#define mpck_1 mpck_

/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__120 = 120;
static integer c__0 = 0;
static integer c__127 = 127;
static real c_b36 = 0.f;
static integer c__150 = 150;
static integer c__160 = 160;
static integer c__170 = 170;
static integer c__140 = 140;
static integer c__130 = 130;
static integer c__210 = 210;
static integer c__2000 = 2000;
static integer c__320 = 320;
static integer c__720 = 720;
static integer c__440 = 440;
static integer c__480 = 480;
static integer c__760 = 760;
static integer c__830 = 830;
static integer c__850 = 850;

/* molecular surface program */
/* ms */

/* December 16, 1983 */

/* Copyright c 1983 */
/* by Michael Connolly */

/* Written by Michael Connolly */

/* References: */

/* M.L. Connolly, "Solvent-accessible surfaces of proteins */
/* and nucleic acids", Science, 221, 709-713 (1983). */

/* M.L. Connolly, "Analytical molecular surface calculation", */
/* Journal of Applied Crystallography, 16, 548-558 (1983). */

/* This program may be freely distributed to anyone. */

/* It is written in fortran 77. */

/* ms calculates the molecular surface of a molecule */
/* given the coordinates of its atoms.  van der waals radii for */
/* the atoms and the probe radius must also be specified. */

/* The term molecular surface was introduced by F.M. Richards */
/* (Annual Reviews of Biophysics and Bioengineering, 1977, */
/* pages 151-176) with a specific meaning.  The surface */
/* Richards defined consists of two parts:  (1) the contact */
/* surface and (2) the reentrant surface.  He defines the */
/* contact surface to be that part of the van der waals */
/* surface of each atom which is accessible to a probe sphere */
/* of a given radius.  He defines the reentrant surface to be */
/* the inward-facing part of the probe sphere when it is */
/* simultaneously in contact with more than one atom. */

/* In implementing this definition I have found that there are */
/* two kinds of reentrant surface:  (1) concave reentrant */
/* surface, which is generated when the probe sphere */
/* simultaneously touches three atoms and (2) saddle-shaped */
/* reentrant surface, which is generated as the probe sphere */
/* rolls along the crevice between two atoms. */
/* I have also found that reentrant surface belonging to one */
/* probe may be contained in the interior volume of an */
/* overlapping probe and so must be removed. */

/* The input to this program consists of three files. */

/* The first file contains one record specifying the requested */
/* density of surface points (the number of surface points */
/* per unit area), the probe radius, the buried surface */
/* flag and the ascii/binary long/short output flag. */
/* Normally the buried surface flag will be blank or zero. */
/* If this flag is equal to 1, then the */
/* only surface calculated will be the surface of each */
/* molecule that is buried by the other molecules. */
/* If this flag is equal to 2, both buried and unburied */
/* surface are calculated, but the buried surface points */
/* are flagged by a 1 in the buried output field to the */
/* right of the surface normal. */
/* The ascii/binary long/short output flag */
/* may have one of these four values: */

/* 0      ascii       long */
/* 1      binary      long */
/* 2      ascii       short */
/* 3      binary      short */

/* The various output formats are discussed below. */

/* The format of the parameter record is: (2f10.5,2i5) */

/* The second file contains atomic radius records. */
/* each record has an integer that is the atom type, */
/* a van der waals radius for that type, */
/* and the point density for that type [added by dac 12/95!]. */
/* The atom types need not be contiguous or in */
/* increasing order. the format is: (i5,2f10.5). */

/* The third file contains the atomic coordinate records. */
/* Each atomic coordinate record has the x, y and z coordinates */
/* of the atoms, followed by the atom type, a surface */
/* request number and a molecule number. */
/* The format is: (3f10.5,3i5). */
/* The surface request number may be 0, 1 or 2. */
/* 0 means that the atom is to be ignored, */
/* except for occupying a place in the sequence of atom numbers. */
/* 1 means that no surface is to be generated for this atom, */
/* but probe positions that collide with this atom will still */
/* be disallowed.  2 means that we are requested to generate */
/* surface for this atom.  In most cases 2 should be specified */
/* for all atoms. The molecule number is an integer that is */
/* used to divide the atoms into groups so that each group */
/* may be given its own surface. These groups need not */
/* correspond to actual molecules, as the program knows */
/* nothing about bonding. The characters to the right of */
/* these six fields are not read and may contain anything. */

/* The output from this program consists of two files, */
/* called 'contact' and 'reentrant', containing the contact */
/* and reentrant surface points.  All lines in both files have */
/* the same format.  The first three fields are the atom */
/* numbers of the atoms the probe was touching when it */
/* generated the given surface point.  The first number is */
/* the atom whose van der waals surface the point is closest to. */
/* The fourth field is the number of atoms the probe was */
/* touching. If this number is less then three, one or both */
/* of the second and third fields will be zero. The fifth, sixth */
/* and seventh fields are the coordinates of the surface point. */
/* The eighth field is the molecular area associated with */
/* the surface point.  The ninth, tenth and eleventh fields */
/* are a unit vector perpendicular to the surface pointing in */
/* the direction of the probe center. */
/* If the buried surface flag equals 2, there is a '1' written */
/* at the end of the record for buried surface points, and */
/* a '0' for accessible points. The output format is: */
/* (3i5,i2,3f9.3,4f7.3,i2) */
/* There is also a short output format: (3i5,i2,3f9.3,i2). */
/* The corresponding binary records are: */
/* (4i*2,7r*4,i*2) and (4i*2,3r*4,i*2) */

/* The contact file is generated in atom number order. */
/* The reentrant file should be sorted */
/* into atom number order and then merged with the contact file */
/* using the sort/merge utility programs available at your */
/* installation.  Then all the surface points belonging to a */
/* given atom will be in a contiguous series of records. */


/* The program has the ability to calculate the */
/* van der waals surface */
/* of a molecule simply by specifying a probe radius of zero. */
/* The reentrant code is bypassed completely and there are no */
/* before and reentrant files. For a probe radius of zero */
/* the van der waals surface and the contact surface */
/* are equivalent. */

/* The flow of the program may be described in general terms */
/* as follows. First all the input is read. Then the contact */
/* and reentrant surface is generated. The contact surface */
/* is in its final form, but the reentrant surface is */
/* written to a temporary file, called 'before'. */
/* Each reentrant probe position is written to this file, */
/* followed by all its surface points. */
/* After all the contact and reentrant surface has been generated, */
/* the 'before' file is read and a final reentrant surface file is */
/* written which contains all the reentrant surface points */
/* not lying within any reentrant probe. */


/* Sometimes it is desirable to calculate the surface of only */
/* part of a molecule. One cannot simply remove the remaining */
/* atoms from the input file as this will generate false */
/* surfaces. One could calculate a surface for the entire */
/* molecule and then edit out that part of the surface belonging */
/* to the atoms of interest, but this is a needless */
/* waste of computer time.  The proper way to accomplish this */
/* task is to place the probe next to the atoms whose surface */
/* is requested and to use the remaining atoms only for */
/* collision checks so that false surfaces will not be */
/* generated. When the probe is placed next to two or three */
/* atoms, only one of these need be an atom whose surface */
/* is requested. only that part of the arc or spherical */
/* triangle that belongs to the atoms whose surface is */
/* requested will be written as output. */

/* This program allows the calculation of individual surfaces */
/* for several interacting molecules. this may be done with */
/* or without the buried surface option. */

/* This program should be compiled with integer*4 as the default. */
/* It uses no include files or libraries. */

/* Main program */ int main( int argc, char* argv[] )
{
    /* Format strings */
    static char fmt_3[] = "(\002Usage: ms -d <#> -rp <#> -pqr <file> \002)";
    static char fmt_150[] = "(4x,i7,20x,3f8.3,2f8.4)";
    static char fmt_1150[] = "(1x,\002atom\002,i5,\002 dropped (same co as"
	    " \002,i5,\002)\002)";
    static char fmt_1650[] = "(1x,i5,1x,\002atoms\002)";
    static char fmt_1700[] = "(1x,i5,1x,\002omitted\002)";
    static char fmt_1750[] = "(1x,i5,1x,\002collision only\002)";
    static char fmt_1800[] = "(1x,i5,1x,\002surface\002)";
    static char fmt_1850[] = "(1x,\002surface point density = \002,f10.5,5x"
	    ",\002probe radius = \002,f10.5)";
    static char fmt_1900[] = "(1x,\002buried surface only\002)";
    static char fmt_1950[] = "(1x,\002buried surface flagged\002)";
    static char fmt_2000[] = "(1x,\002binary contact and reentrant files\002)"
	    ;
    static char fmt_2050[] = "(1x,\002short output records\002)";
    static char fmt_3150[] = "(1x,\002atoms\002,2i5,\002 have the same cen"
	    "ter\002)";
    static char fmt_3500[] = "(1x,\002atoms\002,3i5,\002 have concentric cir"
	    "cles\002)";
    static char fmt_5500[] = "(3i5,i2,3f9.3,4f7.3,i2)";
    static char fmt_5550[] = "(3i5,i2,3f9.3,i2)";
    static char fmt_5700[] = "(1x,i5,\002 neighbors maximum\002)";
    static char fmt_7900[] = "(1x,i5,\002 yon and \002,i5,\002 victim probe"
	    "s\002)";
    static char fmt_7950[] = "(1x,i5,\002 saddle and \002,i5,\002 concave su"
	    "rface points removed during non-symmetry orsr\002)";
    static char fmt_8050[] = "(1x,i5,\002 contact and \002,i5,\002 saddle an"
	    "d \002,i5,\002 concave surface points\002)";
    static char fmt_8100[] = "(1x,i8,\002 total surface points\002)";
    static char fmt_8150[] = "(1x,\002contact area:\002,f10.3,2x,\002reentra"
	    "nt area:\002,f10.3,2x,\002total area:\002,f10.3)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3;
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), f_open(olist 
	    *), f_rew(alist *), s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer do_fio(integer *, char *, ftnlen), s_rsfe(cilist *), e_rsfe(void),
	     f_clos(cllist *);
    double sqrt(doublereal), cos(doublereal), sin(doublereal);
    integer s_wsue(cilist *), do_uio(integer *, char *, ftnlen), e_wsue(void),
	     s_rsue(cilist *), e_rsue(void);

    /* Local variables */
    static real a[3], d__[6001], f, g[9]	/* was [3][3] */, h__[9]	
	    /* was [3][3] */;
    static integer i__, j, k, l, n;
    static real p[3], q[3], s[3000]	/* was [3][1000] */, t[3], x, y, z__, 
	    d2, f1, f2;
    static integer i1, j1, k1, l1, l2, n1[1000], n2[1000], n3[1000], ib;
    static real ci[3], cj[3], ck[3], fi, dk, co[45003]	/* was [3][15001] */, 
	    ar, dp, ua[18000000]	/* was [3][1000][6000] */, av[18000]	
	    /* was [3][6000] */, ay[36000]	/* was [3][12000] */, ri;
    static logical si, sj, sk;
    static integer np;
    static real rj, rp, rk;
    static integer ip;
    static real up[18000000]	/* was [3][1000][6000] */, pv[18000]	/* 
	    was [3][6000] */;
    static integer ny;
    static real ht, py[36000]	/* was [3][12000] */;
    static integer iy, nv, jp;
    static real dp2, rk2, rp2, aij[3], bij[3];
    static integer ici;
    static real rad[15000];
    static integer icj, ick;
    static real eat[3000]	/* was [3][1000] */;
    static integer ico[45000]	/* was [3][15000] */;
    static real eva[18000000]	/* was [3][1000][6000] */;
    static integer ias[15001], iat[15001];
    static char arg[60];
    static integer jck, jcj, jci;
    static real pij[3];
    static integer nua[6000], idx;
    static real uij[3], vij[3], vbs[3];
    static integer lkf;
    static real dij, vpi[3], vpj[3], vpk[3], hij;
    static integer nup[6000];
    static logical sns;
    static real sum, pow[9]	/* was [3][3] */;
    static logical yon[1000], srs[15000];
    extern doublereal det_(real *, real *, real *), dot_(real *, real *);
    static real dsi, dsj, dsk, rij, avh;
    extern /* Subroutine */ int cat_(real *, real *);
    static integer ipt;
    static real vbs0[6000]	/* was [3][1000][2] */, vps0[3000]	/* 
	    was [3][1000] */, arca[1000], area, dfar, aijk[6]	/* was [3][2] 
	    */, bijk[3], cijk[3], dijk;
    static integer iarg;
    static real aijp[6]	/* was [3][2] */, cnbr[6000]	/* was [3][2000] */, 
	    hijk;
    static integer narc, neat;
    static real ghgt[9]	/* was [3][3] */;
    static integer inbr[2000];
    static logical pair[2], both;
    static real pijk[6]	/* was [3][2] */;
    static logical mnbr[2000];
    static integer ivic[6000], imol, nias[3];
    static real uijk[3], rnbr[2000];
    static logical snbr[2000];
    static real vijk[3];
    static integer itnl[2000];
    static real pijp[6]	/* was [3][2] */;
    static logical ayon[1000];
    static real dens;
    static integer indx, isph;
    static real pipt[3];
    static integer nnbr, iptr;
    static real sumi;
    static integer iuse;
    static logical bury;
    static integer jmin, jnbr, knbr;
    static real rijk, sign;
    extern doublereal dist_(real *, real *);
    static integer nrot;
    extern /* Subroutine */ int conj_(real *, real *, real *);
    static real duij;
    static integer ptyp, irot;
    static real vect1, vect2, vect3;
    extern doublereal dist2_(real *, real *);
    static real areac;
    static real angle;
    static integer icube[64000]	/* was [40][40][40] */, iabls;
    static real arear;
    static integer ihash, ifree, ncirc[6000];
    static logical scube[64000]	/* was [40][40][40] */;
    static real comin[3];
    static integer maxnb, jmold;
    static real ernbr[2000];
    static integer lknbr[2000], ivicp[6000];
    static logical found;
    extern /* Subroutine */ int genun_(real *, integer *);
    static integer iatom, natom;
    static real width;
    static integer nimol, jatom;
    extern doublereal anorm_(real *);
    static integer katom;
    extern /* Subroutine */ int imatx_(real *);
    static integer iprev, nyeat;
    static real outco[3];
    static integer ibury, itype[6001];
    extern /* Subroutine */ int error_(integer *, integer *, real *), vperp_(
	    real *, real *), cross_(real *, real *, real *);
    static logical short__;
    static integer ntype, nlost[3];
    extern /* Subroutine */ int vnorm_(real *, real *);
    static real rtype[6001];
    extern /* Subroutine */ int multv_(real *, real *, real *);
    static real torus[1000], circle[18000000]	/* was [3][1000][6000] */;
    extern logical collid_(real *, real *, real *, real *, logical *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *);
    static integer ishape;
    extern logical buried_(real *, real *, real *, real *, logical *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *);
    static real radmax;
    static integer nshape[3];
    static real disnbr[2000];
    static logical binary, sscube[64000]	/* was [40][40][40] */;
    static integer molnbr[2000], molvic[6000], jminbr, iatnum[3];
    static real vector[3];
    static integer ifrlst, icuptr[15000], molnum[15001];
    static real outvec[3];
    static logical yonprb;
    static integer mutual, molyon[12000];

    /* Fortran I/O blocks */
    static icilist io___11 = { 0, arg, 0, 0, 60, 1 };
    static icilist io___12 = { 0, arg, 0, 0, 60, 1 };
    static icilist io___13 = { 0, arg, 0, 0, 60, 1 };
    static icilist io___14 = { 0, arg, 0, 0, 60, 1 };
    static icilist io___16 = { 0, arg, 0, 0, 60, 1 };
    static cilist io___17 = { 0, 6, 0, fmt_3, 0 };
    static cilist io___18 = { 0, 6, 0, "(/,5x,a,a)", 0 };
    static cilist io___19 = { 0, 6, 0, fmt_3, 0 };
    static cilist io___26 = { 0, 2, 1, fmt_150, 0 };
    static cilist io___51 = { 0, 6, 0, fmt_1150, 0 };
    static cilist io___61 = { 0, 6, 0, fmt_1650, 0 };
    static cilist io___62 = { 0, 6, 0, fmt_1700, 0 };
    static cilist io___63 = { 0, 6, 0, fmt_1750, 0 };
    static cilist io___64 = { 0, 6, 0, fmt_1800, 0 };
    static cilist io___65 = { 0, 6, 0, fmt_1850, 0 };
    static cilist io___66 = { 0, 6, 0, fmt_1900, 0 };
    static cilist io___67 = { 0, 6, 0, fmt_1950, 0 };
    static cilist io___68 = { 0, 6, 0, fmt_2000, 0 };
    static cilist io___69 = { 0, 6, 0, fmt_2050, 0 };
    static cilist io___76 = { 0, 7, 0, "(i5)", 0 };
    static cilist io___121 = { 0, 6, 0, fmt_3150, 0 };
    static cilist io___145 = { 0, 6, 0, fmt_3500, 0 };
    static cilist io___170 = { 0, 4, 0, 0, 0 };
    static cilist io___174 = { 0, 4, 0, 0, 0 };
    static cilist io___198 = { 0, 4, 0, 0, 0 };
    static cilist io___199 = { 0, 4, 0, 0, 0 };
    static cilist io___205 = { 0, 7, 0, fmt_5500, 0 };
    static cilist io___206 = { 0, 7, 0, 0, 0 };
    static cilist io___207 = { 0, 7, 0, "(a5,4f8.3)", 0 };
    static cilist io___208 = { 0, 7, 0, 0, 0 };
    static cilist io___209 = { 0, 6, 0, fmt_5700, 0 };
    static cilist io___217 = { 0, 4, 1, 0, 0 };
    static cilist io___221 = { 0, 4, 0, 0, 0 };
    static cilist io___230 = { 0, 4, 1, 0, 0 };
    static cilist io___237 = { 0, 4, 0, 0, 0 };
    static cilist io___238 = { 0, 8, 0, fmt_5500, 0 };
    static cilist io___239 = { 0, 8, 0, 0, 0 };
    static cilist io___240 = { 0, 8, 0, fmt_5550, 0 };
    static cilist io___241 = { 0, 8, 0, 0, 0 };
    static cilist io___242 = { 0, 6, 0, fmt_7900, 0 };
    static cilist io___243 = { 0, 6, 0, fmt_7950, 0 };
    static cilist io___244 = { 0, 6, 0, fmt_8050, 0 };
    static cilist io___245 = { 0, 6, 0, fmt_8100, 0 };
    static cilist io___246 = { 0, 6, 0, fmt_8150, 0 };



/* maxatm     maximum number of atoms */
/* maxtyp     maximum number of atom types */
/* maxnbr     maximum number of neighbors an atom may have */
/* maxsph     maximum number of surface points on a sphere */
/* maxcir     maximum number of surface points on a circle */
/* maxarc     maximum number of surface points on an arc */
/* maxppp     maximum number of surface points per probe */
/* maxyon     maximum number of yon probes */
/* maxvic     maximum number of victim probes */
/* maxeat     maximum number of eaters of a probe's surface */
/* maxcub     maximum number of cubes in one direction */

/* maxatm must be greater than or equal to maxyon */
/* because they share the same cubing arrays */



/* run-time options and parameters */
/* rp	probe radius, set by -rp VALUE */
/* ar	additional radius for solvent acc surface with rp=0 */
/*  , set by -ar VALUE */
/* d	dot density, set by -d VALUE */
/* ibury	switch for turning on buried surface, set by -b or -bonly */
/* binary     binary output records instead of ascii, set by -bin */
/* short      short output records instead of ordinary long, set by -short */


/* atom type arrays */

/* itype     atom type number */
/* rtype     van der waals radius */
/* nua       number of unit vectors on sphere */

/* dimensioned one more in case input file is too long */

/* arrays for all atoms */

/* co        atomic coordinates */
/* ias       surface request number */
/* iat       atom itype */
/* molnum    molecule number */
/* rad       radius */
/* srs       some reentrant surface */

/* dimension arrays 1 more in case input file is too long */

/* cube arrays */

/* ico       integer cube coordinates */
/* icuptr    pointer to next atom in cube */
/* comin     minimum atomic coordinates (cube corner) */
/* icube     pointer to first atom in list for cube */
/* scube     some atom with srn=2 in cube */
/* sscube    some atom with srn=2 in cube or adjoining cubes */


/* neighbor arrays */

/* inbr      atom number */
/* cnbr      coordinates */
/* rnbr      radius */
/* snbr      true if srn = 2 for neighbor */
/* mnbr      mutual neighbor of iatom and jatom */
/* molnbr    molecule number */
/* ernbr     expanded radius (rnbr + rp) */
/* disnbr    distance from neighbor to iatom */
/* lknbr     link to next farthest out neighbor */
/* itnl      temporary neighbor list (before sort) */


/* circle and sphere unit vector arrays */

/* up        unit vectors for probe */
/* ua        unit vectors for atom */
/* eva       extended vectors for atom */
/* circle    points on a circle */


/* ci        coordinates of atom i */
/* cj        coordinates of atom j */
/* ck        coordinates of atom k */
/* si        srn = 2 for atom i */
/* sj        srn = 2 for atom j */
/* sk        srn = 2 for atom k */


/* geometric construction vectors */

/* vij       vector from atom i to atom j */
/* uij       unit vector of vij */
/* q,t       two perpendicular vectors in the saddle plane */
/* cijk      center of circle of intersection of expanded sphere */
/*           of atom k with the saddle plane of atoms i and j */
/* vijk      vector from torus center to cijk */
/* uijk      unit vector of vijk */
/* bij       torus center */
/* aij       starting altitude vector from torus center */
/* bijk      base point for atoms i, j and k */
/* aijk      altitude vector for atoms i, j and k */
/* aijp      altitude vector to probe (rotated) */
/* a         altitude vector, general, used in orsr */
/* p         probe coordinates (general, used in orsr) */
/* pijp      center of probe placed tangent to atoms i and j */
/* pij       starting center of probe tangent to atoms i and j */
/* pijk      probe placed tangent to atoms i, j and k */
/* pipt      probe placed tangent to atom i */
/* vpi       vector from probe center to contact point with atom i */
/* vpj       vector from probe center to contact point with atom j */
/* vpk       vector from probe center to contact point with atom k */
/* vps0      starting arc points relative to probe center */
/* vbs0      starting arc points relative to torus center */
/* vbs       rotated arc point relative to torus center */
/* arca      area of arc point */
/* ayon      arc point on yon side of symmetry element */
/* vector    temporary vector storage */



/* g         uij, q, t frame for torus */
/* h         rotation about x-axis */
/* ghgt      rotation about uij axis */
/* pow       powers of ghgt */

/* rotation matrices */

/* s         surface points for probe */
/* torus     surface points for torus */
/* n1        atom surface point is on or closest to */
/* n2,n3     other atoms probe is touching */
/* yon       whether point lies on yon side of symmetry element */

/* reentrant probe record */

/* both      both probe positions free of collisions */
/* pair      this member of pair free from collisions */
/* yonprb    probe crosses symmetry element */
/* found     search flag */

/* logical variables */

/* orsr for non-symmetry-related probes */
/* the factor of three for the victim arrays */
/* is based upon experience */

/* py        center of yon probe */
/* ay        altitude vector of yon probe */
/* pv        center of victim probe */
/* av        altitude vector of victim probe */
/* ivic      list of probe numbers of victims */
/* ivicp     pointer to next victim with same hash */
/* molyon    molecule number of yon probe */
/* molvic    molecule number of victim probe */
/* eat       coordinates of eaters of probe */


/* multiple-molecule variables */

/* bury      probe position is buried */


/* counters */

/* nias      number of atoms with given srn */
/* nshape    number of surface points with given shape */
/* nlost     number of surface points lost in non-symmetry orsr */


/* output variables */

/* iatnum   atom numbers */
/* ishape   surface shape number */
/* ib       buried surface flag */



/* MAS	common block for optimization */
/* 	real	colpre(3) */
/* 	real 	radpre,r2pre */



/* logical functions */


/* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */



/* Set default values for run options: */

    dens = 8.f;
    ar = 0.f;
    rp = 1.4f;
    ibury = 0;
    iabls = 0;
    short__ = FALSE_;
    binary = FALSE_;

/*    --- process command-line arguments: */

    iarg = 0;
    indx = argc;
    if (indx == 0) {
	goto L2;
    }
L1:
    ++iarg;
    /* getarg_(&iarg, arg, (ftnlen)60); */
    strncpy( arg, argv[iarg], 59 );  arg[59] = '\0';
    if (s_cmp(arg, "-d", (ftnlen)60, (ftnlen)2) == 0) {
	++iarg;
    /* getarg_(&iarg, arg, (ftnlen)60); */
    strncpy( arg, argv[iarg], 59 );  arg[59] = '\0';
	s_rsli(&io___11);
	do_lio(&c__4, &c__1, (char *)&dens, (ftnlen)sizeof(real));
	e_rsli();
    } else if (s_cmp(arg, "-rp", (ftnlen)60, (ftnlen)3) == 0) {
	++iarg;
    /* getarg_(&iarg, arg, (ftnlen)60); */
    strncpy( arg, argv[iarg], 59 );  arg[59] = '\0';
	s_rsli(&io___12);
	do_lio(&c__4, &c__1, (char *)&rp, (ftnlen)sizeof(real));
	e_rsli();
    } else if (s_cmp(arg, "-ar", (ftnlen)60, (ftnlen)3) == 0) {
	++iarg;
    /* getarg_(&iarg, arg, (ftnlen)60); */
    strncpy( arg, argv[iarg], 59 );  arg[59] = '\0';
	s_rsli(&io___13);
	do_lio(&c__4, &c__1, (char *)&ar, (ftnlen)sizeof(real));
	e_rsli();
    } else if (s_cmp(arg, "-dfar", (ftnlen)60, (ftnlen)5) == 0) {
	++iarg;
    /* getarg_(&iarg, arg, (ftnlen)60); */
    strncpy( arg, argv[iarg], 59 );  arg[59] = '\0';
	s_rsli(&io___14);
	do_lio(&c__4, &c__1, (char *)&dfar, (ftnlen)sizeof(real));
	e_rsli();
    } else if (s_cmp(arg, "-short", (ftnlen)60, (ftnlen)6) == 0) {
	short__ = TRUE_;
    } else if (s_cmp(arg, "-bin", (ftnlen)60, (ftnlen)4) == 0) {
	binary = TRUE_;
    } else if (s_cmp(arg, "-ibury", (ftnlen)60, (ftnlen)6) == 0) {
	++iarg;
    /* getarg_(&iarg, arg, (ftnlen)60); */
    strncpy( arg, argv[iarg], 59 );  arg[59] = '\0';
	s_rsli(&io___16);
	do_lio(&c__3, &c__1, (char *)&ibury, (ftnlen)sizeof(integer));
	e_rsli();
    } else if (s_cmp(arg, "-pqr", (ftnlen)60, (ftnlen)4) == 0) {
	++iarg;
    /* getarg_(&iarg, arg, (ftnlen)60); */
    strncpy( arg, argv[iarg], 59 );  arg[59] = '\0';
	o__1.oerr = 0;
	o__1.ounit = 2;
	o__1.ofnmlen = 60;
	o__1.ofnm = arg;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = "formatted";
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = 2;
	f_rew(&al__1);
    } else if (s_cmp(arg, "-help", (ftnlen)60, (ftnlen)5) == 0) {
	s_wsfe(&io___17);
	e_wsfe();
	s_stop("", (ftnlen)0);
    } else {
	if (s_cmp(arg, " ", (ftnlen)60, (ftnlen)1) == 0) {
	    goto L2;
	}
	s_wsfe(&io___18);
	do_fio(&c__1, "unknown flag: ", (ftnlen)14);
	do_fio(&c__1, arg, (ftnlen)60);
	e_wsfe();
	s_wsfe(&io___19);
	e_wsfe();
	s_stop("", (ftnlen)0);
    }
    if (iarg < indx) {
	goto L1;
    }

L2:

/* input value checking */
/* check for negative probe radius */
    if (rp < 0.f) {
	error_(&c__120, &c__0, &rp);
    }
/* check buried surface flag */
    if (ibury < 0 || ibury > 2) {
	error_(&c__127, &ibury, &c_b36);
    }


/* initialize srn counters and coordinate minima */
    for (k = 1; k <= 3; ++k) {
	nias[k - 1] = 0;
/* L350: */
    }
    for (k = 1; k <= 3; ++k) {
	comin[k - 1] = 1e6f;
/* L400: */
    }
/* initialization */
    ntype = 1;
    radmax = 0.f;
    n = 1;
/* read atom types with their radii */
L100:
    i__1 = s_rsfe(&io___26);
    if (i__1 != 0) {
	goto L300;
    }
    i__1 = do_fio(&c__1, (char *)&itype[ntype - 1], (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L300;
    }
    i__1 = do_fio(&c__1, (char *)&co[n * 3 - 3], (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L300;
    }
    i__1 = do_fio(&c__1, (char *)&co[n * 3 - 2], (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L300;
    }
    i__1 = do_fio(&c__1, (char *)&co[n * 3 - 1], (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L300;
    }
    i__1 = do_fio(&c__1, (char *)&d__[ntype - 1], (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L300;
    }
    i__1 = do_fio(&c__1, (char *)&rtype[ntype - 1], (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L300;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L300;
    }
    iat[n - 1] = itype[n - 1];
    ias[n - 1] = 2;
    molnum[n - 1] = 1;
    rtype[ntype - 1] += ar;

/* check for atom overflow */
    if (n > 15000) {
	error_(&c__150, &n, &c_b36);
    }
/* check surface request number */
    if (ias[n - 1] < 0 || ias[n - 1] > 2) {
	error_(&c__160, &n, &c_b36);
    }
/* check for new coordinate minima */
    for (k = 1; k <= 3; ++k) {
	if (co[k + n * 3 - 4] < comin[k - 1]) {
	    comin[k - 1] = co[k + n * 3 - 4];
	}
/* L550: */
    }
/* increment counters for each srn type */
    for (k = 1; k <= 3; ++k) {
	if (ias[n - 1] == k - 1) {
	    ++nias[k - 1];
	}
/* L600: */
    }
/* we don't care whether ignored atoms have a radius */
    if (ias[n - 1] == 0) {
	goto L700;
    }

/* look for this atom type number so we may assign a radius */
    found = FALSE_;
    i__1 = ntype;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iat[n - 1] != itype[i__ - 1]) {
	    goto L650;
	}
	found = TRUE_;
/* transfer radius from atom type to atom radius array */
	rad[n - 1] = rtype[i__ - 1];
/* end of atom type search loop */
L650:
	;
    }
/* check whether atom type found */
    if (! found) {
	error_(&c__170, &n, &c_b36);
    }
L700:
    if (d__[ntype - 1] <= 0.f) {
	d__[ntype - 1] = dens;
    }
/* atom radii must be zero or positive */
    if (rtype[ntype - 1] < 1e-6f) {
	error_(&c__140, &ntype, &rtype[ntype - 1]);
    }
/* check for atom type overflow */
    if (ntype > 6000) {
	error_(&c__130, &ntype, &c_b36);
    }
/* check for new maximum radius */
    if (rtype[ntype - 1] > radmax) {
	radmax = rtype[ntype - 1];
    }
/* number of unit vectors depends on sphere area and input density */
/* Computing 2nd power */
    r__1 = rtype[ntype - 1];
    nua[ntype - 1] = r__1 * r__1 * 12.566370616f * d__[ntype - 1];
/* decrease to array size if too large */
    if (nua[ntype - 1] > 1000) {
	nua[ntype - 1] = 1000;
    }
    if (nua[ntype - 1] < 1) {
	nua[ntype - 1] = 1;
    }
/* create unit vector arrays */
    genun_(&ua[(ntype * 1000 + 1) * 3 - 3003], &nua[ntype - 1]);
/* compute extended vectors for later probe placement */
    i__1 = nua[ntype - 1];
    for (isph = 1; isph <= i__1; ++isph) {
	for (k = 1; k <= 3; ++k) {
	    eva[k + (isph + ntype * 1000) * 3 - 3004] = (rtype[ntype - 1] + 
		    rp) * ua[k + (isph + ntype * 1000) * 3 - 3004];
/* L200: */
	}
/* L250: */
    }
/* one more atom type */
    ++ntype;
    ++n;
    natom = n;
    goto L100;
L300:
    cl__1.cerr = 0;
    cl__1.cunit = 2;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* decrement on end of file */
    --ntype;
    --natom;
/* calculate width of cube from maximum atom radius and probe radius */
    width = (radmax + rp) * 2;
/* ============================================================================== */

/*     set up cube arrays */
/*     first the integer coordinate arrays */
    i__1 = natom;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (k = 1; k <= 3; ++k) {
	    ico[k + i__ * 3 - 4] = (co[k + i__ * 3 - 4] - comin[k - 1]) / 
		    width + 1;
	    if (ico[k + i__ * 3 - 4] < 1) {
		s_stop("cube coordinate too small", (ftnlen)25);
	    }
	    if (ico[k + i__ * 3 - 4] > 40) {
		s_stop("cube coordinate too large", (ftnlen)25);
	    }
/* L800: */
	}
/* L850: */
    }

/* initialize head pointer and srn=2 arrays */
    for (k = 1; k <= 40; ++k) {
	for (j = 1; j <= 40; ++j) {
	    for (i__ = 1; i__ <= 40; ++i__) {
		icube[i__ + (j + k * 40) * 40 - 1641] = 0;
		scube[i__ + (j + k * 40) * 40 - 1641] = FALSE_;
		sscube[i__ + (j + k * 40) * 40 - 1641] = FALSE_;
/* L900: */
	    }
/* L950: */
	}
/* L1000: */
    }

/* initialize linked list pointers */
    i__1 = natom;
    for (i__ = 1; i__ <= i__1; ++i__) {
	icuptr[i__ - 1] = 0;
/* L1050: */
    }

/* set up head and later pointers for each atom */
    i__1 = natom;
    for (iatom = 1; iatom <= i__1; ++iatom) {
/* skip atoms with surface request numbers of zero */
	if (ias[iatom - 1] == 0) {
	    goto L1250;
	}
	i__ = ico[iatom * 3 - 3];
	j = ico[iatom * 3 - 2];
	k = ico[iatom * 3 - 1];
	if (icube[i__ + (j + k * 40) * 40 - 1641] <= 0) {
/*     first atom in this cube */
	    icube[i__ + (j + k * 40) * 40 - 1641] = iatom;
	} else {
/*     add to end of linked list */
	    iptr = icube[i__ + (j + k * 40) * 40 - 1641];
L1100:
/* check for duplicate coordinates */
	    if (molnum[iatom - 1] == molnum[iptr - 1] && dist2_(&co[iatom * 3 
		    - 3], &co[iptr * 3 - 3]) <= 0.f) {
		ias[iatom - 1] = 0;
		s_wsfe(&io___51);
		do_fio(&c__1, (char *)&iatom, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iptr, (ftnlen)sizeof(integer));
		e_wsfe();
		goto L1250;
	    }
	    if (icuptr[iptr - 1] <= 0) {
		goto L1200;
	    }
/* move on down the list */
	    iptr = icuptr[iptr - 1];
	    goto L1100;
L1200:
/* store atom number */
	    icuptr[iptr - 1] = iatom;
	}
/* check for surfaced atom */
	if (ias[iatom - 1] == 2) {
	    scube[i__ + (j + k * 40) * 40 - 1641] = TRUE_;
	}
L1250:
	;
    }

/* check for 3 x 3 x 3 with some srn = 2 */

    for (k = 1; k <= 40; ++k) {
	for (j = 1; j <= 40; ++j) {
	    for (i__ = 1; i__ <= 40; ++i__) {
		if (icube[i__ + (j + k * 40) * 40 - 1641] == 0) {
		    goto L1450;
		}
/* check whether this cube or any adjacent cube has srn = 2 */
		i__1 = k + 1;
		for (k1 = k - 1; k1 <= i__1; ++k1) {
		    if (k1 < 1 || k1 > 40) {
			goto L1400;
		    }
		    i__2 = j + 1;
		    for (j1 = j - 1; j1 <= i__2; ++j1) {
			if (j1 < 1 || j1 > 40) {
			    goto L1350;
			}
			i__3 = i__ + 1;
			for (i1 = i__ - 1; i1 <= i__3; ++i1) {
			    if (i1 < 1 || i1 > 40) {
				goto L1300;
			    }
			    if (scube[i1 + (j1 + k1 * 40) * 40 - 1641]) {
				sscube[i__ + (j + k * 40) * 40 - 1641] = 
					TRUE_;
			    }
L1300:
			    ;
			}
L1350:
			;
		    }
L1400:
		    ;
		}
L1450:
		;
	    }
/* L1500: */
	}
/* L1550: */
    }

/* initialization */
/* maximum number of neighbors any atom has */
    maxnb = 0;
/* numbers of surface points */
    for (k = 1; k <= 3; ++k) {
	nshape[k - 1] = 0;
	nlost[k - 1] = 0;
/* L1600: */
    }
/* number of yon probes */
    ny = 0;
/* contact and reentrant areas */
    areac = 0.f;
    arear = 0.f;

/* write out messages */
    s_wsfe(&io___61);
    do_fio(&c__1, (char *)&natom, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___62);
    do_fio(&c__1, (char *)&nias[0], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___63);
    do_fio(&c__1, (char *)&nias[1], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___64);
    do_fio(&c__1, (char *)&nias[2], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___65);
    do_fio(&c__1, (char *)&dens, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&rp, (ftnlen)sizeof(real));
    e_wsfe();
    if (ibury == 1) {
	s_wsfe(&io___66);
	e_wsfe();
    }
    if (ibury == 2) {
	s_wsfe(&io___67);
	e_wsfe();
    }
    if (binary) {
	s_wsfe(&io___68);
	e_wsfe();
    }
    if (short__) {
	s_wsfe(&io___69);
	e_wsfe();
    }
/* stop if density is not positive */
    if (dens <= 0.f) {
	s_stop("non-positive density", (ftnlen)20);
    }
/* skip probe and circle setup if van der waals surface */
    if (rp == 0.f) {
	goto L2150;
    }
    i__1 = ntype;
    for (ptyp = 1; ptyp <= i__1; ++ptyp) {

/* set up probe sphere and circle */

/* Computing 2nd power */
	r__1 = rp;
	nup[ptyp - 1] = r__1 * r__1 * 12.566370616f * d__[ptyp - 1];
	if (nup[ptyp - 1] < 1) {
	    nup[ptyp - 1] = 1;
	}
	if (nup[ptyp - 1] > 1000) {
	    nup[ptyp - 1] = 1000;
	}
	genun_(&up[(ptyp * 1000 + 1) * 3 - 3003], &nup[ptyp - 1]);
	ncirc[ptyp - 1] = rp * 6.2831853080000002f * sqrt(d__[ptyp - 1]);
	if (ncirc[ptyp - 1] < 1) {
	    ncirc[ptyp - 1] = 1;
	}
	if (ncirc[ptyp - 1] > 1000) {
	    ncirc[ptyp - 1] = 1000;
	}
	i__2 = ncirc[ptyp - 1];
	for (i__ = 1; i__ <= i__2; ++i__) {
	    fi = (i__ - 1) * 6.2831853080000002f / ncirc[ptyp - 1];
	    circle[(i__ + ptyp * 1000) * 3 - 3003] = rp * cos(fi);
	    circle[(i__ + ptyp * 1000) * 3 - 3002] = rp * sin(fi);
	    circle[(i__ + ptyp * 1000) * 3 - 3001] = 0.f;
/* L2100: */
	}
/* L2110: */
    }

/* open before file for writing */
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 6;
    o__1.ofnm = "before";
    o__1.orl = 0;
    o__1.osta = "unknown";
    o__1.oacc = 0;
    o__1.ofm = "unformatted";
    o__1.oblnk = 0;
    f_open(&o__1);
    al__1.aerr = 0;
    al__1.aunit = 4;
    f_rew(&al__1);
/* skip to here if no reentrant surface will be calculated */
L2150:

/* open contact file for writing */
    if (binary) {
	o__1.oerr = 0;
	o__1.ounit = 7;
	o__1.ofnmlen = 7;
	o__1.ofnm = "contact";
	o__1.orl = 0;
	o__1.osta = "unknown";
	o__1.oacc = 0;
	o__1.ofm = "unformatted";
	o__1.oblnk = 0;
	f_open(&o__1);
    } else {
	o__1.oerr = 0;
	o__1.ounit = 7;
	o__1.ofnmlen = 7;
	o__1.ofnm = "contact";
	o__1.orl = 0;
	o__1.osta = "unknown";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }
    al__1.aerr = 0;
    al__1.aunit = 7;
    f_rew(&al__1);
    s_wsfe(&io___76);
    do_fio(&c__1, (char *)&natom, (ftnlen)sizeof(integer));
    e_wsfe();

/* initialize some reentrant surface to false for each atom */
    i__1 = natom;
    for (iatom = 1; iatom <= i__1; ++iatom) {
	srs[iatom - 1] = FALSE_;
/* L2200: */
    }

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* big loop for each atom */
    i__1 = natom;
    for (iatom = 1; iatom <= i__1; ++iatom) {
/* skip ignored atoms */
	if (ias[iatom - 1] == 0) {
	    goto L5650;
	}
/* find the index into the atom type arrays for iatom */
	i__2 = ntype;
	for (idx = 1; idx <= i__2; ++idx) {
	    if (itype[idx - 1] == iat[iatom - 1]) {
		goto L2240;
	    }
/* L2230: */
	}
	s_stop("logic error in ms regarding atom types", (ftnlen)38);
L2240:
	ptyp = idx;

	ici = ico[iatom * 3 - 3];
	icj = ico[iatom * 3 - 2];
	ick = ico[iatom * 3 - 1];
/* skip iatom if its cube and adjoining cubes contain only blockers */
	if (! sscube[ici + (icj + ick * 40) * 40 - 1641]) {
	    goto L5650;
	}
/* transfer values from large arrays to iatom variables */
	ri = rad[iatom - 1];
	si = ias[iatom - 1] == 2;
	for (k = 1; k <= 3; ++k) {
	    ci[k - 1] = co[k + iatom * 3 - 4];
/* L2250: */
	}
	imol = molnum[iatom - 1];

/* gather the neighboring atoms of iatom */
/* initialize number of neighbors, and number of neighbors in the */
/* same molecule as atom i */
	nnbr = 0;
	nimol = 0;
/* initialize srn = 2 for some neighbor to false */
	sns = FALSE_;
/* save a little time for distance check */
	sumi = rp * 2 + ri;
/* check iatom cube and adjacent cubes for neighboring atoms */
	i__2 = ick + 1;
	for (jck = ick - 1; jck <= i__2; ++jck) {
	    if (jck < 1 || jck > 40) {
		goto L2550;
	    }
	    i__3 = icj + 1;
	    for (jcj = icj - 1; jcj <= i__3; ++jcj) {
		if (jcj < 1 || jcj > 40) {
		    goto L2500;
		}
		i__4 = ici + 1;
		for (jci = ici - 1; jci <= i__4; ++jci) {
		    if (jci < 1 || jci > 40) {
			goto L2450;
		    }
		    jatom = icube[jci + (jcj + jck * 40) * 40 - 1641];
L2300:
/* check for end of linked list for this cube */
		    if (jatom <= 0) {
			goto L2400;
		    }
/* distance check */
		    sum = sumi + rad[jatom - 1];
		    vect1 = (r__1 = co[jatom * 3 - 3] - ci[0], dabs(r__1));
		    if (vect1 >= sum) {
			goto L2350;
		    }
		    vect2 = (r__1 = co[jatom * 3 - 2] - ci[1], dabs(r__1));
		    if (vect2 >= sum) {
			goto L2350;
		    }
		    vect3 = (r__1 = co[jatom * 3 - 1] - ci[2], dabs(r__1));
		    if (vect3 >= sum) {
			goto L2350;
		    }
/* Computing 2nd power */
		    r__1 = vect1;
/* Computing 2nd power */
		    r__2 = vect2;
/* Computing 2nd power */
		    r__3 = vect3;
		    d2 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
/* Computing 2nd power */
		    r__1 = sum;
		    if (d2 >= r__1 * r__1) {
			goto L2350;
		    }
/* iatom is not its own neighbor */
		    if (iatom == jatom) {
			goto L2350;
		    }
/* we have a new neighbor */
		    ++nnbr;
/* check for neighbor overflow */
		    if (nnbr > 2000) {
			error_(&c__210, &nnbr, &c_b36);
		    }
/* save atom number in temporary array */
		    itnl[nnbr - 1] = jatom;
/* check whether surfaced neighbor in same molecule */
		    if (ias[jatom - 1] == 2 && molnum[jatom - 1] == imol) {
			sns = TRUE_;
		    }
/* count the number of atoms in the same molecule as iatom */
		    if (imol == molnum[jatom - 1]) {
			++nimol;
		    }
L2350:
/* get number of next atom in cube */
		    jatom = icuptr[jatom - 1];
		    goto L2300;
L2400:
L2450:
		    ;
		}
L2500:
		;
	    }
L2550:
	    ;
	}
/* keep track of maximum number of neighbors */
/* for array-dimensioning purposes */
	if (nnbr > maxnb) {
	    maxnb = nnbr;
	}

/* no surface for atom i if buried only flag set and */
/* there are no neighbors from the other molecules */
	if (ibury == 1 && nimol == nnbr) {
	    goto L5650;
	}

/* no surface if iatom and all neighbor atoms */
/* in the same molecule have surface request numbers < 2 */
	if (! si && ! sns) {
	    goto L5650;
	}

/* set up neighbors arrays with jatom in increasing order */

/* initialize minimum neighbor atom number */
	jmold = 0;
	i__2 = nnbr;
	for (iuse = 1; iuse <= i__2; ++iuse) {
	    jmin = natom + 1;
	    i__3 = nnbr;
	    for (jnbr = 1; jnbr <= i__3; ++jnbr) {
/* don't use ones already sorted */
		if (itnl[jnbr - 1] <= jmold) {
		    goto L2600;
		}
		if (itnl[jnbr - 1] < jmin) {
		    jmin = itnl[jnbr - 1];
		    jminbr = jnbr;
		}
L2600:
		;
	    }
	    jmold = jmin;
	    jnbr = jminbr;
	    jatom = itnl[jnbr - 1];
/* transfer atom number, coordinates, radius, surface request number, */
/* molecule number, expanded radius, distance from iatom */
	    inbr[iuse - 1] = jatom;
	    for (k = 1; k <= 3; ++k) {
		cnbr[k + iuse * 3 - 4] = co[k + jatom * 3 - 4];
/* L2650: */
	    }
	    rnbr[iuse - 1] = rad[jatom - 1];
	    snbr[iuse - 1] = ias[jatom - 1] == 2;
	    molnbr[iuse - 1] = molnum[jatom - 1];
	    ernbr[iuse - 1] = rnbr[iuse - 1] + rp;
	    disnbr[iuse - 1] = dist2_(ci, &cnbr[iuse * 3 - 3]);
/* initialize link to next farthest out neighbor */
	    lknbr[iuse - 1] = 0;
/* L2700: */
	}
/* set up a linked list of neighbors in order of */
/* increasing distance from iatom */
/* initialize pointer to first neighbor to 0 */
	lkf = 0;
/* look for neighbor in same molecule */
/* we want only atoms in same molecule for collision check */
	i__2 = nnbr;
	for (l = 1; l <= i__2; ++l) {
	    if (imol != molnbr[l - 1]) {
		goto L2750;
	    }
	    lkf = l;
	    goto L2800;
L2750:
	    ;
	}
	if (lkf == 0) {
	    goto L3000;
	}
L2800:
/* put remaining neighbors in linked list at proper position */
	i__2 = nnbr;
	for (l = lkf + 1; l <= i__2; ++l) {
	    if (imol != molnbr[l - 1]) {
		goto L2950;
	    }
	    l1 = 0;
	    l2 = lkf;
L2850:
	    if (disnbr[l - 1] < disnbr[l2 - 1]) {
		goto L2900;
	    }
	    l1 = l2;
	    l2 = lknbr[l2 - 1];
	    if (l2 != 0) {
		goto L2850;
	    }
L2900:
/* add to list */
	    if (l1 == 0) {
		lkf = l;
		lknbr[l - 1] = l2;
	    } else {
		lknbr[l1 - 1] = l;
		lknbr[l - 1] = l2;
	    }
L2950:
	    ;
	}
L3000:

/* no reentrant surface will be calculated if we are */
/* calculating the van der waals surface */
/* instead of the molecular surface */
	if (rp == 0.f) {
	    goto L5200;
	}
/* no reentrant surface if iatom has no neighbors */
	if (nimol <= 0) {
	    goto L5200;
	}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* medium loop for each neighbor of iatom */

	i__2 = nnbr;
	for (jnbr = 1; jnbr <= i__2; ++jnbr) {
	    jatom = inbr[jnbr - 1];

/* each pair of atoms is considered only once */
	    if (jatom <= iatom) {
		goto L5150;
	    }
/* each molecule gets a separate surface */
	    if (imol != molnbr[jnbr - 1]) {
		goto L5150;
	    }

/* tranfer from neighbor arrays to jatom variables */
	    rj = rnbr[jnbr - 1];
	    sj = snbr[jnbr - 1];
	    for (k = 1; k <= 3; ++k) {
		cj[k - 1] = cnbr[k + jnbr * 3 - 4];
/* L3050: */
	    }

/* here follow geometric calculations of points, vectors and */
/* distances used for probe placement in both saddle and */
/* concave reentrant surface generation */

/* calculate the intersection */
/* of the expanded spheres of iatom and jatom */
/* this circle is called the saddle circle */
/* the plane it lies in is called the saddle plane */

	    for (k = 1; k <= 3; ++k) {
		vij[k - 1] = cj[k - 1] - ci[k - 1];
/* L3100: */
	    }
/* create an orthonormal frame */
/* with uij pointing along the inter-atomic axis */
/* and q and t defining the saddle plane */
	    if (anorm_(vij) <= 0.f) {
		s_wsfe(&io___121);
		do_fio(&c__1, (char *)&iatom, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&jatom, (ftnlen)sizeof(integer));
		e_wsfe();
		goto L5150;
	    }
	    vnorm_(vij, uij);
	    vperp_(uij, q);
	    cross_(uij, q, t);

/* calculate the saddle circle center and radius */
	    dij = anorm_(vij);
/* Computing 2nd power */
	    r__1 = ri + rp;
/* Computing 2nd power */
	    r__2 = rj + rp;
/* Computing 2nd power */
	    r__3 = dij;
	    f = ((r__1 * r__1 - r__2 * r__2) / (r__3 * r__3) + 1.f) * .5f;
/* base point */
	    for (k = 1; k <= 3; ++k) {
		bij[k - 1] = ci[k - 1] + f * vij[k - 1];
/* L3200: */
	    }
/* Computing 2nd power */
	    r__1 = ri + rj + rp * 2;
/* Computing 2nd power */
	    r__2 = dij;
	    f1 = r__1 * r__1 - r__2 * r__2;
/* skip to bottom of middle loop if atoms are too far apart */
	    if (f1 <= 0.f) {
		goto L5150;
	    }
/* Computing 2nd power */
	    r__1 = dij;
/* Computing 2nd power */
	    r__2 = ri - rj;
	    f2 = r__1 * r__1 - r__2 * r__2;
/* skip to bottom of middle loop if one atom inside the other */
	    if (f2 <= 0.f) {
		goto L5150;
	    }
/* height (radius of saddle circle) */
	    hij = sqrt(f1 * f2) / (dij * 2);
/* a starting altitude */
	    for (k = 1; k <= 3; ++k) {
		aij[k - 1] = hij * q[k - 1];
/* L3250: */
	    }


/* concave reentrant surface */

/* gather mutual neighbors of iatom and jatom */
	    mutual = 0;
	    i__3 = nnbr;
	    for (knbr = 1; knbr <= i__3; ++knbr) {
		d2 = dist2_(cj, &cnbr[knbr * 3 - 3]);
/* Computing 2nd power */
		r__1 = rp * 2 + rj + rnbr[knbr - 1];
		mnbr[knbr - 1] = d2 < r__1 * r__1 && knbr != jnbr;
		if (mnbr[knbr - 1]) {
		    ++mutual;
		}
/* L3300: */
	    }

/* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* inner loop for each mutual neighbor of iatom and jatom */
	    ishape = 3;
	    i__3 = nnbr;
	    for (knbr = 1; knbr <= i__3; ++knbr) {
		if (! mnbr[knbr - 1]) {
		    goto L4200;
		}
		katom = inbr[knbr - 1];
/* iatom < jatom < katom */
		if (katom <= jatom) {
		    goto L4200;
		}
		sk = snbr[knbr - 1];
/* skip neighbor if all three atom not marked to be surfaced */
		if (! (si || sj || sk)) {
		    goto L4200;
		}
/* each molecule gets a separate surface */
		if (imol != molnbr[knbr - 1]) {
		    goto L4200;
		}

/* tranfer from neighbor array to katom variables */
		rk = rnbr[knbr - 1];
		for (k = 1; k <= 3; ++k) {
		    ck[k - 1] = cnbr[k + knbr * 3 - 4];
/* L3350: */
		}

/* calculate intersection of expanded sphere of katom */
/* with saddle plane. we will call this the katom circle. */

/* projection of vector, */
/* from katom to a point on the saddle plane, */
/* onto iatom-jatom axis, */
/* in order to get distance katom is from saddle plane */
		dk = uij[0] * (bij[0] - ck[0]) + uij[1] * (bij[1] - ck[1]) + 
			uij[2] * (bij[2] - ck[2]);

/* calculate radius of katom circle */
/* Computing 2nd power */
		r__1 = rk + rp;
/* Computing 2nd power */
		r__2 = dk;
		rijk = r__1 * r__1 - r__2 * r__2;
/* skip concave calculation if no intersection */
		if (rijk <= 0.f) {
		    goto L4200;
		}
		rijk = sqrt(rijk);
/* calculate center of katom circle */
		for (k = 1; k <= 3; ++k) {
		    cijk[k - 1] = ck[k - 1] + dk * uij[k - 1];
/* L3400: */
		}

/* calculate intersection of the katom circle with the saddle circle */
		for (k = 1; k <= 3; ++k) {
		    vijk[k - 1] = cijk[k - 1] - bij[k - 1];
/* L3450: */
		}
		dijk = anorm_(vijk);
		if (dijk <= 0.f) {
		    s_wsfe(&io___145);
		    do_fio(&c__1, (char *)&iatom, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&jatom, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&katom, (ftnlen)sizeof(integer));
		    e_wsfe();
		    goto L4200;
		}
/* Computing 2nd power */
		r__1 = hij;
/* Computing 2nd power */
		r__2 = rijk;
/* Computing 2nd power */
		r__3 = dijk;
		f = ((r__1 * r__1 - r__2 * r__2) / (r__3 * r__3) + 1.f) * .5f;
/* base point bijk is on symmetry plane and saddle plane */
		for (k = 1; k <= 3; ++k) {
		    bijk[k - 1] = bij[k - 1] + f * vijk[k - 1];
/* L3550: */
		}
/* Computing 2nd power */
		r__1 = hij + rijk;
/* Computing 2nd power */
		r__2 = dijk;
		f1 = r__1 * r__1 - r__2 * r__2;
/* skip to bottom of inner loop if katom too far away */
		if (f1 <= 0.f) {
		    goto L4200;
		}
/* Computing 2nd power */
		r__1 = dijk;
/* Computing 2nd power */
		r__2 = hij - rijk;
		f2 = r__1 * r__1 - r__2 * r__2;
/* skip to bottom of inner loop if katom circle inside saddle circle */
/* or vice-versa */
		if (f2 <= 0.f) {
		    goto L4200;
		}
		hijk = sqrt(f1 * f2) / (dijk * 2);
		vnorm_(vijk, uijk);
/* uij and uijk lie in the symmetry plane passing through the atoms */
/* so their cross product is perpendicular to this plane */
		cross_(uij, uijk, aijk);
/* two altitudes */
		for (k = 1; k <= 3; ++k) {
		    aijk[k - 1] = hijk * aijk[k - 1];
		    aijk[k + 2] = -aijk[k - 1];
/* L3600: */
		}

/* probe placement at ends of altitude vectors */
		for (ip = 1; ip <= 2; ++ip) {
		    for (k = 1; k <= 3; ++k) {
			pijk[k + ip * 3 - 4] = bijk[k - 1] + aijk[k + ip * 3 
				- 4];
/* L3650: */
		    }
/* collision check with mutual neighbors */
		    pair[ip - 1] = ! collid_(&pijk[ip * 3 - 3], &rp, cnbr, 
			    ernbr, mnbr, &nnbr, &c__2000, &ishape, &jnbr, &
			    knbr, molnbr, &imol, &lkf, lknbr);
/* L3700: */
		}
/* if neither probe position is allowed, skip to bottom of inner loop */
		if (! pair[0] && ! pair[1]) {
		    goto L4200;
		}
		both = pair[0] && pair[1];
/* some reentrant surface for all three atoms */
		srs[iatom - 1] = TRUE_;
		srs[jatom - 1] = TRUE_;
		srs[katom - 1] = TRUE_;

/* generate surface points */
/* Computing 2nd power */
		r__1 = rp;
		area = r__1 * r__1 * 12.56636f / nup[ptyp - 1];
		for (ip = 1; ip <= 2; ++ip) {
		    if (! pair[ip - 1]) {
			goto L4150;
		    }
/* give it some kind of value, in case we don't call buried */
		    bury = FALSE_;
/* only call buried if we care what the answer is */
		    if (ibury > 0) {
			bury = buried_(&pijk[ip * 3 - 3], &rp, cnbr, rnbr, 
				mnbr, &nnbr, &c__2000, &ishape, &jnbr, &knbr, 
				molnbr, &imol);
		    }
/* skip if not buried and buried surface only flag set */
		    if (ibury == 1 && ! bury) {
			goto L4150;
		    }
/* determine whether probe has surface on far side of plane */
		    yonprb = hijk < rp && ! both;
/* calculate vectors defining spherical triangle */
/* the vectors are given the probe radius as a length */
/* only for the purpose of making the geometry more clear */
		    for (k = 1; k <= 3; ++k) {
			vpi[k - 1] = (ci[k - 1] - pijk[k + ip * 3 - 4]) * rp /
				 (ri + rp);
			vpj[k - 1] = (cj[k - 1] - pijk[k + ip * 3 - 4]) * rp /
				 (rj + rp);
			vpk[k - 1] = (ck[k - 1] - pijk[k + ip * 3 - 4]) * rp /
				 (rk + rp);
/* L3750: */
		    }
		    sign = det_(vpi, vpj, vpk);
/* initialize number of surface points written */
		    np = 1;
/* gather points on probe sphere lying within triangle */
		    i__4 = nup[ptyp - 1];
		    for (i__ = 1; i__ <= i__4; ++i__) {
/* if the unit vector is pointing away from the symmetry plane */
/* the surface point cannot lie within the inward-facing triangle */
			if (dot_(&up[(i__ + ptyp * 1000) * 3 - 3003], &aijk[
				ip * 3 - 3]) > 0.f) {
			    goto L4000;
			}
			if (sign * det_(&up[(i__ + ptyp * 1000) * 3 - 3003], 
				vpj, vpk) < 0.f) {
			    goto L4000;
			}
			if (sign * det_(vpi, &up[(i__ + ptyp * 1000) * 3 - 
				3003], vpk) < 0.f) {
			    goto L4000;
			}
			if (sign * det_(vpi, vpj, &up[(i__ + ptyp * 1000) * 3 
				- 3003]) < 0.f) {
			    goto L4000;
			}
			if (np > 1000) {
			    error_(&c__320, &np, &c_b36);
			}
/* calculated whether point is on yon side of plane */
			yon[np - 1] = aijk[ip * 3 - 3] * (aijk[ip * 3 - 3] + 
				up[(i__ + ptyp * 1000) * 3 - 3003]) + aijk[ip 
				* 3 - 2] * (aijk[ip * 3 - 2] + up[(i__ + ptyp 
				* 1000) * 3 - 3002]) + aijk[ip * 3 - 1] * (
				aijk[ip * 3 - 1] + up[(i__ + ptyp * 1000) * 3 
				- 3001]) < 0.f;
/* overlapping reentrant surface removal */
/* for symmetry-related probe positions */
			if (yon[np - 1] && both) {
			    goto L4000;
			}
/* calculate coordinates of surface point */
			for (k = 1; k <= 3; ++k) {
			    s[k + np * 3 - 4] = pijk[k + ip * 3 - 4] + up[k + 
				    (i__ + ptyp * 1000) * 3 - 3004] * rp;
/* L3800: */
			}
/* find the closest atom and put the three atom numbers */
/* in the proper order */
/* n1 is closest, n2 < n3 */
			dsi = dist_(&s[np * 3 - 3], ci) - ri;
			dsj = dist_(&s[np * 3 - 3], cj) - rj;
			dsk = dist_(&s[np * 3 - 3], ck) - rk;
			if (dsi <= dsj && dsi <= dsk) {
			    goto L3850;
			}
			if (dsj <= dsi && dsj <= dsk) {
			    goto L3900;
			}
			if (! sk) {
			    goto L4000;
			}
			n1[np - 1] = katom;
			n2[np - 1] = iatom;
			n3[np - 1] = jatom;
			goto L3950;
L3850:
			if (! si) {
			    goto L4000;
			}
			n1[np - 1] = iatom;
			n2[np - 1] = jatom;
			n3[np - 1] = katom;
			goto L3950;
L3900:
			if (! sj) {
			    goto L4000;
			}
			n1[np - 1] = jatom;
			n2[np - 1] = iatom;
			n3[np - 1] = katom;
L3950:
			++np;
/* end of nup loop */
L4000:
			;
		    }
		    --np;
/* skip the write if no points */
		    if (np <= 0) {
			goto L4150;
		    }


/* write the molecule number, shape, number of points, */
/* probe position and */
/* the vector from the base to the probe center */
		    s_wsue(&io___170);
		    do_uio(&c__1, (char *)&imol, (ftnlen)sizeof(integer));
		    do_uio(&c__1, (char *)&ishape, (ftnlen)sizeof(integer));
		    do_uio(&c__1, (char *)&np, (ftnlen)sizeof(integer));
		    for (k = 1; k <= 3; ++k) {
			do_uio(&c__1, (char *)&pijk[k + ip * 3 - 4], (ftnlen)
				sizeof(real));
		    }
		    for (k = 1; k <= 3; ++k) {
			do_uio(&c__1, (char *)&aijk[k + ip * 3 - 4], (ftnlen)
				sizeof(real));
		    }
		    do_uio(&c__1, (char *)&yonprb, (ftnlen)sizeof(logical));
		    do_uio(&c__1, (char *)&bury, (ftnlen)sizeof(logical));
		    e_wsue();
/* save probe in yon probe arrays */
		    if (yonprb) {
/* check for overflow */
			if (ny >= 12000) {
			    error_(&c__720, &ny, &c_b36);
			}
			++ny;
			molyon[ny - 1] = imol;
			for (k = 1; k <= 3; ++k) {
			    py[k + ny * 3 - 4] = pijk[k + ip * 3 - 4];
			    ay[k + ny * 3 - 4] = aijk[k + ip * 3 - 4];
/* L4050: */
			}
		    }

/* write surface points for this probe position */
		    i__4 = np;
		    for (i__ = 1; i__ <= i__4; ++i__) {
			s_wsue(&io___174);
			do_uio(&c__1, (char *)&n1[i__ - 1], (ftnlen)sizeof(
				integer));
			do_uio(&c__1, (char *)&n2[i__ - 1], (ftnlen)sizeof(
				integer));
			do_uio(&c__1, (char *)&n3[i__ - 1], (ftnlen)sizeof(
				integer));
			for (k = 1; k <= 3; ++k) {
			    do_uio(&c__1, (char *)&s[k + i__ * 3 - 4], (
				    ftnlen)sizeof(real));
			}
			do_uio(&c__1, (char *)&area, (ftnlen)sizeof(real));
			do_uio(&c__1, (char *)&yon[i__ - 1], (ftnlen)sizeof(
				logical));
			e_wsue();
/* L4100: */
		    }
/* end of ip loop */
L4150:
		    ;
		}
/* end of concave reentrant loop */
L4200:
		;
	    }

/* saddle-shaped reentrant */
	    ishape = 2;

/* check for neither atom to be surfaces */
	    if (! (si || sj)) {
		goto L5150;
	    }

/* special check for buried tori */

/* if both atoms are marked to be surface, */
/* but neither atom has any reentrant surface so far */
/* (after triangles with all katoms have been checked) */
/* and if there is some mutual neighbor in the same molecule */
/* close enough so that the torus cannot be free, */
/* then we know that this must be a buried torus */

	    if (si && sj && ! srs[iatom - 1] && ! srs[jatom - 1] && mutual > 
		    0) {
		i__3 = nnbr;
		for (knbr = 1; knbr <= i__3; ++knbr) {
		    if (! mnbr[knbr - 1]) {
			goto L4250;
		    }
		    if (imol != molnbr[knbr - 1]) {
			goto L4250;
		    }
		    d2 = dist2_(bij, &cnbr[knbr * 3 - 3]);
/* Computing 2nd power */
		    r__1 = ernbr[knbr - 1];
/* Computing 2nd power */
		    r__2 = hij;
		    rk2 = r__1 * r__1 - r__2 * r__2;
		    if (d2 < rk2) {
			goto L5150;
		    }
L4250:
		    ;
		}
	    }
/* calculate number of rotations of probe pair, */
/* rotation angle and rotation matrix */
	    rij = ri / (ri + rp) + rj / (rj + rp);
	    avh = ((r__1 = hij - rp, dabs(r__1)) + hij * rij) / 3;
	    nrot = sqrt(d__[ptyp - 1]) * 3.141592654f * avh;
	    if (nrot < 1) {
		nrot = 1;
	    }
	    angle = 3.14159f / nrot;
/* set up rotation matrix around x-axis */
	    imatx_(h__);
	    h__[4] = cos(angle);
	    h__[8] = h__[4];
	    h__[5] = sin(angle);
	    h__[7] = -h__[5];
/* calculate matrix to rotate x-axis onto iatom-jatom axis */
	    for (k = 1; k <= 3; ++k) {
		g[k - 1] = uij[k - 1];
		g[k + 2] = q[k - 1];
		g[k + 5] = t[k - 1];
/* L4300: */
	    }
/* make the probe pair rotation matrix be about the iatom-jatom axis */
	    conj_(h__, g, ghgt);

/* arc generation */
	    for (k = 1; k <= 3; ++k) {
		pij[k - 1] = bij[k - 1] + aij[k - 1];
		vpi[k - 1] = (ci[k - 1] - pij[k - 1]) * rp / (ri + rp);
		vpj[k - 1] = (cj[k - 1] - pij[k - 1]) * rp / (rj + rp);
/* L4350: */
	    }

/* rotate circle onto iatom-jatom-probe plane */
/* and select points between probe-iatom and */
/* probe-jatom vector to form the arc */
	    narc = 1;
	    i__3 = ncirc[ptyp - 1];
	    for (i__ = 1; i__ <= i__3; ++i__) {
		if (narc > 1000) {
		    error_(&c__440, &narc, &c_b36);
		}
/* rotation */
		multv_(&circle[(i__ + ptyp * 1000) * 3 - 3003], g, &vps0[narc 
			* 3 - 3]);
/* if the vector is pointing away from the symmetry line */
/* the surface point cannot lie on the inward-facing arc */
		if (dot_(&vps0[narc * 3 - 3], aij) > 0.f) {
		    goto L4500;
		}
		cross_(vpi, &vps0[narc * 3 - 3], vector);
		if (dot_(&g[6], vector) < 0.f) {
		    goto L4500;
		}
		cross_(&vps0[narc * 3 - 3], vpj, vector);
		if (dot_(&g[6], vector) < 0.f) {
		    goto L4500;
		}

/* make arc point vectors originate with saddle circle center bij */
/* rather than probe center because they will be */
/* rotated around the iatom-jatom axis */
		for (k = 1; k <= 3; ++k) {
		    vbs0[k + (narc + 1000) * 3 - 3004] = vps0[k + narc * 3 - 
			    4] + aij[k - 1];
/* L4400: */
		}
/* invert arc through line of symmetry */
		duij = dot_(uij, &vbs0[(narc + 1000) * 3 - 3003]);
		for (k = 1; k <= 3; ++k) {
		    vbs0[k + (narc + 2000) * 3 - 3004] = -vbs0[k + (narc + 
			    1000) * 3 - 3004] + duij * 2 * uij[k - 1];
/* L4450: */
		}

/* check whether the arc point crosses the iatom-jatom axis */
/* and calculate the area associated with the point */
		ht = dot_(aij, &vbs0[(narc + 1000) * 3 - 3003]) / hij;
		ayon[narc - 1] = ht < 0.f;
		arca[narc - 1] = rp * 19.739175456199998f * dabs(ht) / (ncirc[
			ptyp - 1] * nrot);
		++narc;
L4500:
		;
	    }
	    --narc;

/* initialize power matrix to identity */
	    imatx_(pow);

/* set knbr to zero for collision and buried checks */
	    knbr = 0;
/* rotate the probe pair around the pair of atoms */
	    i__3 = nrot;
	    for (irot = 1; irot <= i__3; ++irot) {
/* multiply altitude vector by power matrix */
		multv_(aij, pow, aijp);
/* set up opposing altitude */
		for (k = 1; k <= 3; ++k) {
		    aijp[k + 2] = -aijp[k - 1];
/* L4550: */
		}
/* set up probe sphere positions */
		for (ip = 1; ip <= 2; ++ip) {
		    for (k = 1; k <= 3; ++k) {
			pijp[k + ip * 3 - 4] = bij[k - 1] + aijp[k + ip * 3 - 
				4];
/* L4600: */
		    }
/* check for collisions with neighboring atoms */
		    pair[ip - 1] = ! collid_(&pijp[ip * 3 - 3], &rp, cnbr, 
			    ernbr, mnbr, &nnbr, &c__2000, &ishape, &jnbr, &
			    knbr, molnbr, &imol, &lkf, lknbr);
/* L4650: */
		}
/* no surface generation if neither probe position is allowed */
		if (! pair[0] && ! pair[1]) {
		    goto L5050;
		}
		both = pair[0] && pair[1];
/* some reentrant surface for both atoms */
		srs[iatom - 1] = TRUE_;
		srs[jatom - 1] = TRUE_;
/* skip to bottom of middle loop if iatom and jatom */
/* are close enough and the surface point density is */
/* low enough so that the arc has no points */
		if (narc <= 0) {
		    goto L5050;
		}

/* surface generation */
		for (ip = 1; ip <= 2; ++ip) {
		    if (! pair[ip - 1]) {
			goto L5000;
		    }
/* set default value for bury */
		    bury = FALSE_;
/* don't check for probe collisions against other molecules */
/* unless we need to */
		    if (ibury > 0) {
			bury = buried_(&pijp[ip * 3 - 3], &rp, cnbr, rnbr, 
				mnbr, &nnbr, &c__2000, &ishape, &jnbr, &knbr, 
				molnbr, &imol);
		    }
/* skip if not buried and buried surface only flag set */
		    if (ibury == 1 && ! bury) {
			goto L5000;
		    }
/* determine whether probe has surface on far side of line */
		    yonprb = hij < rp && ! both;
		    np = 1;
/* the saddle-shaped reentrant surface points come from the arc */
		    i__4 = narc;
		    for (i__ = 1; i__ <= i__4; ++i__) {
/* overlapping reentrant surface removal */
/* for symmetry-related probe positions */
			if (both && ayon[i__ - 1]) {
			    goto L4850;
			}
			if (np > 1000) {
			    error_(&c__480, &np, &c_b36);
			}
/* rotate the arc from the xy plane onto the iatom-jatom-probe plane */
			multv_(&vbs0[(i__ + ip * 1000) * 3 - 3003], pow, vbs);
/* make coordinates relative to origin */
			for (k = 1; k <= 3; ++k) {
			    s[k + np * 3 - 4] = bij[k - 1] + vbs[k - 1];
/* L4700: */
			}
/* find the closest atom and set up the atom numbers for the point */
			dsi = dist_(&s[np * 3 - 3], ci) - ri;
			dsj = dist_(&s[np * 3 - 3], cj) - rj;
			if (dsi <= dsj) {
			    goto L4750;
			}
			if (! sj) {
			    goto L4850;
			}
			n1[np - 1] = jatom;
			n2[np - 1] = iatom;
			n3[np - 1] = 0;
			goto L4800;
L4750:
			if (! si) {
			    goto L4850;
			}
			n1[np - 1] = iatom;
			n2[np - 1] = jatom;
			n3[np - 1] = 0;
L4800:

/* we've got a surface point */
			yon[np - 1] = ayon[i__ - 1];
			torus[np - 1] = arca[i__ - 1];
			++np;
/* end of arc point loop */
L4850:
			;
		    }
		    --np;
		    if (np <= 0) {
			goto L5000;
		    }

/* write the molecule number, shape,number of points, */
/* probe position and the vector from the base to the probe center */
		    s_wsue(&io___198);
		    do_uio(&c__1, (char *)&imol, (ftnlen)sizeof(integer));
		    do_uio(&c__1, (char *)&ishape, (ftnlen)sizeof(integer));
		    do_uio(&c__1, (char *)&np, (ftnlen)sizeof(integer));
		    for (k = 1; k <= 3; ++k) {
			do_uio(&c__1, (char *)&pijp[k + ip * 3 - 4], (ftnlen)
				sizeof(real));
		    }
		    for (k = 1; k <= 3; ++k) {
			do_uio(&c__1, (char *)&aijp[k + ip * 3 - 4], (ftnlen)
				sizeof(real));
		    }
		    do_uio(&c__1, (char *)&yonprb, (ftnlen)sizeof(logical));
		    do_uio(&c__1, (char *)&bury, (ftnlen)sizeof(logical));
		    e_wsue();
		    if (yonprb) {
/* save probe in yon probe arrays */
/* check for overflow */
			if (ny >= 12000) {
			    error_(&c__720, &ny, &c_b36);
			}
			++ny;
			molyon[ny - 1] = imol;
			for (k = 1; k <= 3; ++k) {
			    py[k + ny * 3 - 4] = pijp[k + ip * 3 - 4];
			    ay[k + ny * 3 - 4] = aijp[k + ip * 3 - 4];
/* L4900: */
			}
		    }

/* write surface points for this probe position */
		    i__4 = np;
		    for (i__ = 1; i__ <= i__4; ++i__) {
			s_wsue(&io___199);
			do_uio(&c__1, (char *)&n1[i__ - 1], (ftnlen)sizeof(
				integer));
			do_uio(&c__1, (char *)&n2[i__ - 1], (ftnlen)sizeof(
				integer));
			do_uio(&c__1, (char *)&n3[i__ - 1], (ftnlen)sizeof(
				integer));
			for (k = 1; k <= 3; ++k) {
			    do_uio(&c__1, (char *)&s[k + i__ * 3 - 4], (
				    ftnlen)sizeof(real));
			}
			do_uio(&c__1, (char *)&torus[i__ - 1], (ftnlen)sizeof(
				real));
			do_uio(&c__1, (char *)&yon[i__ - 1], (ftnlen)sizeof(
				logical));
			e_wsue();
/* end of arc point loop */
/* L4950: */
		    }
/* end of probe pair loop */
L5000:
		    ;
		}
/* skip to here if both probe positions disallowed or no arc points */
L5050:
/* calculate new power matrix */
		cat_(pow, ghgt);
/* end of rotation loop */
/* L5100: */
	    }
/* end of neighbor loop */
L5150:
	    ;
	}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* skip to here if van der waals surface calculation */
L5200:

/* contact surface */
	ishape = 1;
/* skip atom i if marked no surface requested */
	if (! si) {
	    goto L5650;
	}

/* if we are not calculating buried surface */
/* and the probe radius is greater than zero */
/* and iatom has at least one neighbor, but no reentrant surface, */
/* then iatom must be completely inaccessible to the probe */
	if (! bury && rp > 0.f && nimol > 0 && ! srs[iatom - 1]) {
	    goto L5650;
	}
/* find the index into the atom type arrays for iatom */
	i__2 = ntype;
	for (idx = 1; idx <= i__2; ++idx) {
	    if (itype[idx - 1] == iat[iatom - 1]) {
		goto L5300;
	    }
/* L5250: */
	}
	s_stop("logic error in ms regarding atom types", (ftnlen)38);
L5300:
/* Computing 2nd power */
	r__1 = ri;
	area = r__1 * r__1 * 12.566370616f / nua[idx - 1];
/* set jnbr, knbr to zero for collision, buried checks */
	jnbr = 0;
	knbr = 0;

/* contact probe placement loop */
	i__2 = nua[idx - 1];
	for (i__ = 1; i__ <= i__2; ++i__) {
/* set up probe coordinates */
	    for (k = 1; k <= 3; ++k) {
		pipt[k - 1] = ci[k - 1] + eva[k + (i__ + idx * 1000) * 3 - 
			3004];
/* L5350: */
	    }
/* check for collision with neighboring atoms */
	    if (collid_(pipt, &rp, cnbr, ernbr, mnbr, &nnbr, &c__2000, &
		    ishape, &jnbr, &knbr, molnbr, &imol, &lkf, lknbr)) {
		goto L5600;
	    }
/* go write it out if we don't care about buried surface */
	    if (ibury == 0) {
		ib = 0;
		goto L5400;
	    }
	    bury = buried_(pipt, &rp, cnbr, rnbr, mnbr, &nnbr, &c__2000, &
		    ishape, &jnbr, &knbr, molnbr, &imol);
	    if (ibury == 1 && ! bury) {
		goto L5600;
	    }
	    if (bury) {
		ib = 1;
	    } else {
		ib = 0;
	    }

L5400:
/* increment surface point counter for convex surface */
	    ++nshape[0];
/* add surface point area to contact area */
	    areac += area;
	    iatnum[0] = iatom;
	    iatnum[1] = 0;
	    iatnum[2] = 0;
	    for (k = 1; k <= 3; ++k) {
		outco[k - 1] = ci[k - 1] + ri * ua[k + (i__ + idx * 1000) * 3 
			- 3004];
		outvec[k - 1] = ua[k + (i__ + idx * 1000) * 3 - 3004];
/* L5450: */
	    }
/* four different output formats */
	    if (! binary && ! short__) {
		s_wsfe(&io___205);
		do_fio(&c__3, (char *)&iatnum[0], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ishape, (ftnlen)sizeof(integer));
		do_fio(&c__3, (char *)&outco[0], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&area, (ftnlen)sizeof(real));
		do_fio(&c__3, (char *)&outvec[0], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		e_wsfe();
	    } else if (binary && ! short__) {
		s_wsue(&io___206);
		do_uio(&c__3, (char *)&iatnum[0], (ftnlen)sizeof(integer));
		do_uio(&c__1, (char *)&ishape, (ftnlen)sizeof(integer));
		do_uio(&c__3, (char *)&outco[0], (ftnlen)sizeof(real));
		do_uio(&c__1, (char *)&area, (ftnlen)sizeof(real));
		do_uio(&c__3, (char *)&outvec[0], (ftnlen)sizeof(real));
		do_uio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		e_wsue();
	    } else if (! binary && short__) {
/* 	         write (7,5550) iatnum,ishape,outco,ib */
/*   ---dac change to write records for DeMon dummy atoms: */
		s_wsfe(&io___207);
		do_fio(&c__1, "H    ", (ftnlen)5);
		do_fio(&c__3, (char *)&outco[0], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c_b36, (ftnlen)sizeof(real));
		e_wsfe();
/* L5550: */
	    } else if (binary && short__) {
		s_wsue(&io___208);
		do_uio(&c__3, (char *)&iatnum[0], (ftnlen)sizeof(integer));
		do_uio(&c__1, (char *)&ishape, (ftnlen)sizeof(integer));
		do_uio(&c__3, (char *)&outco[0], (ftnlen)sizeof(real));
		do_uio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		e_wsue();
	    }

/* end of nua loop */
L5600:
	    ;
	}
/* end of iatom loop */
L5650:
	;
    }

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* write out messages */

/* for array dimensioning */
    s_wsfe(&io___209);
    do_fio(&c__1, (char *)&maxnb, (ftnlen)sizeof(integer));
    e_wsfe();

/* close contact and before files */
    cl__1.cerr = 0;
    cl__1.cunit = 7;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* if van der waals surface we are finished */
    if (rp == 0.f) {
	goto L8000;
    }
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);

/* orsr    orsr    orsr    orsr    orsr    orsr    orsr    orsr    orsr */


/* overlapping reentrant surface removal */
/* for non-symmetry-related probes */
/* probe diameter */
    dp = rp * 2;
/* diameter squared */
/* Computing 2nd power */
    r__1 = dp;
    dp2 = r__1 * r__1;
/* radius squared */
/* Computing 2nd power */
    r__1 = rp;
    rp2 = r__1 * r__1;
/* width for cubing algorithm */
/* 	width = dp */

/*     set up cube arrays */
/*     first the integer coordinate arrays */
    i__1 = ny;
    for (iy = 1; iy <= i__1; ++iy) {
	for (k = 1; k <= 3; ++k) {
	    ico[k + iy * 3 - 4] = (py[k + iy * 3 - 4] - comin[k - 1] - radmax 
		    - rp) / width + 1;
	    if (ico[k + iy * 3 - 4] < 1) {
		ico[k + iy * 3 - 4] = 1;
	    }
	    if (ico[k + iy * 3 - 4] > 40) {
		s_stop("cube coordinate too large", (ftnlen)25);
	    }
/* L5800: */
	}
/* L5850: */
    }

/* initialize head pointer array */
    for (k = 1; k <= 40; ++k) {
	for (j = 1; j <= 40; ++j) {
	    for (i__ = 1; i__ <= 40; ++i__) {
		icube[i__ + (j + k * 40) * 40 - 1641] = 0;
/* L5900: */
	    }
/* L5950: */
	}
/* L6000: */
    }

/* initialize linked list pointers */
    for (iy = 1; iy <= 12000; ++iy) {
	icuptr[iy - 1] = 0;
/* L6050: */
    }

/* set up head and later pointers for each yon probe */
    i__1 = ny;
    for (iy = 1; iy <= i__1; ++iy) {
/* skip atoms with surface request numbers of zero */
	i__ = ico[iy * 3 - 3];
	j = ico[iy * 3 - 2];
	k = ico[iy * 3 - 1];
	if (icube[i__ + (j + k * 40) * 40 - 1641] <= 0) {
/*     first atom in this cube */
	    icube[i__ + (j + k * 40) * 40 - 1641] = iy;
	} else {
/*     add to end of linked list */
	    iptr = icube[i__ + (j + k * 40) * 40 - 1641];
L6100:
	    if (icuptr[iptr - 1] <= 0) {
		goto L6150;
	    }
	    iptr = icuptr[iptr - 1];
	    goto L6100;
L6150:
	    icuptr[iptr - 1] = iy;
	}
/* L6200: */
    }

/* reopen before file for reading */
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 6;
    o__1.ofnm = "before";
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = "unformatted";
    o__1.oblnk = 0;
    f_open(&o__1);
    al__1.aerr = 0;
    al__1.aunit = 4;
    f_rew(&al__1);

/* first pass */
/* gather victim probes */
    nv = 1;
/* no victim probes if no yon probes */
    if (ny <= 0) {
	goto L6950;
    }
    al__1.aerr = 0;
    al__1.aunit = 4;
    f_rew(&al__1);
/* initialize victim hashing array */
    for (j = 1; j <= 6000; ++j) {
	ivic[j - 1] = 0;
/* L6250: */
    }
/* initialize index of free slot last used */
    ifrlst = 0;
/* initialize probe record number */
    i__ = 1;
L6300:
/* check for victim overflow */
    if (nv > 6000) {
	error_(&c__760, &nv, &c_b36);
    }
/* read reentrant probe and points */
    i__1 = s_rsue(&io___217);
    if (i__1 != 0) {
	goto L6950;
    }
    i__1 = do_uio(&c__1, (char *)&molvic[nv - 1], (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L6950;
    }
    i__1 = do_uio(&c__1, (char *)&ishape, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L6950;
    }
    i__1 = do_uio(&c__1, (char *)&np, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L6950;
    }
    for (k = 1; k <= 3; ++k) {
	i__1 = do_uio(&c__1, (char *)&pv[k + nv * 3 - 4], (ftnlen)sizeof(real)
		);
	if (i__1 != 0) {
	    goto L6950;
	}
    }
    for (k = 1; k <= 3; ++k) {
	i__1 = do_uio(&c__1, (char *)&av[k + nv * 3 - 4], (ftnlen)sizeof(real)
		);
	if (i__1 != 0) {
	    goto L6950;
	}
    }
    i__1 = do_uio(&c__1, (char *)&yonprb, (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L6950;
    }
    i__1 = do_uio(&c__1, (char *)&bury, (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L6950;
    }
    i__1 = e_rsue();
    if (i__1 != 0) {
	goto L6950;
    }
    i__1 = np;
    for (j = 1; j <= i__1; ++j) {
	s_rsue(&io___221);
	do_uio(&c__1, (char *)&n1[0], (ftnlen)sizeof(integer));
	do_uio(&c__1, (char *)&n2[0], (ftnlen)sizeof(integer));
	do_uio(&c__1, (char *)&n3[0], (ftnlen)sizeof(integer));
	for (k = 1; k <= 3; ++k) {
	    do_uio(&c__1, (char *)&s[k - 1], (ftnlen)sizeof(real));
	}
	do_uio(&c__1, (char *)&area, (ftnlen)sizeof(real));
	do_uio(&c__1, (char *)&yon[0], (ftnlen)sizeof(logical));
	e_rsue();
/* L6350: */
    }
    if (yonprb) {
	goto L6900;
    }
/* check if probe too far from symmetry element for possible overlap */
    if (anorm_(&av[nv * 3 - 3]) > dp) {
	goto L6900;
    }

/* look for overlap with any yon probe in the same molecule */
/* use cubing algorithm to save time */

/* calculate which cube this probe lies in */
    ici = (pv[nv * 3 - 3] - comin[0] - radmax - rp) / width + 1;
    if (ici < 1) {
	ici = 1;
    }
    if (ici > 40) {
	s_stop("cube coordinate too large", (ftnlen)25);
    }
    icj = (pv[nv * 3 - 2] - comin[1] - radmax - rp) / width + 1;
    if (icj < 1) {
	icj = 1;
    }
    if (icj > 40) {
	s_stop("cube coordinate too large", (ftnlen)25);
    }
    ick = (pv[nv * 3 - 1] - comin[2] - radmax - rp) / width + 1;
    if (ick < 1) {
	ick = 1;
    }
    if (ick > 40) {
	s_stop("cube coordinate too large", (ftnlen)25);
    }
/* check for overlap with probes in adjoining cubes */
    i__1 = ick + 1;
    for (jck = ick - 1; jck <= i__1; ++jck) {
	if (jck < 1 || jck > 40) {
	    goto L6850;
	}
	i__2 = icj + 1;
	for (jcj = icj - 1; jcj <= i__2; ++jcj) {
	    if (jcj < 1 || jcj > 40) {
		goto L6800;
	    }
	    i__3 = ici + 1;
	    for (jci = ici - 1; jci <= i__3; ++jci) {
		if (jci < 1 || jci > 40) {
		    goto L6750;
		}
		jp = icube[jci + (jcj + jck * 40) * 40 - 1641];


L6400:
		if (jp <= 0) {
		    goto L6700;
		}
		if (molyon[jp - 1] != molvic[nv - 1]) {
		    goto L6650;
		}
		x = (r__1 = py[jp * 3 - 3] - pv[nv * 3 - 3], dabs(r__1));
		if (x >= dp) {
		    goto L6650;
		}
		y = (r__1 = py[jp * 3 - 2] - pv[nv * 3 - 2], dabs(r__1));
		if (y >= dp) {
		    goto L6650;
		}
		z__ = (r__1 = py[jp * 3 - 1] - pv[nv * 3 - 1], dabs(r__1));
		if (z__ >= dp) {
		    goto L6650;
		}
/* Computing 2nd power */
		r__1 = x;
/* Computing 2nd power */
		r__2 = y;
/* Computing 2nd power */
		r__3 = z__;
		d2 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
		if (d2 >= dp2) {
		    goto L6650;
		}
/* check that probes face each other */
		if (dot_(&ay[jp * 3 - 3], &av[nv * 3 - 3]) >= 0.f) {
		    goto L6650;
		}
/* new victim probe */
/* put into hashing table */
		ihash = i__ % 6000 + 1;
		if (ivic[ihash - 1] == 0) {
/* empty slot */
		    ivic[ihash - 1] = i__;
		    ivicp[ihash - 1] = 0;
		} else {
		    iprev = ihash;
		    iptr = ivicp[ihash - 1];
L6450:
/* check for end of linked list */
		    if (iptr == 0) {
			goto L6500;
		    }
		    iprev = iptr;
		    iptr = ivicp[iptr - 1];
		    goto L6450;
L6500:
/* look for a free slot */
		    for (ifree = ifrlst + 1; ifree <= 6000; ++ifree) {
			if (ivic[ifree - 1] == 0) {
			    goto L6600;
			}
/* L6550: */
		    }
		    s_stop("victim oveflow", (ftnlen)14);
L6600:
/* store record number in free slot */
		    ivic[ifree - 1] = i__;
		    ivicp[iprev - 1] = ifree;
		    ivicp[ifree - 1] = 0;
/* new index to last free slot used */
		    ifrlst = ifree;
		}
		++nv;
/* one overlap makes this probe a victim */
/* we don't need to check any more */
		goto L6900;
L6650:
		jp = icuptr[jp - 1];
		goto L6400;
L6700:
L6750:
		;
	    }
L6800:
	    ;
	}
L6850:
	;
    }
/* end of yon probe loop */
/* skip to here if finished with hunt for overlapping probes */
L6900:
    ++i__;
    goto L6300;
/* skip to here if there are no yon probes and hence no victims */
L6950:
    --nv;

/* open reentrant file for writing */
    if (binary) {
	o__1.oerr = 0;
	o__1.ounit = 8;
	o__1.ofnmlen = 9;
	o__1.ofnm = "reentrant";
	o__1.orl = 0;
	o__1.osta = "unknown";
	o__1.oacc = 0;
	o__1.ofm = "unformatted";
	o__1.oblnk = 0;
	f_open(&o__1);
    } else {
	o__1.oerr = 0;
	o__1.ounit = 8;
	o__1.ofnmlen = 9;
	o__1.ofnm = "reentrant";
	o__1.orl = 0;
	o__1.osta = "unknown";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }
    al__1.aerr = 0;
    al__1.aunit = 8;
    f_rew(&al__1);

/* second pass */
/* read, check and write surface points */
    al__1.aerr = 0;
    al__1.aunit = 4;
    f_rew(&al__1);
    i__ = 1;
L7000:
    i__1 = s_rsue(&io___230);
    if (i__1 != 0) {
	goto L7850;
    }
    i__1 = do_uio(&c__1, (char *)&imol, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L7850;
    }
    i__1 = do_uio(&c__1, (char *)&ishape, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L7850;
    }
    i__1 = do_uio(&c__1, (char *)&np, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L7850;
    }
    for (k = 1; k <= 3; ++k) {
	i__1 = do_uio(&c__1, (char *)&p[k - 1], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L7850;
	}
    }
    for (k = 1; k <= 3; ++k) {
	i__1 = do_uio(&c__1, (char *)&a[k - 1], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L7850;
	}
    }
    i__1 = do_uio(&c__1, (char *)&yonprb, (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L7850;
    }
    i__1 = do_uio(&c__1, (char *)&bury, (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L7850;
    }
    i__1 = e_rsue();
    if (i__1 != 0) {
	goto L7850;
    }
/* no points can be eaten if this probe is neither yon nor a victim */
    neat = 0;
    nyeat = 0;
    if (ny <= 0) {
	goto L7450;
    }
    ipt = 0;
/* determine if probe is a yon or victim probe */
    if (! yonprb) {
	goto L7050;
    }
/* we've got a yon probe here */
    ipt = 2;
    goto L7200;
L7050:
    if (nv <= 0) {
	goto L7450;
    }
/* hash into table of victim probes */
    iptr = i__ % 6000 + 1;
L7100:
    if (iptr == 0) {
	goto L7200;
    }
    if (ivic[iptr - 1] == 0) {
	goto L7200;
    }
    if (ivic[iptr - 1] == i__) {
	goto L7150;
    }
    iptr = ivicp[iptr - 1];
    goto L7100;
L7150:
/* we've got a victim */
    ipt = 1;
L7200:

    if (ipt <= 0) {
	goto L7450;
    }
/* check this victim or yon probe against all yon probes */
    i__1 = ny;
    for (j = 1; j <= i__1; ++j) {
	if (imol != molyon[j - 1]) {
	    goto L7300;
	}
	if (dist2_(p, &py[j * 3 - 3]) >= dp2) {
	    goto L7300;
	}
	if (dot_(a, &ay[j * 3 - 3]) >= 0.f) {
	    goto L7300;
	}
/* this yon probe could eat some of the probe's points */
	++neat;
	++nyeat;
	if (neat > 1000) {
	    error_(&c__830, &neat, &c_b36);
	}
	for (k = 1; k <= 3; ++k) {
	    eat[k + neat * 3 - 4] = py[k + j * 3 - 4];
/* L7250: */
	}
/* end of yon probe loop */
L7300:
	;
    }

/* only yon probes can have their points eaten by victims */
    if (ipt <= 1) {
	goto L7450;
    }
/* check this yon probe against all victim probes */
    i__1 = nv;
    for (j = 1; j <= i__1; ++j) {
	if (imol != molvic[j - 1]) {
	    goto L7400;
	}
	if (dist2_(p, &pv[j * 3 - 3]) >= dp2) {
	    goto L7400;
	}
	if (dot_(a, &av[j * 3 - 3]) >= 0.f) {
	    goto L7400;
	}
/* this victim probe could eat some of the probe's points */
	++neat;
	if (neat > 1000) {
	    error_(&c__850, &neat, &c_b36);
	}
	for (k = 1; k <= 3; ++k) {
	    eat[k + neat * 3 - 4] = pv[k + j * 3 - 4];
/* L7350: */
	}
/* end of victim probe loop */
L7400:
	;
    }

/* skip to here if victim or both probe overlap checks omitted */
L7450:

/* read the surface points belonging to the probe */
    i__1 = np;
    for (j = 1; j <= i__1; ++j) {
	s_rsue(&io___237);
	do_uio(&c__1, (char *)&n1[0], (ftnlen)sizeof(integer));
	do_uio(&c__1, (char *)&n2[0], (ftnlen)sizeof(integer));
	do_uio(&c__1, (char *)&n3[0], (ftnlen)sizeof(integer));
	for (k = 1; k <= 3; ++k) {
	    do_uio(&c__1, (char *)&s[k - 1], (ftnlen)sizeof(real));
	}
	do_uio(&c__1, (char *)&area, (ftnlen)sizeof(real));
	do_uio(&c__1, (char *)&yon[0], (ftnlen)sizeof(logical));
	e_rsue();
	if (neat <= 0) {
	    goto L7550;
	}
/* check surface point against all eaters of this probe */
	i__2 = neat;
	for (k = 1; k <= i__2; ++k) {
/* victim probes cannot eat non-yon points of yon probes */
	    if (yonprb && ! yon[0] && k > nyeat) {
		goto L7500;
	    }
	    if (dist2_(&eat[k * 3 - 3], s) < rp2) {
		goto L7700;
	    }
L7500:
	    ;
	}
/* skip to here if no overlapping probes could eat this point */
L7550:
	for (k = 1; k <= 3; ++k) {
	    outvec[k - 1] = (p[k - 1] - s[k - 1]) / rp;
/* L7600: */
	}
/* reentrant surface point */
	++nshape[ishape - 1];
	arear += area;
/* mark whether buried */
	ib = 0;
	if (ibury > 0 && bury) {
	    ib = 1;
	}
/* four possible output formats */
	iatnum[0] = n1[0];
	iatnum[1] = n2[0];
	iatnum[2] = n3[0];
	for (k = 1; k <= 3; ++k) {
	    outco[k - 1] = s[k - 1];
/* L7650: */
	}
	if (! binary && ! short__) {
	    s_wsfe(&io___238);
	    do_fio(&c__3, (char *)&iatnum[0], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ishape, (ftnlen)sizeof(integer));
	    do_fio(&c__3, (char *)&outco[0], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&area, (ftnlen)sizeof(real));
	    do_fio(&c__3, (char *)&outvec[0], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else if (binary && ! short__) {
	    s_wsue(&io___239);
	    do_uio(&c__3, (char *)&iatnum[0], (ftnlen)sizeof(integer));
	    do_uio(&c__1, (char *)&ishape, (ftnlen)sizeof(integer));
	    do_uio(&c__3, (char *)&outco[0], (ftnlen)sizeof(real));
	    do_uio(&c__1, (char *)&area, (ftnlen)sizeof(real));
	    do_uio(&c__3, (char *)&outvec[0], (ftnlen)sizeof(real));
	    do_uio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    e_wsue();
	} else if (! binary && short__) {
	    s_wsfe(&io___240);
	    do_fio(&c__3, (char *)&iatnum[0], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ishape, (ftnlen)sizeof(integer));
	    do_fio(&c__3, (char *)&outco[0], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else if (binary && short__) {
	    s_wsue(&io___241);
	    do_uio(&c__3, (char *)&iatnum[0], (ftnlen)sizeof(integer));
	    do_uio(&c__1, (char *)&ishape, (ftnlen)sizeof(integer));
	    do_uio(&c__3, (char *)&outco[0], (ftnlen)sizeof(real));
	    do_uio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    e_wsue();
	}
	goto L7750;
L7700:
	++nlost[ishape - 1];
/* end of np loop */
L7750:
	;
    }
/* end of i loop */
/* L7800: */
    ++i__;
    goto L7000;
L7850:
/* close files */
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
    cl__1.cerr = 0;
    cl__1.cunit = 8;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* messages */
    s_wsfe(&io___242);
    do_fio(&c__1, (char *)&ny, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nv, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___243);
    do_fio(&c__1, (char *)&nlost[1], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nlost[2], (ftnlen)sizeof(integer));
    e_wsfe();
L8000:
/* write out how many points */
    s_wsfe(&io___244);
    do_fio(&c__1, (char *)&nshape[0], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nshape[1], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nshape[2], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___245);
    i__1 = nshape[0] + nshape[1] + nshape[2];
    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___246);
    do_fio(&c__1, (char *)&areac, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&arear, (ftnlen)sizeof(real));
    r__1 = areac + arear;
    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
    e_wsfe();
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */



/* subroutines and functions */


/* general vector and matrix routines */

doublereal dist_(real *a, real *b)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

    /* Builtin functions */
    double sqrt(doublereal);

/* distance between a and b */
    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
/* Computing 2nd power */
    r__1 = a[1] - b[1];
/* Computing 2nd power */
    r__2 = a[2] - b[2];
/* Computing 2nd power */
    r__3 = a[3] - b[3];
    ret_val = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    return ret_val;
} /* dist_ */


doublereal dist2_(real *a, real *b)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

/* distance between a and b squared */
    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
/* Computing 2nd power */
    r__1 = a[1] - b[1];
/* Computing 2nd power */
    r__2 = a[2] - b[2];
/* Computing 2nd power */
    r__3 = a[3] - b[3];
    ret_val = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
    return ret_val;
} /* dist2_ */


doublereal anorm_(real *a)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

    /* Builtin functions */
    double sqrt(doublereal);

/* norm of a */
    /* Parameter adjustments */
    --a;

    /* Function Body */
/* Computing 2nd power */
    r__1 = a[1];
/* Computing 2nd power */
    r__2 = a[2];
/* Computing 2nd power */
    r__3 = a[3];
    ret_val = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    return ret_val;
} /* anorm_ */


doublereal dot_(real *a, real *b)
{
    /* System generated locals */
    real ret_val;

/* dot product */
    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    ret_val = a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
    return ret_val;
} /* dot_ */


/* Subroutine */ int cross_(real *a, real *b, real *c__)
{
/* cross product */
    /* Parameter adjustments */
    --c__;
    --b;
    --a;

    /* Function Body */
    c__[1] = a[2] * b[3] - a[3] * b[2];
    c__[2] = a[3] * b[1] - a[1] * b[3];
    c__[3] = a[1] * b[2] - a[2] * b[1];
    return 0;
} /* cross_ */


/* Subroutine */ int multv_(real *v, real *a, real *w)
{
    static integer i__;

/* multiply v by a giving w */
    /* Parameter adjustments */
    --w;
    a -= 4;
    --v;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	w[i__] = a[i__ + 3] * v[1] + a[i__ + 6] * v[2] + a[i__ + 9] * v[3];
/* L50: */
    }
    return 0;
} /* multv_ */


/* Subroutine */ int vnorm_(real *a, real *b)
{
    static integer k;
    static real v;
    extern doublereal anorm_(real *);

/* normalize a giving b */
    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    v = anorm_(&a[1]);
    for (k = 1; k <= 3; ++k) {
	b[k] = a[k] / v;
/* L50: */
    }
    return 0;
} /* vnorm_ */


/* Subroutine */ int vperp_(real *a, real *b)
{
    /* System generated locals */
    real r__1, r__2, r__3;

    /* Local variables */
    static integer k, m;
    static real p[3], dt, small;
    extern /* Subroutine */ int vnorm_(real *, real *);

/* return b perpendicular to a */
/* find smallest component */
    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    small = 1e4f;
    m = 0;
    for (k = 1; k <= 3; ++k) {
	if ((r__1 = a[k], dabs(r__1)) >= small) {
	    goto L50;
	}
	small = (r__1 = a[k], dabs(r__1));
	m = k;
L50:
	;
    }
    for (k = 1; k <= 3; ++k) {
	b[k] = 0.f;
	if (k == m) {
	    b[k] = 1.f;
	}
/* L100: */
    }
/* take projection along a */
/* Computing 2nd power */
    r__1 = a[1];
/* Computing 2nd power */
    r__2 = a[2];
/* Computing 2nd power */
    r__3 = a[3];
    dt = a[m] / (r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    for (k = 1; k <= 3; ++k) {
	p[k - 1] = dt * a[k];
/* subtract projection from b */
	b[k] -= p[k - 1];
/* L150: */
    }
/* renormalize b */
    vnorm_(&b[1], &b[1]);
    return 0;
} /* vperp_ */


/* Subroutine */ int cat_(real *a, real *b)
{
    static integer i__, j;
    static real temp[9]	/* was [3][3] */;

/* concatenate matrix b into matrix a */
    /* Parameter adjustments */
    b -= 4;
    a -= 4;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    temp[i__ + j * 3 - 4] = a[i__ + 3] * b[j * 3 + 1] + a[i__ + 6] * 
		    b[j * 3 + 2] + a[i__ + 9] * b[j * 3 + 3];
/* L50: */
	}
/* L100: */
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    a[i__ + j * 3] = temp[i__ + j * 3 - 4];
/* L150: */
	}
/* L200: */
    }
    return 0;
} /* cat_ */


/* Subroutine */ int conj_(real *h__, real *g, real *ghgt)
{
    static integer k, l;
    static real gt[9]	/* was [3][3] */;
    extern /* Subroutine */ int cat_(real *, real *), imatx_(real *);

/* conjugate matrix g with matrix h giving ghgt */
/* initialize ghgt matrix to identity */
/* concatenate g h gt */
    /* Parameter adjustments */
    ghgt -= 4;
    g -= 4;
    h__ -= 4;

    /* Function Body */
    imatx_(&ghgt[4]);
    cat_(&ghgt[4], &g[4]);
    cat_(&ghgt[4], &h__[4]);
/* calculate gt */
    for (k = 1; k <= 3; ++k) {
	for (l = 1; l <= 3; ++l) {
	    gt[k + l * 3 - 4] = g[l + k * 3];
/* L50: */
	}
/* L100: */
    }
    cat_(&ghgt[4], gt);
    return 0;
} /* conj_ */


/* Subroutine */ int imatx_(real *a)
{
    static integer i__, j;

/* load identity matrix */
    /* Parameter adjustments */
    a -= 4;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    a[i__ + j * 3] = 0.f;
/* L50: */
	}
	a[i__ + i__ * 3] = 1.f;
/* L100: */
    }
    return 0;
} /* imatx_ */


doublereal det_(real *a, real *b, real *c__)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real ab[3];
    extern doublereal dot_(real *, real *);
    extern /* Subroutine */ int cross_(real *, real *, real *);

/* return triple product of the three vectors */
    /* Parameter adjustments */
    --c__;
    --b;
    --a;

    /* Function Body */
    cross_(&a[1], &b[1], ab);
    ret_val = dot_(ab, &c__[1]);
    return ret_val;
} /* det_ */


/* geometric routines */


logical collid_(real *p, real *rp, real *cnbr, real *ernbr, logical *mnbr, 
	integer *nnbr, integer *maxnbr, integer *ishape, integer *jnbr, 
	integer *knbr, integer *molnbr, integer *imol, integer *lkf, integer *
	lknbr)
{
    /* System generated locals */
    real r__1, r__2, r__3;
    logical ret_val;

    /* Local variables */
    static integer i__;
    static real dd2, sr2, vect1, vect2, vect3;

/* collision check of probe with neighboring atoms */
/* belonging to the same molecule */
/* MAS */
/* 	real	colpre(3) */
/* 	real 	radpre,r2pre */
/* MAS */

/* MAS check whether probe is too close to any neighbor */

/* 	check neighbor from previous collision first */
    /* Parameter adjustments */
    --p;
    --lknbr;
    --molnbr;
    --mnbr;
    --ernbr;
    cnbr -= 4;

    /* Function Body */
    vect1 = (r__1 = p[1] - mpck_1.colpre[0], dabs(r__1));
    vect2 = (r__1 = p[2] - mpck_1.colpre[1], dabs(r__1));
    vect3 = (r__1 = p[3] - mpck_1.colpre[2], dabs(r__1));
/* Computing 2nd power */
    r__1 = vect1;
/* Computing 2nd power */
    r__2 = vect2;
/* Computing 2nd power */
    r__3 = vect3;
    dd2 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
    if (dd2 < mpck_1.r2pre) {
	ret_val = TRUE_;
	return ret_val;
    }

/* MAS */

    i__ = *lkf;
    goto L100;
L50:
    i__ = lknbr[i__];
L100:
    if (i__ == 0) {
	goto L150;
    }
    vect1 = (r__1 = p[1] - cnbr[i__ * 3 + 1], dabs(r__1));
    if (vect1 >= ernbr[i__]) {
	goto L50;
    }
    vect2 = (r__1 = p[2] - cnbr[i__ * 3 + 2], dabs(r__1));
    if (vect2 >= ernbr[i__]) {
	goto L50;
    }
    vect3 = (r__1 = p[3] - cnbr[i__ * 3 + 3], dabs(r__1));
    if (vect3 >= ernbr[i__]) {
	goto L50;
    }
    if (i__ == *jnbr || i__ == *knbr) {
	goto L50;
    }
/* Computing 2nd power */
    r__1 = ernbr[i__];
    sr2 = r__1 * r__1;
/* Computing 2nd power */
    r__1 = vect1;
/* Computing 2nd power */
    r__2 = vect2;
/* Computing 2nd power */
    r__3 = vect3;
    dd2 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
    if (dd2 >= sr2) {
	goto L50;
    }
    ret_val = TRUE_;

/* MAS	found a collision */

    mpck_1.colpre[0] = cnbr[i__ * 3 + 1];
    mpck_1.colpre[1] = cnbr[i__ * 3 + 2];
    mpck_1.colpre[2] = cnbr[i__ * 3 + 3];
/* 	radpre = ernbr(i) */
    mpck_1.r2pre = sr2;
/* MAS */
    return ret_val;
L150:
    ret_val = FALSE_;
    return ret_val;
} /* collid_ */


logical buried_(real *p, real *rp, real *cnbr, real *rnbr, logical *mnbr, 
	integer *nnbr, integer *maxnbr, integer *ishape, integer *jnbr, 
	integer *knbr, integer *molnbr, integer *imol)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3;
    logical ret_val;

    /* Local variables */
    static integer i__;
    static real dd2, sr2, vect1, vect2, vect3, sumrad;

/* collision check of probe with neighboring atoms */
/* belonging to a different molecule */

    /* Parameter adjustments */
    --p;
    --molnbr;
    --mnbr;
    --rnbr;
    cnbr -= 4;

    /* Function Body */
    if (*nnbr <= 0) {
	goto L100;
    }
/* check whether probe is too close to any neighbor */
    i__1 = *nnbr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*imol == molnbr[i__]) {
	    goto L50;
	}
	if (*ishape > 1 && i__ == *jnbr) {
	    goto L50;
	}
	if (*ishape == 3 && (i__ == *knbr || ! mnbr[i__])) {
	    goto L50;
	}
	sumrad = *rp + rnbr[i__];
	vect1 = (r__1 = p[1] - cnbr[i__ * 3 + 1], dabs(r__1));
	if (vect1 >= sumrad) {
	    goto L50;
	}
	vect2 = (r__1 = p[2] - cnbr[i__ * 3 + 2], dabs(r__1));
	if (vect2 >= sumrad) {
	    goto L50;
	}
	vect3 = (r__1 = p[3] - cnbr[i__ * 3 + 3], dabs(r__1));
	if (vect3 >= sumrad) {
	    goto L50;
	}
/* Computing 2nd power */
	r__1 = sumrad;
	sr2 = r__1 * r__1;
/* Computing 2nd power */
	r__1 = vect1;
/* Computing 2nd power */
	r__2 = vect2;
/* Computing 2nd power */
	r__3 = vect3;
	dd2 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
	if (dd2 < sr2) {
	    goto L150;
	}
L50:
	;
    }
L100:
    ret_val = FALSE_;
    goto L200;
L150:
    ret_val = TRUE_;
L200:
    return ret_val;
} /* buried_ */


/* Subroutine */ int genun_(real *u, integer *n)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j;
    static real x, y, z__, fi, fj;
    static integer nu;
    static real xy;
    static integer nhor, nvert, nequat;

/* generate unit vectors over sphere */
    /* Parameter adjustments */
    u -= 4;

    /* Function Body */
    nequat = sqrt(*n * 3.14159f);
    nvert = nequat * .5f;
    if (nvert < 1) {
	nvert = 1;
    }
    nu = 0;
    i__1 = nvert;
    for (i__ = 0; i__ <= i__1; ++i__) {
	fi = i__ * 3.14159f / nvert;
	z__ = cos(fi);
	xy = sin(fi);
	nhor = nequat * xy;
	if (nhor < 1) {
	    nhor = 1;
	}
	i__2 = nhor - 1;
	for (j = 0; j <= i__2; ++j) {
	    fj = j * 6.2831799999999998f / nhor;
	    x = cos(fj) * xy;
	    y = sin(fj) * xy;
	    if (nu >= *n) {
		goto L150;
	    }
	    ++nu;
	    u[nu * 3 + 1] = x;
	    u[nu * 3 + 2] = y;
	    u[nu * 3 + 3] = z__;
/* L50: */
	}
/* L100: */
    }
L150:
    *n = nu;
    return 0;
} /* genun_ */


/* error message subroutine */

/* Subroutine */ int error_(integer *number, integer *int__, real *float__)
{
    /* Initialized data */

    static integer list[15] = { 120,127,130,140,150,160,170,210,320,440,480,
	    720,760,830,850 };

    /* Format strings */
    static char fmt_100[] = "(1x,\002error of unidentifiable type\002)";
    static char fmt_250[] = "(1x,\002error\002,i5,2x,\002negative probe radi"
	    "us:\002,f10.5)";
    static char fmt_350[] = "(1x,\002error\002,i5,2x,\002bad buried surface "
	    "flag:\002,i5)";
    static char fmt_450[] = "(1x,\002error\002,i5,2x,\002too few or too many"
	    " atom types:\002,i5)";
    static char fmt_550[] = "(1x,\002error\002,i5,2x,\002negative atom radiu"
	    "s:\002,f10.5,\002 atom\002,i5)";
    static char fmt_650[] = "(1x,\002error\002,i5,2x,\002too many atoms:\002"
	    ",i5)";
    static char fmt_750[] = "(1x,\002error\002,i5,2x,\002invalid surface req"
	    "uest number for atom:\002,i5)";
    static char fmt_850[] = "(1x,\002error\002,i5,2x,\002invalid atom type f"
	    "or atom:\002,i5)";
    static char fmt_950[] = "(1x,\002error\002,i5,2x,\002too many neighbors"
	    ":\002,i5)";
    static char fmt_1050[] = "(1x,\002error\002,i5,2x,\002too many points fo"
	    "r reentrant probe:\002,i5)";
    static char fmt_1150[] = "(1x,\002error\002,i5,2x,\002too many points fo"
	    "r arc:\002,i5)";
    static char fmt_1250[] = "(1x,\002error\002,i5,2x,\002too many points fo"
	    "r reentrant probe:\002,i5)";
    static char fmt_1350[] = "(1x,\002error\002,i5,2x,\002too many yon probe"
	    "s:\002,i5)";
    static char fmt_1450[] = "(1x,\002error\002,i5,2x,\002too many victim pr"
	    "obes:\002,i5)";
    static char fmt_1550[] = "(1x,\002error\002,i5,2x,\002too many eaters"
	    ":\002,i5)";
    static char fmt_1650[] = "(1x,\002error\002,i5,2x,\002too many eaters"
	    ":\002,i5)";

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__;

    /* Fortran I/O blocks */
    static cilist io___291 = { 0, 6, 0, fmt_100, 0 };
    static cilist io___292 = { 0, 6, 0, fmt_250, 0 };
    static cilist io___293 = { 0, 6, 0, fmt_350, 0 };
    static cilist io___294 = { 0, 6, 0, fmt_450, 0 };
    static cilist io___295 = { 0, 6, 0, fmt_550, 0 };
    static cilist io___296 = { 0, 6, 0, fmt_650, 0 };
    static cilist io___297 = { 0, 6, 0, fmt_750, 0 };
    static cilist io___298 = { 0, 6, 0, fmt_850, 0 };
    static cilist io___299 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___300 = { 0, 6, 0, fmt_1050, 0 };
    static cilist io___301 = { 0, 6, 0, fmt_1150, 0 };
    static cilist io___302 = { 0, 6, 0, fmt_1250, 0 };
    static cilist io___303 = { 0, 6, 0, fmt_1350, 0 };
    static cilist io___304 = { 0, 6, 0, fmt_1450, 0 };
    static cilist io___305 = { 0, 6, 0, fmt_1550, 0 };
    static cilist io___306 = { 0, 6, 0, fmt_1650, 0 };




    for (i__ = 1; i__ <= 15; ++i__) {
	if (list[i__ - 1] == *number) {
	    goto L150;
	}
/* L50: */
    }
    s_wsfe(&io___291);
    e_wsfe();
    s_stop("", (ftnlen)0);
L150:

    switch (i__) {
	case 1:  goto L200;
	case 2:  goto L300;
	case 3:  goto L400;
	case 4:  goto L500;
	case 5:  goto L600;
	case 6:  goto L700;
	case 7:  goto L800;
	case 8:  goto L900;
	case 9:  goto L1000;
	case 10:  goto L1100;
	case 11:  goto L1200;
	case 12:  goto L1300;
	case 13:  goto L1400;
	case 14:  goto L1500;
	case 15:  goto L1600;
    }

L200:
    s_wsfe(&io___292);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*float__), (ftnlen)sizeof(real));
    e_wsfe();
    s_stop("", (ftnlen)0);
L300:
    s_wsfe(&io___293);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*int__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
L400:
    s_wsfe(&io___294);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*int__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
L500:
    s_wsfe(&io___295);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*float__), (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&(*int__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
L600:
    s_wsfe(&io___296);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*int__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
L700:
    s_wsfe(&io___297);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*int__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
L800:
    s_wsfe(&io___298);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*int__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
L900:
    s_wsfe(&io___299);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*int__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
L1000:
    s_wsfe(&io___300);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*int__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
L1100:
    s_wsfe(&io___301);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*int__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
L1200:
    s_wsfe(&io___302);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*int__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
L1300:
    s_wsfe(&io___303);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*int__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
L1400:
    s_wsfe(&io___304);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*int__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
L1500:
    s_wsfe(&io___305);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*int__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
L1600:
    s_wsfe(&io___306);
    do_fio(&c__1, (char *)&(*number), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*int__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
    return 0;
} /* error_ */

