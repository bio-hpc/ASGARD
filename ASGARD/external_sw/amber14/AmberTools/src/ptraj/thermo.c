/* thermo.f -- translated by f2c (version 20030306).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "../f2c/f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__9 = 9;
static integer c__3 = 3;
static integer c__5 = 5;
static integer c__2 = 2;


/* ************************************************************************ */
/*                              AMBER                                   ** */
/*                                                                      ** */
/*                  Copyright (c) 1986, 1991, 1995                      ** */
/*             Regents of the University of California                  ** */
/*                       All Rights Reserved.                           ** */
/*                                                                      ** */
/*  This software provided pursuant to a license agreement containing   ** */
/*  restrictions on its disclosure, duplication, and use. This software ** */
/*  contains confidential and proprietary information, and may not be   ** */
/*  extracted or distributed, in whole or in part, for any purpose      ** */
/*  whatsoever, without the express written permission of the authors.  ** */
/*  This notice, and the associated author list, must be attached to    ** */
/*  all copies, or extracts, of this software. Any additional           ** */
/*  restrictions set forth in the license agreement also apply to this  ** */
/*  software.                                                           ** */
/* ************************************************************************ */

/* Subroutine */ int thermo_(integer *natoms, integer *nvecs, integer *ilevel,
	 doublereal *c__, doublereal *amass, doublereal *freq, doublereal *
	vtemp, doublereal *evibn, doublereal *cvibn, doublereal *svibn, 
	doublereal *t, doublereal *patm)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal akilo = 1e3;
    static doublereal thresh = 900.;
    static integer iout = 6;
    static doublereal pstd = 101325.;
    static doublereal pt2 = .2;
    static doublereal half = .5;
    static doublereal one = 1.;
    static doublereal onept5 = 1.5;
    static doublereal two = 2.;
    static doublereal twopt5 = 2.5;
    static doublereal four = 4.;
    static doublereal eight = 8.;

    /* Format strings */
    static char fmt_1000[] = "(/20x,19(\002*\002),/20x,\002- Thermochemistry"
	    " -\002,/20x,19(\002*\002),//)";
    static char fmt_1010[] = "(1x,\002molecular mass (principal isotopes)"
	    " \002,f11.5,\002 amu\002)";
    static char fmt_1020[] = "(/1x,\002temperature \002,f9.3,\002 kelvin\002"
	    ",/1x,\002pressure    \002,f9.5,\002 atm\002)";
    static char fmt_1040[] = "(/1x,\002internal energy:   \002,f10.3,\002 jo"
	    "ule/mol\002,9x,f10.3,\002 kcal/mol\002/1x,\002entropy:           "
	    "\002,f10.3,\002 joule/k-mol\002,7x,f10.3,\002 cal/k-mol\002/1x"
	    ",\002heat capacity cv:  \002,f10.3,\002 joule/k-mol\002,7x,f10.3,"
	    "\002 cal/k-mol\002)";
    static char fmt_1050[] = "(/1x,\002rotational symmetry number \002,f3.0)";
    static char fmt_1060[] = "(/1x,\002Warning-- assumption of classical beh"
	    "avior for \002,\002rotation\002,/1x,\002          may cause sign"
	    "ificant error\002/)";
    static char fmt_1070[] = "(/1x,\002rotational temperatures (kelvin) \002"
	    ",3f12.5)";
    static char fmt_1080[] = "(/1x,\002rotational temperature (kelvin) \002,"
	    "f12.5)";
    static char fmt_1090[] = "(/1x,\002zero point vibrational energy \002,f1"
	    "2.1,\002 (joules/mol) \002,/1x,30x,f12.5,\002 (kcal/mol)\002,/1x"
	    ",30x,f12.7,\002 (hartree/particle)\002)";
    static char fmt_1100[] = "(/1x,\002Warning-- \002,i3,\002 vibrations hav"
	    "e low frequencies\002,\002 and may represent hindered \002,/1x"
	    ",\002        internal rotations.  The contributions \002,\002pri"
	    "nted below assume that these \002,/1x,\002        really are vib"
	    "rations.\002)";
    static char fmt_1130[] = "(//1x,10x,\002freq.\002,9x,\002E\002,9x,9x,"
	    "\002Cv\002,8x,9x,\002S\002)";
    static char fmt_1150[] = "(1x,\002Total\002,10x,3(4x,f11.3,4x))";
    static char fmt_1160[] = "(1x,\002translational\002,2x,3(4x,f11.3,4x))";
    static char fmt_1170[] = "(1x,\002rotational\002,5x,3(4x,f11.3,4x))";
    static char fmt_1180[] = "(1x,\002vibrational\002,4x,3(4x,f11.3,4x))";
    static char fmt_1190[] = "(1x,i5,f10.3,3(4x,f11.3,4x))";
    static char fmt_1200[] = "(1x,9x,\002cm**-1\002,6x,\002kcal/mol\002,5x,3"
	    "x,\002cal/mol-kelvin\002,2x,2x,\002cal/mol-kelvin\002,/8(\002---"
	    "-------\002))";
    static char fmt_1220[] = "(/1x,\002principal moments of inertia (nuclei "
	    "only) in \002,\002amu-A**2:\002,/1x,5x,3f12.2)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), exp(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), log(
	    doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static doublereal e;
    static integer i__;
    static doublereal p, s, pi, cv, sn, rt, em1;
    static integer iff;
    static doublereal arg, gas;
    static integer iat;
    static doublereal con, dum, ezj, dum1, dum2, argd, cvib, evib;
    static integer ndof;
    extern /* Subroutine */ int mofi_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal avog, ezkc, pipi, ezpe, tokg, ezau, svib, crot, pmom[10]
	    , erot, etot, ctot, srot, stot, tovt, jpcal, tocal, ccont, ctran, 
	    econt, etran, scont, stran, tomet, rtemp, boltz, etovt, rtemp1, 
	    rtemp2, rtemp3, planck;
    static logical linear;
    static doublereal tokcal, hartre, weight;
    static integer lofreq;
    extern /* Subroutine */ int symnum_(integer *, doublereal *, doublereal *,
	     logical *);

    /* Fortran I/O blocks */
    static cilist io___28 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_1020, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_1040, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_1220, 0 };
    static cilist io___48 = { 0, 0, 0, fmt_1050, 0 };
    static cilist io___51 = { 0, 0, 0, fmt_1060, 0 };
    static cilist io___52 = { 0, 0, 0, fmt_1080, 0 };
    static cilist io___56 = { 0, 0, 0, fmt_1060, 0 };
    static cilist io___57 = { 0, 0, 0, fmt_1070, 0 };
    static cilist io___68 = { 0, 0, 0, fmt_1090, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_1100, 0 };
    static cilist io___81 = { 0, 6, 0, 0, 0 };
    static cilist io___85 = { 0, 0, 0, fmt_1130, 0 };
    static cilist io___86 = { 0, 0, 0, fmt_1200, 0 };
    static cilist io___87 = { 0, 0, 0, fmt_1150, 0 };
    static cilist io___88 = { 0, 0, 0, fmt_1160, 0 };
    static cilist io___89 = { 0, 0, 0, fmt_1170, 0 };
    static cilist io___90 = { 0, 0, 0, fmt_1180, 0 };
    static cilist io___91 = { 0, 0, 0, fmt_1190, 0 };
    static cilist io___92 = { 0, 0, 0, fmt_1190, 0 };



/*     given the structure of a molecule and its normal mode vibrational */
/*     frequencies this routine uses standard statistical mechanical */
/*     formulas for an ideal gas (in the canonical ensemble, see, */
/*     for example, d. a. mcquarrie, "statistical thermodynamics", */
/*     harper & row, new york, 1973, chapters 5, 6, and 8) to compute */
/*     the entropy, heat capacity, and internal energy. */

/*     the si system of units is used internally.  conversion to units */
/*     more familiar to most chemists is made for output. */


/*     amass:   atomic weights, in amu. */
/*     pmom:    principal moments of inertia, in amu-bohr**2 and */
/*              in ascending order. */
/*     freq:    vibrational frequencies, in cm**-1 and in ascending */
/*              order */
/*     c    :   coordinates in Angstroms */
/*     vtemp:   vibrational temperatures, in kelvin. */
/*     evibn:   contribution to e from the vibration n. */
/*     cvibn:   contribution to cv from the vibration n. */
/*     svibn:   contribution to s from the vibration n. */
/*     t:       temperature */
/*     patm:    pressure, in atmospheres */


    /* Parameter adjustments */
    --svibn;
    --cvibn;
    --evibn;
    --vtemp;
    --freq;
    --amass;
    --c__;

    /* Function Body */

/* L1000: */
/* L1010: */
/* L1020: */
/* 1030 format(/1x,'Warning-- assumptions made about the electronic ', */
/*    +           'partition function', */
/*    +       /1x,'          are not valid for multiplets!'/) */
/* L1040: */
/* L1050: */
/* L1060: */
/* L1070: */
/* L1080: */
/* L1090: */
/* L1100: */
/* L1110: */
/* L1120: */
/* L1125: */
/* L1130: */
/* L1140: */
/* L1150: */
/* L1160: */
/* L1170: */
/* L1180: */
/* L1190: */
/* L1200: */
/* L1210: */
/* L1220: */

/*     tokg:    kilograms per amu. */
/*     boltz:   boltzman constant, in joules per kelvin. */
/*     planck:  planck constant, in joule-seconds. */
/*     avog:    avogadro constant, in mol**(-1). */
/*     jpcal:   joules per calorie. */
/*     tomet:   metres per Angstrom. */
/*     hartre:  joules per hartree. */

    tokg = 1.660531e-27;
    boltz = 1.380622e-23;
    planck = 6.626196e-34;
    avog = 6.022169e23;
    jpcal = 4.18674;
    tomet = 1e-10;
    hartre = 4.35981e-18;

/*     compute the gas constant, pi, pi**2, and e. */
/*     compute the conversion factors cal per joule and kcal per joule. */

    gas = avog * boltz;
    pi = four * atan(one);
    pipi = pi * pi;
    e = exp(one);
    tocal = one / jpcal;
    tokcal = tocal / akilo;

/*     print the temperature and pressure. */

    p = pstd * *patm;
    io___28.ciunit = iout;
    s_wsfe(&io___28);
    e_wsfe();
    io___29.ciunit = iout;
    s_wsfe(&io___29);
    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*patm), (ftnlen)sizeof(doublereal));
    e_wsfe();
    rt = gas * *t;

/*     compute and print the molecular mass in amu, then convert to */
/*     kilograms. */

    weight = zero;
    i__1 = *natoms;
    for (iat = 1; iat <= i__1; ++iat) {
	weight += amass[iat];
/* L20: */
    }
    io___33.ciunit = iout;
    s_wsfe(&io___33);
    do_fio(&c__1, (char *)&weight, (ftnlen)sizeof(doublereal));
    e_wsfe();
    weight *= tokg;

/*     trap non-unit multiplicities. */

/*     if (multip .ne. 1) write(iout,1030) */

/*     compute contributions due to translation: */
/*        etran-- internal energy */
/*        ctran-- constant v heat capacity */
/*        stran-- entropy */

    dum1 = boltz * *t;
    d__1 = two * pi;
    dum2 = pow_dd(&d__1, &onept5);
    arg = pow_dd(&dum1, &onept5) / planck;
    arg = arg / p * (dum1 / planck);
    arg = arg * dum2 * (weight / planck);
    arg = arg * sqrt(weight) * pow_dd(&e, &twopt5);
    stran = gas * log(arg);
    etran = onept5 * rt;
    ctran = onept5 * gas;

/*     Compute contributions due to electronic motion: */
/*        It is assumed that the first electronic excitation energy */
/*        is much greater than kt and that the ground state has a */
/*        degeneracy of one.  Under these conditions the electronic */
/*        partition function can be considered to be unity.  The */
/*        ground electronic state is taken to be the zero of */
/*        electronic energy. */

/*  40 continue */

/*     for monatomics print and return. */

    if (*natoms <= 1) {
	s = stran * tocal;
	e = etran * tokcal;
	cv = ctran * tocal;
	io___42.ciunit = iout;
	s_wsfe(&io___42);
	do_fio(&c__1, (char *)&etran, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&stran, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ctran, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&cv, (ftnlen)sizeof(doublereal));
	e_wsfe();
	return 0;
    }

/*     compute contributions due to rotation. */

/*     Compute the principal moments of inertia, get the rotational */
/*     symmetry number, see if the molecule is linear, and compute */
/*     the rotational temperatures.  Note the imbedded conversion */
/*     of the moments to SI units. */

    mofi_(natoms, &c__[1], &amass[1], pmom);
    io___44.ciunit = iout;
    s_wsfe(&io___44);
    for (i__ = 1; i__ <= 3; ++i__) {
	do_fio(&c__1, (char *)&pmom[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    linear = FALSE_;
    symnum_(natoms, &amass[1], &sn, &linear);
    io___48.ciunit = iout;
    s_wsfe(&io___48);
    do_fio(&c__1, (char *)&sn, (ftnlen)sizeof(doublereal));
    e_wsfe();
    con = planck / (boltz * eight * pipi);
    con = con / tokg * (planck / (tomet * tomet));
    if (linear) {
	rtemp = con / pmom[2];
	if (rtemp < pt2) {
	    io___51.ciunit = iout;
	    s_wsfe(&io___51);
	    e_wsfe();
	}
	io___52.ciunit = iout;
	s_wsfe(&io___52);
	do_fio(&c__1, (char *)&rtemp, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else {
	rtemp1 = con / pmom[0];
	rtemp2 = con / pmom[1];
	rtemp3 = con / pmom[2];
	if (rtemp1 < pt2) {
	    io___56.ciunit = iout;
	    s_wsfe(&io___56);
	    e_wsfe();
	}
	io___57.ciunit = iout;
	s_wsfe(&io___57);
	do_fio(&c__1, (char *)&rtemp1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rtemp2, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rtemp3, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*         erot-- rotational contribution to internal energy. */
/*         crot-- rotational contribution to cv. */
/*         srot-- rotational contribution to entropy. */

    if (linear) {
	erot = rt;
	crot = gas;
	arg = *t / rtemp * (e / sn);
	srot = gas * log(arg);
    } else {
	erot = onept5 * rt;
	crot = onept5 * gas;
	arg = sqrt(pi * e * e * e) / sn;
	dum = *t / rtemp1 * (*t / rtemp2) * (*t / rtemp3);
	arg *= sqrt(dum);
	srot = gas * log(arg);
    }

/*     compute contributions due to vibration. */

/*     compute vibrational temperatures and zero point vibrational */
/*     energy.  only real frequencies are included in the analysis. */

/*     ndof = 3*natoms - 6 - nimag */
/*     if (nimag .ne. 0) write(iout,1210) nimag */
/*     if (linear) ndof = ndof + 1 */
    ndof = *nvecs;

/*       (---iff is the first frequency to include in thermo:) */

    if (*ilevel != 0) {
	iff = 0;
    } else if (linear) {
	iff = 5;
    } else {
	iff = 6;
    }
    con = planck / boltz;
    ezpe = zero;
    i__1 = ndof;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vtemp[i__] = freq[i__ + iff] * con * 3e10;
	ezpe += freq[i__ + iff] * 3e10;
/* L160: */
    }
    ezpe = half * planck * ezpe;
    ezj = ezpe * avog;
    ezkc = ezpe * tokcal * avog;
    ezau = ezpe / hartre;
    io___68.ciunit = iout;
    s_wsfe(&io___68);
    do_fio(&c__1, (char *)&ezj, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ezkc, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ezau, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     compute the number of vibrations for which more than 5% of an */
/*     assembly of molecules would exist in vibrational excited states. */
/*     special printing for these modes is done to allow the user to */
/*     easily take internal rotations into account.  the criterion */
/*     corresponds roughly to a low frequency of 1.9(10**13) hz, or */
/*     625 cm**(-1), or a vibrational temperature of 900 k. */

    lofreq = 0;
    i__1 = ndof;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (vtemp[i__] < thresh) {
	    ++lofreq;
	}
/* L180: */
    }
    if (lofreq != 0) {
	io___70.ciunit = iout;
	s_wsfe(&io___70);
	do_fio(&c__1, (char *)&lofreq, (ftnlen)sizeof(integer));
	e_wsfe();
    }


/*     compute: */
/*        evib-- the vibrational component of the internal energy. */
/*        cvib-- the vibrational component of the heat capacity. */
/*        svib-- the vibrational component of the entropy. */

    evib = zero;
    cvib = zero;
    svib = zero;
    i__1 = ndof;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*       compute some common factors. */

	tovt = vtemp[i__] / *t;
	etovt = exp(tovt);
	em1 = etovt - one;

/*       compute contributions due to the i'th vibration. */

	econt = tovt * (half + one / em1);
/* Computing 2nd power */
	d__1 = tovt / em1;
	ccont = etovt * (d__1 * d__1);
	argd = one - one / etovt;
	if (argd > 1e-7) {
	    scont = tovt / em1 - log(argd);
	} else {
	    scont = 0.f;
	    s_wsle(&io___81);
	    do_lio(&c__9, &c__1, "warning: setting vibrational entropy to ze"
		    "ro ", (ftnlen)45);
	    do_lio(&c__9, &c__1, "for mode ", (ftnlen)9);
	    do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " with vtemp = ", (ftnlen)14);
	    do_lio(&c__5, &c__1, (char *)&vtemp[i__], (ftnlen)sizeof(
		    doublereal));
	    e_wsle();
	}
/*       if (lofreq .ge. i) then */
	evibn[i__] = econt * rt;
	cvibn[i__] = ccont * gas;
	svibn[i__] = scont * gas;
/*       end if */
	evib += econt;
	cvib += ccont;
	svib += scont;
/* L220: */
    }
    evib *= rt;
    cvib *= gas;
    svib *= gas;

/*     the units are now: */
/*         e-- joules/mol */
/*         c-- joules/mol-kelvin */
/*         s-- joules/mol-kelvin */

    etot = etran + erot + evib;
    ctot = ctran + crot + cvib;
    stot = stran + srot + svib;

/*     print the sum of the hartree-fock energy and the thermal energy. */

/*     call tread(501,gen,47,1,47,1,0) */
/*     esum = gen(32) + etot/avog/hartre */
/*     write(iout,1230) esum */


/*     convert to the following and print */
/*         e-- kcal/mol */
/*         c-- cal/mol-kelvin */
/*         s-- cal/mol-kelvin */

/* L240: */
    etran *= tokcal;
    ctran *= tocal;
    stran *= tocal;
    erot *= tokcal;
    crot *= tocal;
    srot *= tocal;
    evib *= tokcal;
    cvib *= tocal;
    svib *= tocal;
    etot = etran + erot + evib;
    ctot = ctran + crot + cvib;
    stot = stran + srot + svib;
    i__1 = ndof;
    for (i__ = 1; i__ <= i__1; ++i__) {
	evibn[i__] *= tokcal;
	cvibn[i__] *= tocal;
	svibn[i__] *= tocal;
/* L280: */
    }

    io___85.ciunit = iout;
    s_wsfe(&io___85);
    e_wsfe();
    io___86.ciunit = iout;
    s_wsfe(&io___86);
    e_wsfe();
    io___87.ciunit = iout;
    s_wsfe(&io___87);
    do_fio(&c__1, (char *)&etot, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ctot, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&stot, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___88.ciunit = iout;
    s_wsfe(&io___88);
    do_fio(&c__1, (char *)&etran, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ctran, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&stran, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___89.ciunit = iout;
    s_wsfe(&io___89);
    do_fio(&c__1, (char *)&erot, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&crot, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&srot, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___90.ciunit = iout;
    s_wsfe(&io___90);
    do_fio(&c__1, (char *)&evib, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&cvib, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&svib, (ftnlen)sizeof(doublereal));
    e_wsfe();
    i__1 = iff;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___91.ciunit = iout;
	s_wsfe(&io___91);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&freq[i__], (ftnlen)sizeof(doublereal));
	e_wsfe();
/* L318: */
    }
    i__1 = ndof;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___92.ciunit = iout;
	s_wsfe(&io___92);
	i__2 = i__ + iff;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&freq[i__ + iff], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&evibn[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&cvibn[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&svibn[i__], (ftnlen)sizeof(doublereal));
	e_wsfe();
/* L320: */
    }

    return 0;
} /* thermo_ */

/* ----------------------------------------------------------------------- */
/* Subroutine */ int symnum_(integer *natoms, doublereal *amass, doublereal *
	sn, logical *linear)
{

/*     ----- routine to give the symmetry number. only for linear */
/*           molecules. for others symmetry number is unity ----- */


    /* Parameter adjustments */
    --amass;

    /* Function Body */
    *sn = 1.;
    if (*natoms <= 2) {
	*linear = TRUE_;
	if (amass[1] == amass[2]) {
	    *sn = 2.;
	}
    }
    return 0;
} /* symnum_ */

/* ----------------------------------------------------------------------- */
/* Subroutine */ int mofi_(integer *natoms, doublereal *c__, doublereal *
	amass, doublereal *pmom)
{
    /* Initialized data */

    static char uplo[1] = "U";
    static char jobz[1] = "V";
    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal t[9], x, y, z__, e2[30], wt;
    static integer iat;
    static doublereal com[3];
    static integer ier, iaind;
    extern /* Subroutine */ int dspev_(char *, char *, integer *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static doublereal totwt, eigvec[9];


/*     compute the principal moments of inertia. */
/*     units are amu-bohr**2 */

    /* Parameter adjustments */
    --pmom;
    --amass;
    --c__;

    /* Function Body */


/*     stetement function to get ccom. */



/*     compute the position of the center of mass and translate */
/*     it to the origin. */

    com[0] = zero;
    com[1] = zero;
    com[2] = zero;

    totwt = zero;
    i__1 = *natoms;
    for (iat = 1; iat <= i__1; ++iat) {
	iaind = (iat - 1) * 3;
	wt = amass[iat];
	totwt += wt;
	com[0] += wt * c__[iaind + 1];
	com[1] += wt * c__[iaind + 2];
	com[2] += wt * c__[iaind + 3];
/* L20: */
    }

    com[0] /= totwt;
    com[1] /= totwt;
    com[2] /= totwt;

/*     compute the principal moments. */

    for (i__ = 1; i__ <= 9; ++i__) {
	t[i__ - 1] = zero;
/* L60: */
    }

    i__1 = *natoms;
    for (iat = 1; iat <= i__1; ++iat) {
	wt = amass[iat];
	x = c__[c__1 + 3 * (iat - 1)] - com[c__1 - 1];
	y = c__[c__2 + 3 * (iat - 1)] - com[c__2 - 1];
	z__ = c__[c__3 + 3 * (iat - 1)] - com[c__3 - 1];
	t[0] += wt * (y * y + z__ * z__);
	t[2] += wt * (x * x + z__ * z__);
	t[5] += wt * (x * x + y * y);
	t[1] -= wt * x * y;
	t[3] -= wt * x * z__;
	t[4] -= wt * y * z__;
/* L80: */
    }
    ier = 0;
    dspev_(jobz, uplo, &c__3, t, &pmom[1], eigvec, &c__3, e2, &ier, (ftnlen)1,
	     (ftnlen)1);
    return 0;
} /* mofi_ */

