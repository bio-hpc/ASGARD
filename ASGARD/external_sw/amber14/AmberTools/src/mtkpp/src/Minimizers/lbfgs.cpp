/*! 
   \file lbfgs.cpp
   \brief This class contains code for the limited-memory Broyden-Fletcher-Goldfarb-Shanno (LBFGS) algorithm
    for large-scale multidimensional unconstrained minimization problems. This file is a translation of 
    Fortran code written by Jorge Nocedal. It is distributed as part of the RISO project.

    This code is derived from the Fortran program lbfgs.f.
    The C++ translation was effected mostly mechanically, with some manual clean-up.

    Here's some information on the original LBFGS Fortran source code, 
    available at http://www.netlib.org/opt/lbfgs_um.shar.
    This info is taken verbatim from the Netlib blurb on the Fortran source.

    LBFGS minimizer from the original lbfgs.f, converted to c++ by Kenneth Ayers

   $Date: 2008/01/10 20:56:34 $
   $Revision: 1.4 $

   ----------------------------------------------------------------------------
*/

#include "lbfgs.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

/* Table of constant values */
namespace MTKpp
{

static int c__1 = 1;

/*     ---------------------------------------------------------------------- */
/*     This file contains the LBFGS algorithm and supporting routines */

/*     **************** */
/*     LBFGS SUBROUTINE */
/*     **************** */

/* Subroutine */ 
void lbfgs_(int *n, int *m, double *x, double *f, double *g, 
    int *diagco, double *diag, int *iprint, double *eps, 
    double *xtol, double *w, int *iflag) {
    /* Initialized data */

    static double one = 1.0;
    static double zero = 0.0;

    /* System generated locals */
    int i__1;
    double d__1;

    /* Local variables */
    static double beta;
    static int inmc;
    static int info, iscn, nfev, iycn, iter;
    static double ftol;
    static int nfun, ispt, iypt, i__, bound;
    static double gnorm;
    static int point;
    static double xnorm;
    static int cp;
    static double sq, yr, ys;
    static int finish;
    static double yy;
    static int maxfev;
    static int npt;
    static double stp, stp1;

/*        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION */
/*                          JORGE NOCEDAL */
/*                        *** July 1990 *** */


/*     This subroutine solves the unconstrained minimization problem */

/*                      min F(x),    x= (x1,x2,...,xN), */

/*      using the limited memory BFGS method. The routine is especially */
/*      effective on problems involving a large number of variables. In */
/*      a typical iteration of this method an approximation Hk to the */
/*      inverse of the Hessian is obtained by applying M BFGS updates to */
/*      a diagonal matrix Hk0, using information from the previous M steps. */
/*      The user specifies the number M, which determines the amount of */
/*      storage required by the routine. The user may also provide the */
/*      diagonal matrices Hk0 if not satisfied with the default choice. */
/*      The algorithm is described in "On the limited memory BFGS method */
/*      for large scale optimization", by D. Liu and J. Nocedal, */
/*      Mathematical Programming B 45 (1989) 503-528. */

/*      The user is required to calculate the function value F and its */
/*      gradient G. In order to allow the user complete control over */
/*      these computations, reverse  communication is used. The routine */
/*      must be called repeatedly under the control of the parameter */
/*      IFLAG. */

/*      The steplength is determined at each iteration by means of the */
/*      line search routine MCVSRCH, which is a slight modification of */
/*      the routine CSRCH written by More' and Thuente. */

/*      The calling statement is */

/*          CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG) */

/*      where */

/*     N       is an INTEGER variable that must be set by the user to the */
/*             number of variables. It is not altered by the routine. */
/*             Restriction: N>0. */

/*     M       is an INTEGER variable that must be set by the user to */
/*             the number of corrections used in the BFGS update. It */
/*             is not altered by the routine. Values of M less than 3 are */
/*             not recommended; large values of M will result in excessive */
/*             computing time. 3<= M <=7 is recommended. Restriction: M>0. */

/*     X       is a DOUBLE PRECISION array of length N. On initial entry */
/*             it must be set by the user to the values of the initial */
/*             estimate of the solution vector. On exit with IFLAG=0, it */
/*             contains the values of the variables at the best point */
/*             found (usually a solution). */

/*     F       is a DOUBLE PRECISION variable. Before initial entry and on */
/*             a re-entry with IFLAG=1, it must be set by the user to */
/*             contain the value of the function F at the point X. */

/*     G       is a DOUBLE PRECISION array of length N. Before initial */
/*             entry and on a re-entry with IFLAG=1, it must be set by */
/*             the user to contain the components of the gradient G at */
/*             the point X. */

/*     DIAGCO  is a LOGICAL variable that must be set to .TRUE. if the */
/*             user  wishes to provide the diagonal matrix Hk0 at each */
/*             iteration. Otherwise it should be set to .FALSE., in which */
/*             case  LBFGS will use a default value described below. If */
/*             DIAGCO is set to .TRUE. the routine will return at each */
/*             iteration of the algorithm with IFLAG=2, and the diagonal */
/*              matrix Hk0  must be provided in the array DIAG. */


/*     DIAG    is a DOUBLE PRECISION array of length N. If DIAGCO=.TRUE., */
/*             then on initial entry or on re-entry with IFLAG=2, DIAG */
/*             it must be set by the user to contain the values of the */
/*             diagonal matrix Hk0.  Restriction: all elements of DIAG */
/*             must be positive. */

/*     IPRINT  is an INTEGER array of length two which must be set by the */
/*             user. */

/*             IPRINT(1) specifies the frequency of the output: */
/*                IPRINT(1) < 0 : no output is generated, */
/*                IPRINT(1) = 0 : output only at first and last iteration, */
/*                IPRINT(1) > 0 : output every IPRINT(1) iterations. */

/*             IPRINT(2) specifies the type of output generated: */
/*                IPRINT(2) = 0 : iteration count, number of function */
/*                                evaluations, function value, norm of the */
/*                                gradient, and steplength, */
/*                IPRINT(2) = 1 : same as IPRINT(2)=0, plus vector of */
/*                                variables and  gradient vector at the */
/*                                initial point, */
/*                IPRINT(2) = 2 : same as IPRINT(2)=1, plus vector of */
/*                                variables, */
/*                IPRINT(2) = 3 : same as IPRINT(2)=2, plus gradient vector. */


/*     EPS     is a positive DOUBLE PRECISION variable that must be set by */
/*             the user, and determines the accuracy with which the solution */
/*             is to be found. The subroutine terminates when */

/*                         ||G|| < EPS max(1,||X||), */

/*             where ||.|| denotes the Euclidean norm. */

/*     XTOL    is a  positive DOUBLE PRECISION variable that must be set by */
/*             the user to an estimate of the machine precision (e.g. */
/*             10**(-16) on a SUN station 3/60). The line search routine will */
/*             terminate if the relative width of the interval of uncertainty */
/*             is less than XTOL. */

/*     W       is a DOUBLE PRECISION array of length N(2M+1)+2M used as */
/*             workspace for LBFGS. This array must not be altered by the */
/*             user. */

/*     IFLAG   is an INTEGER variable that must be set to 0 on initial entry */
/*             to the subroutine. A return with IFLAG<0 indicates an error, */
/*             and IFLAG=0 indicates that the routine has terminated without */
/*             detecting errors. On a return with IFLAG=1, the user must */
/*             evaluate the function F and gradient G. On a return with */
/*             IFLAG=2, the user must provide the diagonal matrix Hk0. */

/*             The following negative values of IFLAG, detecting an error, */
/*             are possible: */

/*              IFLAG=-1  The line search routine MCSRCH failed. The */
/*                        parameter INFO provides more detailed information */
/*                        (see also the documentation of MCSRCH): */

/*                       INFO = 0  IMPROPER INPUT PARAMETERS. */

/*                       INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF */
/*                                 UNCERTAINTY IS AT MOST XTOL. */

/*                       INFO = 3  MORE THAN 20 FUNCTION EVALUATIONS WERE */
/*                                 REQUIRED AT THE PRESENT ITERATION. */

/*                       INFO = 4  THE STEP IS TOO SMALL. */

/*                       INFO = 5  THE STEP IS TOO LARGE. */

/*                       INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS. */
/*                                 THERE MAY NOT BE A STEP WHICH SATISFIES */
/*                                 THE SUFFICIENT DECREASE AND CURVATURE */
/*                                 CONDITIONS. TOLERANCES MAY BE TOO SMALL. */


/*              IFLAG=-2  The i-th diagonal element of the diagonal inverse */
/*                        Hessian approximation, given in DIAG, is not */
/*                        positive. */

/*              IFLAG=-3  Improper input parameters for LBFGS (N or M are */
/*                        not positive). */



/*    ON THE DRIVER: */

/*    The program that calls LBFGS must contain the declaration: */

/*                       EXTERNAL LB2 */

/*    LB2 is a BLOCK DATA that defines the default values of several */
/*    parameters described in the COMMON section. */



/*    COMMON: */

/*     The subroutine contains one common area, which the user may wish to */
/*    reference: */


/*    MP  is an INTEGER variable with default value 6. It is used as the */
/*        unit number for the printing of the monitoring information */
/*        controlled by IPRINT. */

/*    LP  is an INTEGER variable with default value 6. It is used as the */
/*        unit number for the printing of error messages. This printing */
/*        may be suppressed by setting LP to a non-positive value. */

/*    GTOL is a DOUBLE PRECISION variable with default value 0.9, which */
/*        controls the accuracy of the line search routine MCSRCH. If the */
/*        function and gradient evaluations are inexpensive with respect */
/*        to the cost of the iteration (which is sometimes the case when */
/*        solving very large problems) it may be advantageous to set GTOL */
/*        to a small value. A typical small value is 0.1.  Restriction: */
/*        GTOL should be greater than 1.D-04. */

/*    STPMIN and STPMAX are non-negative DOUBLE PRECISION variables which */
/*        specify lower and uper bounds for the step in the line search. */
/*        Their default values are 1.D-20 and 1.D+20, respectively. These */
/*        values need not be modified unless the exponents are too large */
/*        for the machine being used, or unless the problem is extremely */
/*        badly scaled (in which case the exponents should be increased). */


/*  MACHINE DEPENDENCIES */

/*        The only variables that are machine-dependent are XTOL, */
/*        STPMIN and STPMAX. */


/*  GENERAL INFORMATION */

/*    Other routines called directly:  DAXPY, DDOT, LB1, MCSRCH */

/*    Input/Output  :  No input; diagnostic messages on unit MP and */
/*                     error messages on unit LP. */


/*     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */


    /* Parameter adjustments */
    --diag;
    --g;
    --x;
    --w;
    --iprint;

    /* Function Body */

/*     INITIALIZE */
/*     ---------- */

    if (*iflag == 0) {
      goto L10;
    }
    switch (*iflag) {
      case 1:  goto L172;
      case 2:  goto L100;
    }
L10:
    iter = 0;
    if (*n <= 0 || *m <= 0) {
      goto L196;
    }
/*
    if (lb3_1.gtol <= 1e-4) {
      if (lb3_1.lp > 0) {
        std::cout << "  GTOL IS LESS THAN OR EQUAL TO 1.D-04" << std::endl
                  << " IT HAS BEEN RESET TO 9.D-01" << std::endl;
        //io___4.ciunit = lb3_1.lp;
        //s_wsfe(&io___4);
        //e_wsfe();
      }
      lb3_1.gtol = .9;
    }
*/
    nfun = 1;
    point = 0;
    finish = 0;
    if (*diagco) {
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
        if (diag[i__] <= zero) {
          goto L195;
        }
      }
    } else {
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
      diag[i__] = 1.;
      }
    }

/*     THE WORK VECTOR W IS DIVIDED AS FOLLOWS: */
/*     --------------------------------------- */
/*     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND */
/*         OTHER TEMPORARY INFORMATION. */
/*     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO. */
/*     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED */
/*         IN THE FORMULA THAT COMPUTES H*G. */
/*     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH */
/*         STEPS. */
/*     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M */
/*         GRADIENT DIFFERENCES. */

/*     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A */
/*     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT. */

    ispt = *n + (*m << 1);
    iypt = ispt + *n * *m;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L50: */
      w[ispt + i__] = -g[i__] * diag[i__];
    }
    gnorm = sqrt(ddot_(n, &g[1], &c__1, &g[1], &c__1));
    stp1 = one / gnorm;

/*     PARAMETERS FOR LINE SEARCH ROUTINE */

    ftol = 1e-4;
    maxfev = 20;

    if (iprint[1] >= 0) {
	lb1_(&iprint[1], &iter, &nfun, &gnorm, n, m, &x[1], f, &g[1], &stp, &
		finish);
    }

/*    -------------------- */
/*     MAIN ITERATION LOOP */
/*    -------------------- */

L80:
    ++iter;
    info = 0;
    bound = iter - 1;
    if (iter == 1) {
	goto L165;
    }
    if (iter > *m) {
	bound = *m;
    }

    ys = ddot_(n, &w[iypt + npt + 1], &c__1, &w[ispt + npt + 1], &c__1);
    if (! (*diagco)) {
	yy = ddot_(n, &w[iypt + npt + 1], &c__1, &w[iypt + npt + 1], &c__1);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    diag[i__] = ys / yy;
	}
    } else {
	*iflag = 2;
	return;
    }
L100:
    if (*diagco) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	    if (diag[i__] <= zero) {
		goto L195;
	    }
	}
    }

/*     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980, */
/*     "Updating quasi-Newton matrices with limited storage", */
/*     Mathematics of Computation, Vol.24, No.151, pp. 773-782. */
/*     --------------------------------------------------------- */

    cp = point;
    if (point == 0) {
	cp = *m;
    }
    w[*n + cp] = one / ys;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L112: */
	w[i__] = -g[i__];
    }
    cp = point;
    i__1 = bound;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--cp;
	if (cp == -1) {
	    cp = *m - 1;
	}
	sq = ddot_(n, &w[ispt + cp * *n + 1], &c__1, &w[1], &c__1);
	inmc = *n + *m + cp + 1;
	iycn = iypt + cp * *n;
	w[inmc] = w[*n + cp + 1] * sq;
	d__1 = -w[inmc];
	daxpy_(n, &d__1, &w[iycn + 1], &c__1, &w[1], &c__1);
/* L125: */
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L130: */
	w[i__] = diag[i__] * w[i__];
    }

    i__1 = bound;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yr = ddot_(n, &w[iypt + cp * *n + 1], &c__1, &w[1], &c__1);
	beta = w[*n + cp + 1] * yr;
	inmc = *n + *m + cp + 1;
	beta = w[inmc] - beta;
	iscn = ispt + cp * *n;
	daxpy_(n, &beta, &w[iscn + 1], &c__1, &w[1], &c__1);
	++cp;
	if (cp == *m) {
	    cp = 0;
	}
/* L145: */
    }

/*     STORE THE NEW SEARCH DIRECTION */
/*     ------------------------------ */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L160: */
	w[ispt + point * *n + i__] = w[i__];
    }

/*     OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION */
/*     BY USING THE LINE SEARCH ROUTINE MCSRCH */
/*     ---------------------------------------------------- */
L165:
    nfev = 0;
    stp = one;
    if (iter == 1) {
	stp = stp1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L170: */
	w[i__] = g[i__];
    }
L172:
    mcsrch_(n, &x[1], f, &g[1], &w[ispt + point * *n + 1], &stp, &ftol, xtol, 
	    &maxfev, &info, &nfev, &diag[1]);
    if (info == -1) {
	*iflag = 1;
	return;
    }
    if (info != 1) {
	goto L190;
    }
    nfun += nfev;

/*     COMPUTE THE NEW STEP AND GRADIENT CHANGE */
/*     ----------------------------------------- */

    npt = point * *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[ispt + npt + i__] = stp * w[ispt + npt + i__];
/* L175: */
	w[iypt + npt + i__] = g[i__] - w[i__];
    }
    ++point;
    if (point == *m) {
	point = 0;
    }

/*     TERMINATION TEST */
/*     ---------------- */

    gnorm = sqrt(ddot_(n, &g[1], &c__1, &g[1], &c__1));
    xnorm = sqrt(ddot_(n, &x[1], &c__1, &x[1], &c__1));
    xnorm = std::max(1.0,xnorm);
    if (gnorm / xnorm <= *eps) {
	finish = 1;
    }

    if (iprint[1] >= 0) {
	lb1_(&iprint[1], &iter, &nfun, &gnorm, n, m, &x[1], f, &g[1], &stp, &
		finish);
    }
    if (finish) {
	*iflag = 0;
	return;
    }
    goto L80;

/*     ------------------------------------------------------------ */
/*     END OF MAIN ITERATION LOOP. ERROR EXITS. */
/*     ------------------------------------------------------------ */

L190:
    *iflag = -1;
	std::cout << " IFLAG= -1 " << std::endl << " LINE SEARCH FAILED. SEE"
	<< " DOCUMENTATION OF ROUTINE MCSRCH" << std::endl << " ERROR RETURN"
	<< " OF LINE SEARCH: INFO= " << info << std::endl
	<< " POSSIBLE CAUSES: FUNCTION OR GRADIENT ARE INCORRECT" << std::endl
	<< " OR INCORRECT TOLERANCES" << std::endl;
	/*io___30.ciunit = lb3_1.lp;
	s_wsfe(&io___30);
	do_fio(&c__1, (char *)&info, (int)sizeof(int));
	e_wsfe();*/
    return;
L195:
    *iflag = -2;
    std::cout << " IFLAG= -2" << std::endl << " THE " << i__ << "-TH DIAGONAL ELEMENT OF THE"
	<< std::endl << " INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE" << std::endl;
    //if (lb3_1.lp > 0) {
	//io___31.ciunit = lb3_1.lp;
	//s_wsfe(&io___31);
	//do_fio(&c__1, (char *)&i__, (int)sizeof(int));
	//e_wsfe();
    //}
    return;
L196:
    *iflag = -3;
    std::cout << " IFLAG= -3" << std::endl << " IMPROPER INPUT PARAMETERS (N OR M"
    << " ARE NOT POSITIVE)" << std::endl;

    return;
} /* lbfgs_ */

/* Subroutine */ 
void lb1_(int *iprint, int *iter, int *nfun, double *gnorm, int *n, int *m, 
	  double *x, double *f, double *g, double *stp, int *finish) {
    
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;

/*     ------------------------------------------------------------- */
/*     THIS ROUTINE PRINTS MONITORING INFORMATION. THE FREQUENCY AND */
/*     AMOUNT OF OUTPUT ARE CONTROLLED BY IPRINT. */
/*     ------------------------------------------------------------- */


    /* Parameter adjustments */
    --iprint;
    --g;
    --x;

    /* Function Body */
    if (*iter == 0) {
	std::cout << "*************************************************" << std::endl;
	std::cout << "  N=" << *n << "   NUMBER OF CORRECTIONS=" << *m << std::endl
	<< "       INITIAL VALUES" << std::endl;
	std::cout << " F= " << *f << "   GNORM= " << *gnorm << std::endl;	
	if (iprint[2] >= 1) {
	    std::cout << " VECTOR X= " << std::endl;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		std::cout << x[i__] << std::endl;
	    }
	    std::cout << " GRADIENT VECTOR G= " << std::endl;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		std::cout << g[i__] << std::endl;
	    }
	}
	std::cout << "*************************************************" << std::endl;
	std::cout << "   I   NFN" << "    " << "FUNC" << "        "<< "GNORM" << 
	    "       " << "STEPLENGTH" << std::endl;
    } else {
	if (iprint[1] == 0 && (*iter != 1 && ! (*finish))) {
	    return;
	}
	if (iprint[1] != 0) {
	    if ((*iter - 1) % iprint[1] == 0 || *finish) {
		if (iprint[2] > 1 && *iter > 1) {
		    std::cout << "   I   NFN" << "    " << "FUNC" << "        "
		    << "GNORM" << "       " << "STEPLENGTH" << std::endl;
		}
		std::cout << *iter << " " << *nfun << " " << *f << " " 
		<< *gnorm << " " << *stp << std::endl;
	    } else {
		return;
	    }
	} else {
	    if (iprint[2] > 1 && *finish) {
		std::cout << "   I   NFN" << "    " << "FUNC" << "        "
		<< "GNORM" << "       " << "STEPLENGTH" << std::endl;
	    }
	    std::cout << *iter << " " << *nfun << " " << *f << " " 
	    << *gnorm << " " << *stp << std::endl;
	}
	if (iprint[2] == 2 || iprint[2] == 3) {
	    if (*finish) {
		std::cout << " FINAL POINT X= " << std::endl;
	    } else {
		std::cout << " VECTOR X= " << std::endl;
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		std::cout << x[i__] << std::endl;
	    }
	    if (iprint[2] == 3) {
		std::cout << " GRADIENT VECTOR G= " << std::endl;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    std::cout << g[i__] << std::endl;
		}
	    }
	}
	if (*finish) {
	    std::cout << " THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS." << std::endl
	    << " IFLAG = 0" << std::endl;
	}
    }


    return;
} /* lb1_ */

/* Subroutine */ 
void daxpy_(int *n, double *da, double *dx, int *incx, double *dy, int *incy) {
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, m, ix, iy, mp1;


/*     constant times a vector plus a vector. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return;
    }
    if (*da == 0.) {
	return;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] += *da * dx[i__];
/* L30: */
    }
    if (*n < 4) {
	return;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
	dy[i__] += *da * dx[i__];
	dy[i__ + 1] += *da * dx[i__ + 1];
	dy[i__ + 2] += *da * dx[i__ + 2];
	dy[i__ + 3] += *da * dx[i__ + 3];
/* L50: */
    }
    return;
} /* daxpy_ */



/*   ---------------------------------------------------------- */

double ddot_(int *n, double *dx, int *incx, double *dy, int *incy) {
    /* System generated locals */
    int i__1;
    double ret_val;

    /* Local variables */
    static int i__, m;
    static double dtemp;
    static int ix, iy, mp1;


/*     forms the dot product of two vectors. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[ix] * dy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[i__] * dy[i__];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
	dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + dx[
		i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[i__ + 
		4] * dy[i__ + 4];
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* ddot_ */

/*    ------------------------------------------------------------------ */

/*     ************************** */
/*     LINE SEARCH ROUTINE MCSRCH */
/*     ************************** */

/* Subroutine */ 
void mcsrch_(int *n, double *x, double *f, double *g, double *s, 
	     double *stp, double *ftol, double *xtol, int *maxfev, int *info, 
	     int *nfev, double *wa) {
    /* Initialized data */

    static double p5 = .5;
    static double p66 = .66;
    static double xtrapf = 4.;
    static double zero = 0.;

    /* System generated locals */
    int i__1;
    double d__1;

    /* Builtin functions */
    //int s_wsfe(), e_wsfe();

    /* Local variables */
    static double dgxm, dgym;
    static int j, infoc;
    static double finit, width, stmin, stmax;
    static int stage1;
    static double width1, ftest1, dg, fm, fx, fy;
    static int brackt;
    static double dginit, dgtest;
    static double dgm, dgx, dgy, fxm, fym, stx, sty;

/*                     SUBROUTINE MCSRCH */

/*     A slight modification of the subroutine CSRCH of More' and Thuente. */
/*     The changes are to allow reverse communication, and do not affect */
/*     the performance of the routine. */

/*     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES */
/*     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION. */

/*     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF */
/*     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF */
/*     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A */
/*     MINIMIZER OF THE MODIFIED FUNCTION */

/*          F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S). */

/*     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION */
/*     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE, */
/*     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT */
/*     CONTAINS A MINIMIZER OF F(X+STP*S). */

/*     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES */
/*     THE SUFFICIENT DECREASE CONDITION */

/*           F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S), */

/*     AND THE CURVATURE CONDITION */

/*           ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S). */

/*     IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION */
/*     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES */
/*     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH */
/*     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING */
/*     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY */
/*     SATISFIES THE SUFFICIENT DECREASE CONDITION. */

/*     THE SUBROUTINE STATEMENT IS */

/*        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA) */
/*     WHERE */

/*       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF VARIABLES. */

/*       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE */
/*         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS */
/*         X + STP*S. */

/*       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F */
/*         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S. */

/*       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE */
/*         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT */
/*         OF F AT X + STP*S. */

/*       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE */
/*         SEARCH DIRECTION. */

/*       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN */
/*         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT */
/*         STP CONTAINS THE FINAL ESTIMATE. */

/*       FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse */
/*         communication implementation GTOL is defined in a COMMON */
/*         statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE */
/*         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE */
/*         SATISFIED. */

/*       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS */
/*         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY */
/*         IS AT MOST XTOL. */

/*       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH */
/*         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse */
/*         communication implementatin they are defined in a COMMON */
/*         statement). */

/*       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION */
/*         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST */
/*         MAXFEV BY THE END OF AN ITERATION. */

/*       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS: */

/*         INFO = 0  IMPROPER INPUT PARAMETERS. */

/*         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT. */

/*         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE */
/*                   DIRECTIONAL DERIVATIVE CONDITION HOLD. */

/*         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY */
/*                   IS AT MOST XTOL. */

/*         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV. */

/*         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN. */

/*         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX. */

/*         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS. */
/*                   THERE MAY NOT BE A STEP WHICH SATISFIES THE */
/*                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS. */
/*                   TOLERANCES MAY BE TOO SMALL. */

/*       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF */
/*         CALLS TO FCN. */

/*       WA IS A WORK ARRAY OF LENGTH N. */

/*     SUBPROGRAMS CALLED */

/*       MCSTEP */

/*       FORTRAN-SUPPLIED...ABS,MAX,MIN */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983 */
/*     JORGE J. MORE', DAVID J. THUENTE */

/*     ********** */
    /* Parameter adjustments */
    --wa;
    --s;
    --g;
    --x;

    /* Function Body */
    if (*info == -1) {
	goto L45;
    }
    infoc = 1;

/*     CHECK THE INPUT PARAMETERS FOR ERRORS. */

    //if (*n <= 0 || *stp <= zero || *ftol < zero || lb3_1.gtol < zero || *xtol 
	//    < zero || lb3_1.stpmin < zero || lb3_1.stpmax < lb3_1.stpmin || *
    if (*n <= 0 || *stp <= zero || *ftol < zero ||  *xtol < zero || maxfev <= 0) {
	return;
    }

/*     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION */
/*     AND CHECK THAT S IS A DESCENT DIRECTION. */

    dginit = zero;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dginit += g[j] * s[j];
/* L10: */
    }
    if (dginit >= zero) {
	std::cout << "  THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION" << std::endl;
	return;
    }

/*     INITIALIZE LOCAL VARIABLES. */

    brackt = 0;
    stage1 = 1;
    *nfev = 0;
    finit = *f;
    dgtest = *ftol * dginit;
    //width = lb3_1.stpmax - lb3_1.stpmin;
    width = 1.0e20 - 1.0e-20;
    width1 = width / p5;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa[j] = x[j];
    }

/*     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP, */
/*     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP. */
/*     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP, */
/*     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF */
/*     THE INTERVAL OF UNCERTAINTY. */
/*     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP, */
/*     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP. */

    stx = zero;
    fx = finit;
    dgx = dginit;
    sty = zero;
    fy = finit;
    dgy = dginit;

/*     START OF ITERATION. */

L30:

/*        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND */
/*        TO THE PRESENT INTERVAL OF UNCERTAINTY. */

    if (brackt) {
	stmin = std::min(stx,sty);
	stmax = std::max(stx,sty);
    } else {
	stmin = stx;
	stmax = *stp + xtrapf * (*stp - stx);
    }

/*        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN. */

    *stp = std::max(*stp,1.0e-20);//lb3_1.stpmin)
    *stp = std::min(*stp,1.0e20);//lb3_1.stpmax)

/*        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET */
/*        STP BE THE LOWEST POINT OBTAINED SO FAR. */

    if (brackt && (*stp <= stmin || *stp >= stmax) || *nfev >= *maxfev - 1 || 
	    infoc == 0 || brackt && stmax - stmin <= *xtol * stmax) {
	*stp = stx;
    }

/*        EVALUATE THE FUNCTION AND GRADIENT AT STP */
/*        AND COMPUTE THE DIRECTIONAL DERIVATIVE. */
/*        We return to main program to obtain F and G. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x[j] = wa[j] + *stp * s[j];
    }
    *info = -1;
    return;

L45:
    *info = 0;
    ++(*nfev);
    dg = zero;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dg += g[j] * s[j];
    }
    ftest1 = finit + *stp * dgtest;

/*        TEST FOR CONVERGENCE. */

    if (brackt && (*stp <= stmin || *stp >= stmax) || infoc == 0) {
	*info = 6;
    }
    //if (*stp == lb3_1.stpmax && *f <= ftest1 && dg <= dgtest) {
    if (*stp == 1.0e20 && *f <= ftest1 && dg <= dgtest) {
	*info = 5;
    }
    //if (*stp == lb3_1.stpmin && (*f > ftest1 || dg >= dgtest)) {
    if (*stp == 1.0e-20 && (*f > ftest1 || dg >= dgtest)) {
	*info = 4;
    }
    if (*nfev >= *maxfev) {
	*info = 3;
    }
    if (brackt && stmax - stmin <= *xtol * stmax) {
	*info = 2;
    }
    //if (*f <= ftest1 && fabs(dg) <= lb3_1.gtol * (-dginit)) {
    if (*f <= ftest1 && fabs(dg) <= 9e-1 * (-dginit)) {
	*info = 1;
    }

/*        CHECK FOR TERMINATION. */

    if (*info != 0) {
	return;
    }

/*        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED */
/*        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE. */

    //if (stage1 && *f <= ftest1 && dg >= std::min(*ftol,lb3_1.gtol) * dginit) {
    if (stage1 && *f <= ftest1 && dg >= std::min(*ftol,0.9) * dginit) {
	stage1 = 0;
    }

/*        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF */
/*        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED */
/*        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE */
/*        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN */
/*        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT. */

    if (stage1 && *f <= fx && *f > ftest1) {

/*           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES. */

	fm = *f - *stp * dgtest;
	fxm = fx - stx * dgtest;
	fym = fy - sty * dgtest;
	dgm = dg - dgtest;
	dgxm = dgx - dgtest;
	dgym = dgy - dgtest;

/*           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY */
/*           AND TO COMPUTE THE NEW STEP. */

	mcstep_(&stx, &fxm, &dgxm, &sty, &fym, &dgym, stp, &fm, &dgm, &brackt,
		 &stmin, &stmax, &infoc);

/*           RESET THE FUNCTION AND GRADIENT VALUES FOR F. */

	fx = fxm + stx * dgtest;
	fy = fym + sty * dgtest;
	dgx = dgxm + dgtest;
	dgy = dgym + dgtest;
    } else {

/*           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY */
/*           AND TO COMPUTE THE NEW STEP. */

	mcstep_(&stx, &fx, &dgx, &sty, &fy, &dgy, stp, f, &dg, &brackt, &
		stmin, &stmax, &infoc);
    }

/*        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE */
/*        INTERVAL OF UNCERTAINTY. */

    if (brackt) {
	if ((d__1 = sty - stx, fabs(d__1)) >= p66 * width1) {
	    *stp = stx + p5 * (sty - stx);
	}
	width1 = width;
	width = (d__1 = sty - stx, fabs(d__1));
    }

/*        END OF ITERATION. */

    goto L30;

/*     LAST LINE OF SUBROUTINE MCSRCH. */

} /* mcsrch_ */

/* Subroutine */ 
void mcstep_(double *stx, double *fx, double *dx, double *sty, double *fy, double *dy, 
	     double *stp, double *fp, double *dp, int *brackt, double  *stpmin, 
	     double *stpmax, int *info) {
    /* System generated locals */
    double d__1, d__2, d__3;

    /* Local variables */
    static double sgnd, stpc, stpf, stpq, p, q, gamma, r__, s, theta;
    static int bound;


/*     SUBROUTINE MCSTEP */

/*     THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR */
/*     A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR */
/*     A MINIMIZER OF THE FUNCTION. */

/*     THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION */
/*     VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS */
/*     ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE */
/*     DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A */
/*     MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY */
/*     WITH ENDPOINTS STX AND STY. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT, */
/*                        STPMIN,STPMAX,INFO) */

/*     WHERE */

/*       STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP, */
/*         THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED */
/*         SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION */
/*         OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE */
/*         SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY. */

/*       STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP, */
/*         THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF */
/*         THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE */
/*         UPDATED APPROPRIATELY. */

/*       STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP, */
/*         THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP. */
/*         IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE */
/*         BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP. */

/*       BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER */
/*         HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED */
/*         THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER */
/*         IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE. */

/*       STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER */
/*         AND UPPER BOUNDS FOR THE STEP. */

/*       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS: */
/*         IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED */
/*         ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE */
/*         INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS. */

/*     SUBPROGRAMS CALLED */

/*       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983 */
/*     JORGE J. MORE', DAVID J. THUENTE */

    *info = 0;

/*     CHECK THE INPUT PARAMETERS FOR ERRORS. */

    if (*brackt && (*stp <= std::min(*stx,*sty) || *stp >= std::max(*stx,*sty)) || *dx *
	     (*stp - *stx) >= 0.0 || *stpmax < *stpmin) {
	return;
    }

/*     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN. */

    sgnd = *dp * (*dx / fabs(*dx));

/*     FIRST CASE. A HIGHER FUNCTION VALUE. */
/*     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER */
/*     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN, */
/*     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN. */

    if (*fp > *fx) {
	*info = 1;
	bound = 1;
	theta = (*fx - *fp) * 3 / (*stp - *stx) + *dx + *dp;
/* Computing MAX */
	d__1 = fabs(theta), d__2 = fabs(*dx), d__1 = std::max(d__1,d__2), d__2 = fabs(*dp);
	s = std::max(d__1,d__2);
/* Computing 2nd power */
	d__1 = theta / s;
	gamma = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
	if (*stp < *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dx + theta;
	q = gamma - *dx + gamma + *dp;
	r__ = p / q;
	stpc = *stx + r__ * (*stp - *stx);
	stpq = *stx + *dx / ((*fx - *fp) / (*stp - *stx) + *dx) / 2 * (*stp - 
		*stx);
	if ((d__1 = stpc - *stx, fabs(d__1)) < (d__2 = stpq - *stx, fabs(d__2)))
		 {
	    stpf = stpc;
	} else {
	    stpf = stpc + (stpq - stpc) / 2;
	}
	*brackt = 1;

/*     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF */
/*     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC */
/*     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP, */
/*     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN. */

    } else if (sgnd < (float)0.) {
	*info = 2;
	bound = 0;
	theta = (*fx - *fp) * 3 / (*stp - *stx) + *dx + *dp;
/* Computing MAX */
	d__1 = fabs(theta), d__2 = fabs(*dx), d__1 = std::max(d__1,d__2), d__2 = fabs(*dp);
	s = std::max(d__1,d__2);
/* Computing 2nd power */
	d__1 = theta / s;
	gamma = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
	if (*stp > *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dp + theta;
	q = gamma - *dp + gamma + *dx;
	r__ = p / q;
	stpc = *stp + r__ * (*stx - *stp);
	stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
	if ((d__1 = stpc - *stp, fabs(d__1)) > (d__2 = stpq - *stp, fabs(d__2)))
		 {
	    stpf = stpc;
	} else {
	    stpf = stpq;
	}
	*brackt = 1;

/*     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE */
/*     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES. */
/*     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY */
/*     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC */
/*     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE */
/*     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO */
/*     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP */
/*     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN. */

    } else if (fabs(*dp) < fabs(*dx)) {
	*info = 3;
	bound = 1;
	theta = (*fx - *fp) * 3 / (*stp - *stx) + *dx + *dp;
/* Computing MAX */
	d__1 = fabs(theta), d__2 = fabs(*dx), d__1 = std::max(d__1,d__2), d__2 = fabs(*dp);
	s = std::max(d__1,d__2);

/*        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND */
/*        TO INFINITY IN THE DIRECTION OF THE STEP. */

/* Computing MAX */
/* Computing 2nd power */
	d__3 = theta / s;
	d__1 = 0., d__2 = d__3 * d__3 - *dx / s * (*dp / s);
	gamma = s * sqrt((std::max(d__1,d__2)));
	if (*stp > *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dp + theta;
	q = gamma + (*dx - *dp) + gamma;
	r__ = p / q;
	if (r__ < (float)0. && gamma != (float)0.) {
	    stpc = *stp + r__ * (*stx - *stp);
	} else if (*stp > *stx) {
	    stpc = *stpmax;
	} else {
	    stpc = *stpmin;
	}
	stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
	if (*brackt) {
	    if ((d__1 = *stp - stpc, fabs(d__1)) < (d__2 = *stp - stpq, fabs(
		    d__2))) {
		stpf = stpc;
	    } else {
		stpf = stpq;
	    }
	} else {
	    if ((d__1 = *stp - stpc, fabs(d__1)) > (d__2 = *stp - stpq, fabs(
		    d__2))) {
		stpf = stpc;
	    } else {
		stpf = stpq;
	    }
	}

/*     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE */
/*     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES */
/*     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP */
/*     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN. */

    } else {
	*info = 4;
	bound = 0;
	if (*brackt) {
	    theta = (*fp - *fy) * 3 / (*sty - *stp) + *dy + *dp;
/* Computing MAX */
	    d__1 = fabs(theta), d__2 = fabs(*dy), d__1 = std::max(d__1,d__2), d__2 = 
		    fabs(*dp);
	    s = std::max(d__1,d__2);
/* Computing 2nd power */
	    d__1 = theta / s;
	    gamma = s * sqrt(d__1 * d__1 - *dy / s * (*dp / s));
	    if (*stp > *sty) {
		gamma = -gamma;
	    }
	    p = gamma - *dp + theta;
	    q = gamma - *dp + gamma + *dy;
	    r__ = p / q;
	    stpc = *stp + r__ * (*sty - *stp);
	    stpf = stpc;
	} else if (*stp > *stx) {
	    stpf = *stpmax;
	} else {
	    stpf = *stpmin;
	}
    }

/*     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT */
/*     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE. */

    if (*fp > *fx) {
	*sty = *stp;
	*fy = *fp;
	*dy = *dp;
    } else {
	if (sgnd < 0.0) {
	    *sty = *stx;
	    *fy = *fx;
	    *dy = *dx;
	}
	*stx = *stp;
	*fx = *fp;
	*dx = *dp;
    }

/*     COMPUTE THE NEW STEP AND SAFEGUARD IT. */

    stpf = std::min(*stpmax,stpf);
    stpf = std::max(*stpmin,stpf);
    *stp = stpf;
    if (*brackt && bound) {
	if (*sty > *stx) {
/* Computing MIN */
	    d__1 = *stx + (*sty - *stx) * 0.66;
	    *stp = std::min(d__1,*stp);
	} else {
/* Computing MAX */
	    d__1 = *stx + (*sty - *stx) * 0.66;
	    *stp = std::max(d__1,*stp);
	}
    }
    return;

} /* mcstep_ */

} // MTK++ namespace
