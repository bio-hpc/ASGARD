/*//////////////////////////////////////////////////////////////////////////
// 
// elsize.c Calculates the electrostatic size, A, of a molecule using 
// atomic coordinates and radii  from its PQR (or simply xyz) file.
// Author: Grigori Sigalov <sigalov@vt.edu>, <greg_sigalov@yahoo.com>
// Conversion to ANSI-C(89): Richard Shadrach
// Algorithm Details are in: 
//
// "Analytical Linearized Poisson--Boltzmann Approach 
//   for Molecular Dynamics Applications", by Sigalov, Fenley and Onufriev
//  J. Chem. Phys., 124, 124902  (2006)
//
////////////////////////////////////////////////////////////////////////////

// Brief algorithm: First, the moments of inertia of the molecule are
// calculated. If the program is run without second parameter or, which is the
// same, with second parameter being "-det", the electrostatic size is found
// from the third invariant (determinant) of the tensor of inertia.  Otherwise
// some more work has to be done. The _principal_ moments of inertia are found
// by solving the cubic equation det(I - lambda E) = 0. The semiaxes of the
// ellipsoid that has same principal moments of inertia and mass as the
// original molecule are found. Then, there are two possibilities: option
// "-ell" leads to calculation of the exact electrostatic size, A_ell, of the
// effective ellipoid just found, which involves numerical summation by
// Simpson method of an elliptic intergal of the first kind (to handle the
// infinite limit, the integral is transformed to map infinity to zero). With
// option "-elf", the value A_elf, which is an approximation to A_ell
// expressed using elementary functions only, is calculated. Option "-abc"
// prints the semiaxes of the effective ellipsoid.  Option "-tab" prints all
// of the above into a table; option -hea prints table with a header.
// Finally, option -deb prints additionally some extra information. Any of
// these options, if used at all, MUST follow the input file name. 

// By default, the input file is in PQR format. If you have only PDB file
// available, convert it to a plain XYZ format and use an extra option,
// "-xyz". In this case all atomic sizes are set to the same value,
// DEFAULTATOMSIZE.

// Shortly, the following hierarchy can be built:

// 1. -det (default) gives a reasonable estimate of the electrostatic size. In
// most cases, you will want nothing more. Moreover, A_det has simple
// analytical derivatives, so if you are going to use derivatives of the
// electrostatic size, you should choose A_det. This option involves
// calculation of the moments of inertia of the molecule, which is a very
// straighforward procedure.  

// 2. -elf adds some accuracy if the molecule's shape is close to an
// ellipsoid, which rarely happens, at the price of a cubic equation to solve.

// 3. -ell gives an exact solution for an ellipsoidal molecule, and involves
// numerical integration in addition to the cubic equation. Options -abc,
// -tab, -hea, and -deb include -ell and therefore take the same amount of
// calculation.
//
// Previous versions allowed for 64 bit ints.  In order to conform with ANSI C
// (C89), this must be taken out. It will allow for ~4 billion atoms, which 
// should be more than sufficient*/


#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/*ANSI C does not require M_PI to be defined.*/
#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif

#define DEFAULTATOMSIZE 1.5
#define ACCURACY 1.E-09 /*convergence of Simpson's method, stopping condition*/
#define MAXLINESIZE 1024 /* For the input file */
#define TOL 1.E-12 /* values < TOL are treated as zeros */

/* this setting (default) suppresses all warnings; 
   all with key -deb if unsure about your input file */
int debugEnabled = 0;   

unsigned long GetNumberOfAtoms(const char *fileName, int XYZ_format);
unsigned long ReadAtomicCoordinates(const char *fileName, int XYZ_format, 
  unsigned long natoms, double *x, double *y, double *z, double *r2, double *r3);
double EllipticIntegral(double a, double b, double c, double accuracy);
double A_elementary_functions(double a, double b, double c);

/* a couple of auxiliary functions used in numerical integration procedure */
double function1(double a, double b, double c, double x);
double function2(double a, double b, double c, double x);

void usage();
#ifdef __GNUC__
static void debug(const char*, ...) __attribute__ ((format(printf, 1, 2)));
static void fatal(const char*, ...) __attribute__ ((format(printf, 1, 2)))
                                    __attribute__ ((__noreturn__));
#else
static void debug(const char*, ...);
static void fatal(const char*, ...);
#endif /* __GNUC__ */

int main( int argc, char ** argv )
{
    int option = -1; 
    int header = 0;
    int XYZ_format = 0; /*by default, input file has PQR format, not XYZ*/

    unsigned long natoms, i, natoms_check;
    double  *x, *y, *z, *r2, *r3; /* Position & radius squared and cubed */ 
                                  /* r^2 for the moment about the centroid */
                                  /* r^3 (an estimate of atomic mass) */
    double I11 = 0., I12 = 0., I13 = 0., I22 = 0., I23 = 0., I33 = 0.;
    double molecule_mass = 0.;
    double Ixx, Iyy, Izz;
    double A, B, C, P, Q;
    double a, b, c, d;
    double A_det, A_ell, A_elf;

    if (argc < 2) /* there is NO argument argv[1] */
    {
        usage();
    }

    /* First parameter is the input file (hopefully).  From the 2nd on,
       we accept the first valid argument as the argument to use.  If
       "-xyz" ever appears, use that as well. */
    for (i = 2; i < argc; i++) 
    {
        if (option == -1) {
            if (!strcmp(argv[i], "-det")) {
                option = 0;
            } else if (!strcmp(argv[i], "-ell")) {
                option = 1;
            } else if (!strcmp(argv[i], "-elf")) {
                option = 2;
            } else if (!strcmp(argv[i], "-abc")) {
                option = 3;
            } else if (!strcmp(argv[i], "-tab")) {
                option = 4;
            } else if (!strcmp(argv[i], "-hea")) {
                option = 4; header = 1;
            } else if (!strcmp(argv[i], "-deb")) {
                option = 4; debugEnabled = 1;
            }
            else {
            	fprintf(stderr, "WARNING: Invalid option %s provided.  "
            	  "Ignoring it and continuing execution.\n", argv[i]);
            }
        }
				else if (!strcmp(argv[i], "-xyz"))
        {
            XYZ_format = 1;
        }
				else {
					  fprintf(stderr,"WARNING: The option \"%s\" is either invalid or \n"
              "         a mutually exlcusive option was already provided.\n"
					    "         Ignoring it and continuing execution.\n", argv[i]);
        }
    }
    if (option == -1) {
        option = 0; /*Default is det*/
    }

    /* Get the number of atoms, allocate the proper size, 
       then read in data */
    natoms = GetNumberOfAtoms(argv[1], XYZ_format);
    x  = (double *)calloc(natoms, sizeof(double));
    y  = (double *)calloc(natoms, sizeof(double));
    z  = (double *)calloc(natoms, sizeof(double));
    r2 = (double *)calloc(natoms, sizeof(double));
    r3 = (double *)calloc(natoms, sizeof(double));
    natoms_check = ReadAtomicCoordinates(argv[1], XYZ_format, natoms, x, 
      y, z, r2, r3);

    if ( natoms_check != natoms )
    {
        fatal("First time %lu data lines were read, now it's %lu\n"
          "... exiting ...\n", natoms, natoms_check);
    }
    
    { /*Finds molecule center of mass, centers it, and moments of interia*/
        double xav = 0;
        double yav = 0;
        double zav = 0;

        for (i = 0; i<natoms; i++)
        {
            /* Calculating molecule mass and center of mass */
            xav += x[i] * r3[i]; /* atom's "mass" is simply r^3 */
            yav += y[i] * r3[i]; 
            zav += z[i] * r3[i];
            molecule_mass += r3[i];
        }
    
        xav /= molecule_mass; /* now it's the center */
        yav /= molecule_mass;
        zav /= molecule_mass;
    
        debug("Molecule %s: %lu atoms, center (%.3f, %.3f, %.3f)\n",
          argv[1], natoms, xav, yav, zav);

        for (i=0; i<natoms; i++)
        {
            double atom_mass = r3[i];
            double x2, y2, z2, atoms_MI;
            
            /* moving the center to 0 */
            x[i] -= xav;
            y[i] -= yav;
            z[i] -= zav;

            /* Calculate the moments of inertia */

            /* half atom's moment of ineria about
               diameter (will be added twice!) */
            atoms_MI = r2[i] / 5.;
        
            x2 = x[i] * x[i] + atoms_MI;
            y2 = y[i] * y[i] + atoms_MI;
            z2 = z[i] * z[i] + atoms_MI;
        
            I11 += atom_mass * (y2 + z2);
            I22 += atom_mass * (z2 + x2);
            I33 += atom_mass * (x2 + y2);
            
            /*atoms moments do not contribute because of symmetry*/
            I12 -= atom_mass * x[i] * y[i]; 
            I13 -= atom_mass * x[i] * z[i];
            I23 -= atom_mass * y[i] * z[i];
        }
    }
    

    debug("Original tensor of inertia: %12.4g %12.4g %12.4g\n",I11,I12,I13);
    debug("                            %12.4g %12.4g %12.4g\n",I12,I22,I23);
    debug("                            %12.4g %12.4g %12.4g\n",I13,I23,I33);
    
    if ( option == 0 || option == 4 )   
    {
        double det_I = I11 * I22 * I33 + 2. * I12 * I23 * I13 
          - I11 * I23 * I23  - I22 * I13 * I13 - I33 * I12 * I12;
        if (det_I <= 0.)
        {
            fatal("Problem with the determinant of the tensor of"
              " intertia: Det I = %g\n", det_I);
        }
        A_det = sqrt( 2.5 / molecule_mass ) * pow( det_I, 1./6. );
    }
 
    if ( option == 0 ) /* by default, no need to calculate anything else */
    {
        printf("%.3f\n", A_det);
        return 0;
    }

    /* Now let's find the coefficients of the cubic equation for 
       eigenvalues of the matrix I */
    
    /* To avoid artifacts like complex roots for nearly symmetrical, 
       molecules small enough non-diagonal moments should be set to
       exact zeros */
    {
        double TrI = I11 + I22 + I33; /* Trace I */

        debug("Checking if non-diagonal moments are negligible: "
          "I12/Tr I = %g, I13/Tr I = %g, I23/Tr I = %g\n", 
          I12/TrI, I13/TrI, I23/TrI);

        if ( fabs( I12 / TrI ) < TOL ) I12 = 0.;
        if ( fabs( I13 / TrI ) < TOL ) I13 = 0.;
        if ( fabs( I23 / TrI ) < TOL ) I23 = 0.;
    }
    
    A = - (I11 + I22 + I33);
    B = I11 * I22 + I22 * I33 + I33 * I11 - I12 * I12 - I13 * I13 
      - I23 * I23;
    C = - (I11 * I22 * I33 + 2. * I12 * I23 * I13 - 
        I11 * I23 * I23 - I22 * I13 * I13 - I33 * I12 * I12);

    debug("Original cubic equation's coefficients: A = %g, B = %g, C = %g\n",
      A, B, C);

    /* reducing the cubic equation to an incomplete one */

    P = - A * A / 3. + B;
    Q = A * A * A * 2. / 27. - A * B / 3. + C;  
    
    debug("Reduced cubic equation's coefficients: P = %g, Q = %g\n", P, Q);
    
    if (P >= 0. && Q != 0 && fabs(P/Q) < TOL)
        P = 0.;

    /* if molecule is exactly spherical, extra care should 
       be taken because both P and Q == 0 */

    if ( I12 == 0. && I13 == 0. && I23 == 0. )
    {
        P = 0.;
        Q = 0.;
    }

    if (P > 0.)
    {
        fatal("Problem with the cubic equation: at least some of the"
          "roots are not real!\nOriginal cubic equation's coefficients:"
          " A = %g, B = %g, C = %g\nReduced cubic equation's"
          " coefficients: P = %g, Q = %g\n", A, B, C, P, Q);
    }
    else if (P == 0.) /* degeneration a = b = c = -Q^{1/3} */
    {
        if ( Q > 0 )
            Ixx = Iyy = Izz = - pow(Q, 1./3.) - A/3.;
        else if ( Q < 0 )
            Ixx = Iyy = Izz = pow(-Q, 1./3.) - A/3.;
        else
            Ixx = Iyy = Izz = - A/3.;
    }
    else
    {
        double paux = 2. * sqrt( - P / 3. );
        double cos_alpha = - 4. * Q / paux / paux / paux;
        double alpha;

        if (cos_alpha < -1 && cos_alpha > -1 - TOL) {
            cos_alpha = -1;
        }
        if (cos_alpha > 1 && cos_alpha < 1 + TOL) {
            cos_alpha = 1;
        }
        if (cos_alpha <= -1-TOL || cos_alpha >= 1+TOL) {
            fatal("cos_alpha was %g, can't take the arccos.",
              cos_alpha);
        }
        alpha = acos(cos_alpha);
    
        Ixx = paux * cos(alpha/3.) - A/3.;
        Iyy = paux * cos(alpha/3. + M_PI/3.) - A/3.;
        Izz = paux * cos(alpha/3. - M_PI/3.) - A/3.;
    }

    /* finally, let's find the axes! */
    a = sqrt( 2.5 * (-Ixx + Iyy + Izz) / molecule_mass );
    b = sqrt( 2.5 * ( Ixx - Iyy + Izz) / molecule_mass );
    c = sqrt( 2.5 * ( Ixx + Iyy - Izz) / molecule_mass );

    if (b > a) { d = a; a = b; b = d; } /* to make sure that a >= b >= c */
    if (c > b) { d = b; b = c; c = d; } 
    if (b > a) { d = a; a = b; b = d; } 

    if ( c <= 0. )
    {
        fatal("Problem with the semiaxes: a = %g, b = %g, c = %g\n",
          a, b, c);
    }

    if ( option == 3 )  /* print a, b, c */
    {
        printf("%.3f %9.3f %9.3f   ", a, b, c);
        printf("\n");
        return 0;
    }

    if ( option == 1 || option == 4 )
    {
        if ( b - c < TOL ) /*b==c, its exact in this particular case*/
            A_ell = A_elementary_functions(a, b, c);
        else
            A_ell = 2. / EllipticIntegral(a, b, c, ACCURACY);
    }

    if ( option == 1 )
    {
        printf("%.3f", A_ell);
        printf("\n");
        return 0;
    }

    if ( option == 2 || option == 4 )
        A_elf = A_elementary_functions(a, b, c);

    if ( option == 2 )
    {
        printf("%.3f", A_elf);
        printf("\n");
        return 0;
    }
    
    if (header || debugEnabled) {
        printf("Filename                 a        b        c"
          "       A_ell    A_elf    A_det\n");
    }

    printf("%-20s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f",
      argv[1], a, b, c, A_ell, A_elf, A_det);
    printf("\n");

    return 0;
}

void usage() {
    fprintf(stderr,
      "Usage: esize your_structure.pqr [ arg ]"
      "The second argument is optional: \n"
      "    -det (default) gives A_det\n"
      "    -ell gives A_ell (elliptic integral)\n"
      "    -elf gives A_elf (elementary functions approximation to\n"
      "                      A_ell, normally less than 0.1A apart)\n"
      "    -abc prints a, b, c (semiaxes of the effective ellipsoid,\n"
      "                         just out of curiousity)\n");
    fprintf(stderr,
      "    -tab prints PQR file name and all of the above into a\n"
      "                                     table without header\n"
      "    -hea prints same table as -tab but with a header"
      "    -deb prints same as -tab with some extra (debugging) information"
      "    -xyz uses a file containing only XYZ coordinates as input.\n");
    exit(1);
}

static void fatal(const char* fmt, ...)
{
    va_list ap;

    assert(fmt && *fmt);

    va_start(ap, fmt);
    fputs(" ** Error ** : ", stderr);
    vfprintf(stderr, fmt, ap);
    fflush(stderr);
    va_end(ap);

    exit(EXIT_FAILURE);
}

static void debug(const char* fmt, ...)
{
    assert(fmt && *fmt);

    if (debugEnabled != 0) {
        va_list ap;

        va_start(ap, fmt);
        fputs("DEBUG : ", stdout);
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
        va_end(ap);
    }
}

/*---------------------------------------------------------------------------*/
/*Calculates:
   $\int_0^\infty \frac{d\theta}{\sqrt{(a^2+\theta)(b^2+\theta)(c^2+\theta)}}$
  by 1-4-1 Simpson method */

double EllipticIntegral(double a, double b, double c, double accuracy)
{
    double d;
    int i, n, n0=16;
    double step, lim, sum1=0., sum2=0., total = 0., prev=1.;

    if ( a == b && b == c )
        return a;

    d = c;  /* wild guess */
        
    /* calculating 1st integral from 0 to d/c */

    lim = d / c;
    for (n=n0; fabs(sum1 - prev) > accuracy; n*=2)
    {
        double sum_aux;
        prev = sum1;

        step = lim / n;
        sum1 = function1(a, b, c, 0.) + function1(a, b, c, lim);

        sum_aux = 0.;
        for (i = 1; i < n; i+=2)
            sum_aux += function1(a, b, c, i * step);
        sum1 += 4. * sum_aux;

        sum_aux = 0.;
        for (i = 2; i < n; i+=2)
            sum_aux += function1(a, b, c, i * step);
        sum1 += 2. * sum_aux;

        /* because I transformed the integral
           and 1/c appeared in front of it */
        sum1 *= step / 3. / c; 

    }

    /* calculating 2nd integral from d to \infty. 
       It is transformed to integral from 0 to 1/d,
       but the first part from 0 to some eps must be
       calculated analytically. */
        
    lim = 1./d;
    prev=1.;    /* different enough from sum2 = 0 (initially) */

    for (n=n0; fabs(sum2 - prev) > accuracy; n*=2)
    {
        double sum_aux;
        prev = sum2;
        
        step = lim / n;
        sum2 = function2(a, b, c, 0.) + function2(a, b, c, lim);

        sum_aux = 0.;
        for (i = 1; i < n; i+=2)
            sum_aux += function2(a, b, c, i * step);
        sum2 += 4. * sum_aux;

        sum_aux = 0.;
        for (i = 2; i < n; i+=2)
            sum_aux += function2(a, b, c, i * step);
        sum2 += 2. * sum_aux;

        sum2 *= 2. * step / 3.;
    }

    total = sum1 + sum2;
    
    return total;   
}

/*---------------------------------------------------------------------------*/

double function1(double a, double b, double c, double x)
{
    double ac = a / c;
    double bc = b / c;
    double f = (x + ac*ac) * (x + bc*bc) * (x + 1.);
    f = 1. / sqrt(f);
    return f;
}

/*---------------------------------------------------------------------------*/

double function2(double a, double b, double c, double x)
{
    double xx = x * x;
    double f = (1. + a * a * xx) * (1. + b * b * xx) * (1. + c * c * xx);
    f = 1. / sqrt(f);
    return f;
}

/*---------------------------------------------------------------------------*/

double A_elementary_functions(double a, double b, double c)
{
    double gamma, A_elf;
    if ( a - b < TOL && b - c < TOL )   /* means that a == b && b == c */
        return a;
    
    gamma = sqrt( 1. - (b+c)*(b+c)/a/a/4. );
    A_elf = 2. * a * gamma / log( (1+gamma) / (1-gamma) );
    return A_elf;
}

/*---------------------------------------------------------------------------*/

/* Opens input file, scans is to find out the number of 
   valid data lines, and closes it. No data are actually read here.
   If flag XYZ_format is set, every non-empty line is considered a data line */

unsigned long GetNumberOfAtoms(const char *fileName, int XYZ_format)
{
    FILE *fp = fopen(fileName, "r");
    unsigned long i, j, j_prev=0, natoms=0;
    char line[MAXLINESIZE];

    if (fp == NULL)
    {
        fatal("Error opening file %s to read...\n...exiting ...\n", 
          fileName);
    }

    /* count number of atoms and HETATMs in the file */
    /* and scan for inconsistencies                  */

    /* i is the line count; not every line is a data line! */
    for (i=0; !feof(fp); i++)
    {
        fgets(line, MAXLINESIZE, fp);

        /* catches that trailing newline */
        if (feof(fp)) break;

        if ( XYZ_format )
        {
            natoms++;
        }     
        else
        {
            char firstField[MAXLINESIZE];
            sscanf(line, "%s %lu", firstField, &j);
            if (strcmp(firstField, "ATOM") == 0 || 
              strcmp(firstField, "HETATM") == 0)
            {
                /* shouldn't check it for the first line */
                if ( j_prev+1 != j && i && debugEnabled )
                {
                    fprintf(stderr, "Error in %s:\n there "
                      "appears to be an inconsistency in "
                      "atom numbering between lines %lu"
                      "and %lu\n", fileName, i, i+1);
                } /* this inconsistency may exist for a good reason */

                j_prev = j;
                natoms++;
            }
            else
            {
                debug("Ignoring a line, tag is %s, line %lu\n",
                  firstField, i);
            }
        }
    }
    fclose(fp);
    return natoms;
}

/*---------------------------------------------------------------------------*/

/* Actually reads the coordinates of atoms and, 
   if it's a PQR format file, atomic sizes */
    
unsigned long ReadAtomicCoordinates(const char *fileName, int XYZ_format,
  unsigned long natoms, double *x, double *y, double *z, double *r2, double *r3)
{
    FILE *fp = fopen(fileName, "r");
    unsigned long i;
    double r;
    char line[MAXLINESIZE];
    char junk[MAXLINESIZE];

    if (fp == NULL)
    {
        fatal("Error opening file %s to read...\n...exiting ...\n",
          fileName);
    }

    for (i = 0; !feof(fp) && i<natoms;)
    {
        if ( XYZ_format )
        {
            fscanf(fp, "%lf %lf %lf", &x[i], &y[i], &z[i]);
            /* atoms are same size if it's unknown*/
            r2[i] = DEFAULTATOMSIZE * DEFAULTATOMSIZE;
            r3[i] = r2[i] * DEFAULTATOMSIZE;            
            i++;
        }
        else
        {
            fgets(line, MAXLINESIZE, fp);
            sscanf(line, "%s", junk);
            if (strcmp(junk, "ATOM") == 0 || 
              strcmp(junk, "HETATM") == 0)
            {
                /* Four data fields not needed */
                sscanf(line,"%s %s %s %s %s %lf %lf %lf %s %lf",
                  junk, junk, junk, junk, junk, &x[i], &y[i],
                  &z[i], junk, &r);
                  
                r2[i] = r * r;
                r3[i] = r2[i] * r;
                i++;
            }
        }
    }
    fclose(fp);
    return i;   /* to check if the same number of data lines were read */
}
