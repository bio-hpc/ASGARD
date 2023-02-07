/* elsize.cc
   
   Calculates the electrostatic size, A, of a molecule using 
   atomic coordinates and radii from its PQR (or simply xyz) file.
   Author: Grigori Sigalov <sigalov@vt.edu>, <greg_sigalov@yahoo.com>
  
   Algorithm Details are in: 
  
   "Analytical Linearized Poisson--Boltzmann Approach 
    for Molecular Dynamics Applications", by Sigalov, Fenley and Onufriev
    J. Chem. Phys., 124, 124902  (2006)                                    */
 
/* Brief algorithm: First, the moments of inertia of the molecule are
   calculated. If the program is run without second parameter or, which is the
   same, with second parameter being "-det", the electrostatic size is found
   from the third invariant (determinant) of the tensor of inertia.  Otherwise
   some more work has to be done. The _principal_ moments of inertia are found
   by solving the cubic equation det(I - lambda E) = 0. The semiaxes of the
   ellipsoid that has same principal moments of inertia and mass as the
   original molecule are found. Then, there are two possibilities: option
   "-ell" leads to calculation of the exact electrostatic size, A_ell, of the
   effective ellipoid just found, which involves numerical summation by
   Simpson method of an elliptic intergal of the first kind (to handle the
   infinite limit, the integral is transformed to map infinity to zero). With
   option "-elf", the value A_elf, which is an approximation to A_ell
   expressed using elementary functions only, is calculated. Option "-abc"
   prints the semiaxes of the effective ellipsoid.  Option "-tab" prints all
   of the above into a table; option -hea prints table with a header.
   Finally, option -deb prints additionally some extra information. Any of
   these options, if used at all, MUST follow the input file name. 
    
   By default, the input file is in PQR format. If you have only PDB file
   available, convert it to a plain XYZ format and use an extra option,
   "-xyz". In this case all atomic sizes are set to the same value,
   DEFAULTATOMSIZE.
   
   Shortly, the following hierarchy can be built:
   
   1. -det (default) gives a reasonable estimate of the electrostatic size. In
   most cases, you will want nothing more. Moreover, A_det has simple
   analytical derivatives, so if you are going to use derivatives of the
   electrostatic size, you should choose A_det. This option involves
   calculation of the moments of inertia of the molecule, which is a very
   straighforward procedure.
   
   2. -elf adds some accuracy if the molecule's shape is close to an
   ellipsoid, which rarely happens, at the price of a cubic equation to solve.
   
   3. -ell gives an exact solution for an ellipsoidal molecule, and involves
   numerical integration in addition to the cubic equation. Options -abc,
   -tab, -hea, and -deb include -ell and therefore take the same amount of
   calculation. */

/* Converted to C by vbabin-at-ncsu-dot-edu */

#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* makes gcc happier when it is invoked with -pedantic -ansi */
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif /* M_PI */
/* macros */

#define DEFAULT_ATOM_SIZE 1.5
#define TINY 1.0E-12

/* variables */

static int  do_debug  = 0;

/* prototypes */

#ifdef __GNUC__
static void debug(const char*, ...) __attribute__ ((format(printf, 1, 2)));
static void fatal(const char*, ...) __attribute__ ((format(printf, 1, 2)))
                                    __attribute__ ((__noreturn__));
#else
static void debug(const char*, ...);
static void fatal(const char*, ...);
#endif /* __GNUC__ */

static void usage(void);

static void elsize_file(const char*, int option, int header, int xyzfmt);
static void eigenvalues(double, double, double, double, double, double,
                        double*, double*, double*);

static double elliptic_integral(double, double, double);
static double A_elementary_functions(double, double, double);

/* main */

int main(int argc, char** argv)
{
    int option = 0;
    int header = 0;
    int xyzfmt = 0;

    if (argc < 2) {
        usage();
        fatal("I need a valid PQR (or XYZ) file as the first "
              "command-line argument!\n\n");
    } /* argc < 2 */

    if (argc > 4) {
        usage();
        fatal("I am confused by too many arguments!\n\n");
    } /* arc > 4 */

    if (argc == 4) {
        if (strcmp(argv[3], "-xyz") == 0) {
            xyzfmt = 1;
            --argc;
        } else {
            usage();
            fatal("unexpected value for 3rd argument ('%s')\n\n", argv[3]);
        }
    }

    if (argc == 3) {
        if (strcmp(argv[2], "-det") == 0) {
            option = 0;
        } else if (strcmp(argv[2], "-ell") == 0) {
            option = 1;
        } else if (strcmp(argv[2], "-elf") == 0) {
            option = 2;
        } else if (strcmp(argv[2], "-abc") == 0) {
            option = 3;
        } else if (strcmp(argv[2], "-tab") == 0) {
            option = 4;
        } else if (strcmp(argv[2], "-hea") == 0) {
            option = 4;
            header = 1;
        } else if (strcmp(argv[2], "-deb") == 0) {
            option = 4;
            do_debug = 1;
        } else {
            usage();
            fatal("unexpected value for 2nd argument ('%s')\n\n", argv[2]);
        }
    } /* argc == 3 */

    elsize_file(argv[1], option, header, xyzfmt);

    exit(EXIT_SUCCESS);
}

/* functions */

#define LINE_LENGTH 256

static void elsize_file
    (const char* filename, int option, int header, int xyzfmt)
{
    FILE* file;

    char line[LINE_LENGTH];
    size_t natoms, i;

    double *x, *y, *z, *r2, *r3;
    double I11, I12, I13, I22, I23, I33;
    double molecule_mass, Ixx, Iyy, Izz;
    double A_det, A_ell, A_elf, a, b, c;

    file = fopen(filename, "r");
    if (file == NULL) {
        fatal("could not open '%s' for reading : %s\n",
              filename, strerror(errno));
    }

    debug("opened '%s' for reading\n", filename);

    /* pass 1 : count the atoms */

    natoms = 0;

    while (fgets(line, LINE_LENGTH, file) != NULL) {
        if (xyzfmt != 0) {
            double tmp;
            if (sscanf(line, "%lf %lf %lf\n", &tmp, &tmp, &tmp) == 3) {
                ++natoms;
            }
        } else {
            if (strlen(line) > 6
                && (strncmp(line, "ATOM", 4) == 0
                 || strncmp(line, "HETATM", 6) == 0)) {
                ++natoms;
            }
        }
    }

    if (!feof(file)) {
        fatal("reading error from '%s' : %s\n", filename, strerror(errno));
    }

    rewind(file);
    debug("%lu atoms found\n", (unsigned long) natoms);

    if (natoms == 0) {
        fatal("no atoms found in '%s'\n", filename);
    }

    /* pass 2 : load the numbers */

    x = (double*) malloc(5*sizeof(double)*natoms);
    if (x == NULL) {
        fatal("could not malloc() %lu bytes of memory\n",
              (unsigned long) 5*sizeof(double)*natoms);
    }

    y = x + natoms;
    z = y + natoms;

    r2 = z + natoms;
    r3 = r2 + natoms;

    i = 0; /* loaded atoms */

    while (fgets(line, LINE_LENGTH, file) != NULL) {
        if (xyzfmt != 0) {
            if (sscanf(line, "%lf %lf %lf\n", x + i, y + i, z + i) == 3) {
                r2[i] = DEFAULT_ATOM_SIZE*DEFAULT_ATOM_SIZE;
                r3[i] = r2[i]*DEFAULT_ATOM_SIZE;
                ++i;
            }
        } else {
            if (strlen(line) > 6
                && (strncmp(line, "ATOM", 4) == 0
                 || strncmp(line, "HETATM", 6) == 0)) {

                char tmp[LINE_LENGTH];
                double r;

                if (sscanf(line, "%s %s %s %s %s %lf %lf %lf %s %lf\n", tmp,
                    tmp, tmp, tmp, tmp, x + i, y + i, z + i, tmp, &r) == 10) {
                    r2[i] = r*r;
                    r3[i] = r2[i]*r;
                    ++i;
                }
            }
        }
    }

    fclose(file);

    debug("%lu atoms loaded\n", (unsigned long) i);
    assert(i == natoms);

    /* find the "center of mass" */

    {
        double x0, y0, z0;

        x0 = y0 = z0 = molecule_mass = 0.0;

        for (i = 0; i < natoms; ++i) {
            x0 += x[i]*r3[i];
            y0 += y[i]*r3[i];
            z0 += z[i]*r3[i];

            molecule_mass += r3[i];
        }

        assert(molecule_mass > 0.0);
        x0 /= molecule_mass; y0 /= molecule_mass; z0 /= molecule_mass;

        debug("center = (%.3f, %.3f, %.3f)\n", x0, y0, z0);
        for (i = 0; i < natoms; ++i) {
            x[i] -= x0;
            y[i] -= y0;
            z[i] -= z0;
        }
    }

    /* moments of inertia */

    I11 = I12 = I13 = I22 = I23 = I33 = 0.0;

    for (i = 0; i < natoms; ++i) {
        const double atom_mass = r3[i];

        /* half atom's moment of ineria */
        /* about diameter (will be added twice!) */
        const double atom_MI = r2[i]/5.0;

        const double x2 = x[i]*x[i] + atom_MI;
        const double y2 = y[i]*y[i] + atom_MI;
        const double z2 = z[i]*z[i] + atom_MI;

        I11 += atom_mass*(y2 + z2);
        I22 += atom_mass*(z2 + x2);
        I33 += atom_mass*(x2 + y2);

        I12 -= atom_mass*x[i]*y[i]; /* atoms' moments */
        I13 -= atom_mass*x[i]*z[i]; /* do not contribute */
        I23 -= atom_mass*y[i]*z[i]; /* because of symmetry */
    }

    debug("                     | %12.4g %12.4g %12.4g |\n", I11, I12, I13);
    debug(" tensor of inertia = | %12.4g %12.4g %12.4g |\n", I12, I22, I23);
    debug("                     | %12.4g %12.4g %12.4g |\n", I13, I23, I33);

    A_det = 0.0;
    if (option == 0 || option == 4) {
        const double det_I = I11*I22*I33 + 2*I12*I23*I13 - I11*I23*I23
                           - I22*I13*I13 - I33*I12*I12;
        if (det_I > 0.0) {
            const double sqrt_25_mass = sqrt(2.5/molecule_mass);
            const double sixthroot = pow(det_I, 1.0/6.0);
            A_det = sqrt_25_mass*sixthroot;
        } else {
            fatal("problem with the determinant of the "
                  "tensor of inertia: Det I = %g\n", det_I);
        }
    }

    if (option == 0) {
        printf("%.3f\n", A_det);
        goto done;
    }

    eigenvalues(I11, I12, I13, I22, I23, I33, &Ixx, &Iyy, &Izz);

    /* finally, let's find the axes! */
        
    a = sqrt(2.5*(-Ixx + Iyy + Izz)/molecule_mass);
    b = sqrt(2.5*( Ixx - Iyy + Izz)/molecule_mass);
    c = sqrt(2.5*( Ixx + Iyy - Izz)/molecule_mass);

    {
        double d;

        /* to make sure that a >= b >= c */

        if (b > a) {
            d = a; a = b; b = d;
        }

        if (c > b) {
            d = b; b = c; c = d;
        }

        if (b > a) {
            d = a; a = b; b = d;
        }
    }

    if (c <= 0.0) {
        fatal("problem with the semiaxes: a = %g, b = %g, c = %g\n", a, b, c);
    }

    if (option == 3) {
        printf("%.3f %9.3f %9.3f   \n", a, b, c);
        goto done;
    }

    A_ell = A_elf = 0.0;

    if (option == 1 || option == 4)  {
        if (b - c < TINY) { /* basically, b==c */
            /* it's exact in this particular case */
            A_ell = A_elementary_functions(a, b, c);
        } else {
            A_ell = 2/elliptic_integral(a, b, c);
        }
    }

    if (option == 1) {
        printf("%.3f\n", A_ell);
        goto done;
    }

    if (option == 2 || option == 4)
        A_elf = A_elementary_functions(a, b, c);

    if (option == 2) {
        printf("%.3f\n", A_elf);
        goto done;
    }

    if (header)
        puts("Filename                 "
             "a        b        c       A_ell    A_elf    A_det");

    if (do_debug) {
        debug("a = %g, b = %g, c = %g\n", a, b, c);
        debug("A_ell = %g, A_elf = %g, A_det = %g\n", A_ell, A_elf, A_det);
    } else {
        printf("%-20s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
               filename, a, b, c, A_ell, A_elf, A_det);
    }

done:
    free(x);
}

static void usage(void)
{
    fputs(
    "Usage: elsize your_structure.pqr [ arg ]\n\n"
    "The second argument is optional: \n"
    "    -det (default) gives A_det,\n"
    "    -ell gives A_ell (elliptic integral),\n"
    "    -elf gives A_elf (elementary functions approximation to A_ell,\n"
    "                      normally less than 0.1A apart),\n"
    "    -abc prints a, b, c (semiaxes of the effective ellipsoid,\n"
    "                         just out of curiousity)\n"
    "    -tab prints PQR file name and all of the above into a table"
    " without header\n",
    stderr);
    fputs(
    "    -hea prints same table as -tab but with a header\n"
    "    -deb prints same as -tab with some extra (debugging) information\n"
    "    -xyz uses a file containing only XYZ coordinates as input.\n\n",
    stderr);
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

    if (do_debug != 0) {
        va_list ap;

        va_start(ap, fmt);
        fputs("DEBUG : ", stdout);
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
        va_end(ap);
    }
}

static void eigenvalues
    (double I11, double I12, double I13, double I22, double I23, double I33,
     double* root1, double* root2, double* root3)
{
    /* To avoid artifacts like complex roots for nearly symmetrical molecules, 
       small enough non-diagonal moments should be set to exact zeros */

    const double TrI = I11 + I22 + I33; /* Trace I */
    double A, B, C, P, Q;

    assert(root1 && root2 && root3);

    debug("non-diagonal moments: I12/Tr I = %g, I13/Tr I = %g, I23/Tr I = %g\n", 
          I12/TrI, I13/TrI, I23/TrI);

    if (fabs(I12/TrI) < TINY)
        I12 = 0.0;

    if (fabs(I13/TrI) < TINY)
        I13 = 0.0;

    if (fabs(I23/TrI) < TINY)
        I23 = 0.0;
        
    A = - (I11 + I22 + I33);
    B = I11*I22 + I22*I33 + I33*I11 - I12*I12 - I13*I13 - I23*I23;
    C = - (I11*I22*I33 + 2*I12*I23*I13 - I11*I23*I23
         - I22*I13*I13 - I33*I12*I12);

    debug("cubic equation's coefficients: A = %g, B = %g, C = %g\n", A, B, C);

    /* reducing the cubic equation to an incomplete one */

    P = - (1.0/3.0)*A*A + B;
    Q = (2.0/27.0)*A*A*A - (1.0/3.0)*A*B + C;      

    debug("reduced cubic equation's coefficients: P = %g, Q = %g\n", P, Q);

    if (P >= 0.0 && fabs(Q) > 0.0 && fabs(P/Q) < TINY)
        P = 0.0;

    /* if molecule is exactly spherical, extra care should
       be taken because both P and Q == 0 */

    if (fabs(I12/TrI) < TINY && fabs(I13/TrI) < TINY && fabs(I23/TrI) < TINY) {
        P = 0.0;
        Q = 0.0;
    }

    if (P > 0.0) {
        fatal("problem with the cubic equation: at least some of the "
              "roots are not real!\n"
              "original cubic equation's coefficients:"
              " A = %g, B = %g, C = %g\n"
              "reduced cubic equation's coefficients: "
              "P = %g, Q = %g\n", A, B, C, P, Q);
    } else if (P < 0.0) {
        const double paux = 2*sqrt(-P/3);
        double cos_alpha = -4*Q/paux/paux/paux;
        double alpha;

        if (fabs(cos_alpha) > 1.0 + TINY)
            fatal("cannot compute acos(%g)\n", cos_alpha);

        if (cos_alpha > 1.0)
            cos_alpha = 1.0;

        if (cos_alpha < -1.0)
            cos_alpha = -1.0;

        alpha = acos(cos_alpha);

        *root1 = paux*cos(alpha/3) - A/3;
        *root2 = paux*cos(alpha/3 + M_PI/3) - A/3;
        *root3 = paux*cos(alpha/3 - M_PI/3) - A/3;
    } else { /* P is 0. Degeneration: a = b = c = -Q^{1/3} */
        if (Q > 0.0) {
            *root1 = *root2 = *root3 = -pow(Q, 1.0/3.0) - (1.0/3.0)*A;
        } else if (Q < 0.0) {
            *root1 = *root2 = *root3 = pow(-Q, 1.0/3.0) - (1.0/3.0)*A;
        } else {
            *root1 = *root2 = *root3 = - (1.0/3.0)*A;
        }
    }
}

/* $\int_0^\infty\frac{d\theta}{\sqrt{(a^2+\theta)(b^2+\theta)(c^2+\theta)}}$ */
/* by 1-4-1 Simpson method */

#define SIMPSON_ACCURACY 1.0E-9

static double function1(double, double, double, double);
static double function2(double, double, double, double);

static double elliptic_integral(double a, double b, double c)
{
    double d, step, lim, sum1, sum2, prev, sum_aux;
    size_t i, n, n0;

    if (a - b <= TINY && b - c <= TINY)
        return a;

    d = c; /* wild guess */

    n0 = 16;

    sum1 = sum2 = 0.0;
    prev = 1.0;

    /* calculating 1st integral from 0 to d/c */

    lim = d/c;

    for (n = n0; fabs(sum1 - prev) > SIMPSON_ACCURACY; n <<= 1) {
        prev = sum1;
        step = lim/n;
        const double f1zero = function1(a, b, c, 0.0);
        const double f1lim = function1(a, b, c, lim);
        sum1 = f1zero + f1lim;

        sum_aux = 0.0;
        for (i = 1; i < n; i += 2)
            sum_aux += function1(a, b, c, i*step);
        sum1 += 4*sum_aux;

        sum_aux = 0.0;
        for (i = 2; i < n; i += 2)
            sum_aux += function1(a, b, c, i*step);
        sum1 += 2*sum_aux;

        /* because I transformed the integral and 1/c appeared in front of it */
        sum1 *= step/3.0/c;
    }

    /* calculating 2nd integral from d to \infty. It is transformed to
       integral from 0 to 1/d, but the first part from 0 to some eps must
       be calculated analytically. */
                
    lim = 1/d;
    prev = 1.0; /* different enough from sum2 = 0 (initially) */

    for (n = n0; fabs(sum2 - prev) > SIMPSON_ACCURACY; n <<= 1) {
        prev = sum2;

        step = lim/n;
        const double f2zero = function2(a, b, c, 0.0);
        const double f2lim = function2(a, b, c, lim);
        sum2 = f2zero + f2lim;

        sum_aux = 0.0;
        for (i = 1; i < n; i += 2)
            sum_aux += function2(a, b, c, i*step);
        sum2 += 4*sum_aux;

        sum_aux = 0.0;
        for (i = 2; i < n; i += 2)
            sum_aux += function2(a, b, c, i*step);
        sum2 += 2*sum_aux;

        sum2 *= 2*step/3.0;
    }

    return sum1 + sum2;
}

static double function1(double a, double b, double c, double x)
{
    const double ac = a/c;
    const double bc = b/c;
    const double f = (x + ac*ac)*(x + bc*bc)*(x + 1.0);

    return 1.0/sqrt(f);
}

static double function2(double a, double b, double c, double x)
{
    const double xx = x*x;
    const double f = (1.0 + a*a*xx)*(1.0 + b*b*xx)*(1.0 + c*c*xx);

    return 1.0/sqrt(f);
}

static double A_elementary_functions(double a, double b, double c)
{
    if (a - b <= TINY && b - c <= TINY) /* means that a == b && b == c */
        return a;

    {
        const double gamma = sqrt(1.0 - (b + c)*(b + c)/a/a/4);
        const double A_elf = 2*a*gamma/log((1.0 + gamma)/(1.0 - gamma));
        return A_elf;
    }
}
