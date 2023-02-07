/* A program to generate SAXS and ASAXS profiles from 3D-RISM grids 
   More details are in Hung Nguyen et al, JCP 141, 2014; DOI: 10.1063/1.4896220

   Currently only support ASAXS calculation for Rb+ Sr2+ and Br-, and single salt solution (will add mixed soon?)
 
   Written by Hung Nguyen, Case's group, 2013 */

#include "sphere_lebedev_rule.hpp"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <math.h>
#include <cmath>
#include <complex>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <dirent.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

double avogadro27 = 6.02214129e-4;
double PI = 3.14159265359;

struct dx_type {
    string type;
    double conc;
    vector<double> value;
    vector<size_t> ngrid;
    vector<double> origin;
    vector<double> delta;
};

struct coeff_f {
    double a1, b1, a2, b2, a3, b3, a4, b4, c;
} H, C, O, N, P, S, Fe;

struct coordinate {
    string type;        // Atom type
    double x, y, z;
    double r;
    double B_factor;
    size_t nHyd;        // Number of hydrogen atoms attached
};

vector<dx_type> dx;
vector<coordinate> pdb_coord;

string dx_dir, pdb_file, outfile, exper;

double qcut = 0.5;      // Cutoff q
double anom_f = 0;      // f' in anomalous scattering
double df_on = 0;
double df_off = 0;
double conc_salt = 0;   // Bulk concentration of salt
double conc_wat = 55.34;// Water concentration
double cutoff = 20;
bool expli = 0;         // Account for explicit H atoms in pdb file
bool tight = 0;         // Control Lebedev quadrature tight or loose convergence
bool flex = 0;          // Account for flexibility using B-factor in the PDB file
//bool corr = 0;            // Using corrected atomic factor for water as in J Chem Phys 2000, 113, 9149
bool decomp = 0;        // Decomposing intensity into site contributions
bool off_cutoff = 0;
double dq = 0.01;
unsigned ncpus;

// Declare functions to use
static void usage ();

static void read_dx (const string &dx_dir,
                     vector<dx_type> &dx);

static bool check_dx (const vector<dx_type> &dx);

static void read_pdb (const string &pdb_file,
                      vector<coordinate> &pdb_coord);

static void mergeH (vector<coordinate> &pdb_coord);

static void extreme_pdb (const vector<coordinate> &pdb_coord,
                         coordinate &max,
                         coordinate &min);

static void read_exp (const string &exper,
                      vector<double> &q);

static coordinate get_coordinate (const dx_type &dx,
                                  size_t index);

static size_t cal_index (const dx_type &dx,
                         size_t x,
                         size_t y,
                         size_t z);

static double atom_fact (const coeff_f &atom,
                         double q);

static void assign_val (const string &type,
                        double q,
                        size_t nHyd,
                        bool expli,
                        double &atomic_factor);

static vector<size_t> list_cutoff (const dx_type &dx,
                                   const vector<coordinate> &pdb_coord,
                                   const coordinate &max,
                                   const coordinate &min,
                                   double cutoff);

static complex<double> form_factor (const vector<coordinate> &pdb_coord,
                                    const coordinate &q_vector,
                                    double q,
                                    bool flex,
                                    bool expli);

static complex<double> grid_factor (const dx_type &dx,
                                    bool excess,
                                    double atomic_factor,
                                    const vector<size_t> &list_cutoff,
                                    const coordinate &q_vector);

static dx_type ex_elec (const vector<dx_type> &dx,
                        double anom_f,
                        bool expli);

static vector<double> cal_I (const vector<dx_type> &dx,
                             bool excess,
                             const vector<coordinate> &pdb_coord,
                             const vector< vector<size_t> > &list_cutoff,
                             double q,
                             size_t rule,
                             bool flex,
                             double anom_f,
                             bool expli);

//////////////////////////////////////////////////////////////////////////
int main (int argc,
          char *argv[]) {
    int option_char;
    do {
        static struct option long_options[] =
        {
            {"anom_f",      required_argument,  0,  'a'},
            {"grid_dir",    required_argument,  0,  'g'},
            {"solute",      required_argument,  0,  's'},
            {"bfactor",     no_argument,        0,  'b'},
            {"conc_salt",   required_argument,  0,  'm'},
            {"conc_wat",    required_argument,  0,  'w'},
            {"cutoff",      required_argument,  0,  'c'},
            {"qcut",        required_argument,  0,  'q'},
            {"dq",          required_argument,  0,  'i'},
            {"expli",       no_argument,        0,  'e'},
            {"off_cutoff",  no_argument,        0,  'f'},
            {"decomp",      no_argument,        0,  'd'},
            {"tight",       no_argument,        0,  't'},
            {"output",      required_argument,  0,  'o'},
            {"ncpus",       required_argument,  0,  'n'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        option_char = getopt_long (argc, argv, "a:g:s:bm:w:c:q:i:efdto:n:", long_options, &option_index);
        // followed by 1 colon - require an argument; no colon - not argument required

        if (option_char == -1)
            usage();
        else switch (option_char) {
            case '?':
                usage();
                exit (0);
            case 'a':
                anom_f = atof (optarg);
                break;
            case 'g':
                dx_dir = optarg;
                break;
            case 's':
                pdb_file = optarg;
                break;
            case 'm':
                conc_salt = atof (optarg);
                break;
            case 'w':
                conc_wat = atof (optarg);
                break;
            case 'c':
                cutoff = atof (optarg);
                break;
            case 'q':
                qcut = atof (optarg);
                break;
            case 'i':
                dq = atof (optarg);
                break;
            case 'e':
                expli = 1;
                break;
            case 'f':
                off_cutoff = 1;
                break;
            case 't':
                tight = 1;
                break;
            case 'd':
                decomp = 1;
                break;
            case 'b':
                flex = 1;
                break;
            case 'o':
                outfile = optarg;
                break;
            case 'n':
                ncpus = atoi (optarg);
                break;
        }
    } while (option_char != -1);

#ifdef _OPENMP
    omp_set_dynamic(0);
    if (ncpus == 0)
        ncpus = omp_get_max_threads();
    omp_set_num_threads(ncpus);
#endif
    // Coefficient of atomic scattering factor taken from the "Intensity of diffracted intensities", International Tables for Crystallography, p554.
    // H
    H.a1 = 0.493002;        H.b1 = 10.5109;
    H.a2 = 0.322912;        H.b2 = 26.1257;
    H.a3 = 0.140191;        H.b3 = 3.14236;
    H.a4 = 0.040810;        H.b4 = 57.7997;
    H.c  = 0.003038;
    // C
    C.a1 = 2.31000;         C.b1 = 20.8439;
    C.a2 = 1.02000;         C.b2 = 10.2075;
    C.a3 = 1.58860;         C.b3 = 0.56870;
    C.a4 = 0.86500;         C.b4 = 51.6512;
    C.c  = 0.21560;
    // O
    O.a1 = 3.04850;         O.b1 = 13.2771;
    O.a2 = 2.28680;         O.b2 = 5.70110;
    O.a3 = 1.54630;         O.b3 = 0.32390;
    O.a4 = 0.86700;         O.b4 = 32.9089;
    O.c  = 0.25080;
    // N
    N.a1 = 12.2126;         N.b1 = 0.0057;
    N.a2 = 3.13220;         N.b2 = 9.8933;
    N.a3 = 2.01250;         N.b3 = 28.9975;
    N.a4 = 1.16630;         N.b4 = 0.5826;
    N.c  = -11.529;
    // P
    P.a1 = 6.43450;         P.b1 = 1.9067;
    P.a2 = 4.17910;         P.b2 = 27.157;
    P.a3 = 1.78000;         P.b3 = 0.5260;
    P.a4 = 1.49080;         P.b4 = 68.1645;
    P.c  = 1.11490;
    // S
    S.a1 = 6.90530;         S.b1 = 1.46790;
    S.a2 = 5.20340;         S.b2 = 22.2151;
    S.a3 = 1.43790;         S.b3 = 0.25360;
    S.a4 = 1.58630;         S.b4 = 56.1720;
    S.c  = 0.86690;
    // Fe
    Fe.a1 = 11.7695;        Fe.b1 = 4.76110;
    Fe.a2 = 7.35730;        Fe.b2 = 0.30720;
    Fe.a3 = 3.52220;        Fe.b3 = 15.3535;
    Fe.a4 = 2.30450;        Fe.b4 = 76.8805;
    Fe.c  = 1.03690;

    ofstream OUTPUT (outfile.c_str());
    if (OUTPUT.is_open()) {
        cout << "Reading input files ...\n";
        read_dx (dx_dir, dx);
        read_pdb (pdb_file, pdb_coord);
        if (not expli)
           mergeH (pdb_coord);

        vector<double> q;
        for (size_t i = 0; i <= floor (qcut/dq); i++)
            q.push_back(i*dq);

        cout << "Making lists ...\n";
        bool check = check_dx (dx);
        coordinate max_pdb, min_pdb;
        extreme_pdb (pdb_coord, max_pdb, min_pdb);

        vector< vector<size_t> > list;
        // Making cutoff list
        // Heterogeneous grid
        if (check) {
            if (not decomp) {   // Currently don't support merging heterogeneous grids into a single excess electron grid
                cout << "Using heterogeneous grids -> Auto switch to decomp = 1\n";
                decomp = 1;
            }
            list.resize(dx.size());
/*            list_exclV.resize(dx.size());
            for (size_t i = 0; i < dx.size(); i++) {
                bool flag = 0;
                if (i > 0)
                    for (int j = i-1; j >= 0; j--)
                        if ((dx[i].ngrid == dx[j].ngrid) and (dx[i].origin == dx[j].origin)) {  // These two grids are essential the same
                            flag = 1;
//                            list_exclV[i] = list_exclV[j];
                        }
//                if (not flag)
//                    list_exclV[i] = list_cutoff (dx[i], pdb_coord, max_pdb, min_pdb, -1);
            }*/

            if (off_cutoff)
                for (size_t i = 0; i < dx.size(); i++)
                    for (size_t j = 0; j < dx[i].value.size(); j++)
                        list[i].push_back(j);
            else
                for (size_t i = 0; i < dx.size(); i++) {
                    bool flag = 0;
                    if (i > 0)
                        for (int j = i-1; j >=0; j--)
                            if ((dx[i].ngrid == dx[j].ngrid) and (dx[i].origin == dx[j].origin)) {  // These two grids are essential the same
                                flag = 1;
                                list[i] = list[j];
                            }
                    if (not flag)
                        list[i] = list_cutoff (dx[i], pdb_coord, max_pdb, min_pdb, cutoff);
                }

        // Homogeneous grids
        } else {
            list.resize(1);
//            list_exclV.resize(1);
//            list_exclV[0] = list_cutoff (dx[0], pdb_coord, max_pdb, min_pdb, -1);
            if (off_cutoff)
                for (size_t i = 0; i < dx[0].value.size(); i++)
                    list[0].push_back(i);
            else
                list[0] = list_cutoff (dx[0], pdb_coord, max_pdb, min_pdb, cutoff);
        }

        cout << "List size = " << list[0].size() << endl;

        OUTPUT << "# SAXS ---- A Program for Computing Small Angle X-ray Scattering Intensity from 3D-RISM\n";
        OUTPUT << "#                     Author -- Hung Nguyen, tienhung@rutgers.edu\n";
        OUTPUT << "#                                  Casegroup 2013\n\n\n";
        OUTPUT << "#  Program options:\n\n";
        OUTPUT << "#   + Grid                   " << dx_dir << "/guv*\n";
        OUTPUT << "#   + Pdb                    " << pdb_file << endl;
        OUTPUT << "#   + B-factor               ";
        if (flex)
            OUTPUT << "ON\n";
        else
            OUTPUT << "OFF\n";
        OUTPUT << "#   + Salt  conc.            " << conc_salt << " mol/l\n";
        OUTPUT << "#   + Water conc.            " << conc_wat << " mol/l\n";
        OUTPUT << "#   + Use all grid points    ";
        if (off_cutoff)
            OUTPUT << "ON\n";
        else {
            OUTPUT << "OFF\n";
            OUTPUT << "#   + Space cutoff           " << cutoff << " Angstrom\n";
        }
        OUTPUT << "#   + Anomalous f'           " << anom_f << endl;
        OUTPUT << "#   + Explicit hydrogen      ";
        if (expli)
            OUTPUT << "ON\n";
        else
            OUTPUT << "OFF\n";
        OUTPUT << "#   + Tight convergence      ";
        if (tight)
            OUTPUT << "ON\n";
        else
            OUTPUT << "OFF\n";

        OUTPUT << "\n####     q        Solute           ";
        if (decomp)
            for (size_t i = 0; i < dx.size(); i++)
                OUTPUT << dx[i].type << "              ";
        else
            OUTPUT << "Hyd            ";
        OUTPUT << "Total\n";

        vector<dx_type> excess_e;
        if (not decomp) {
            dx_type ex_e = ex_elec (dx, anom_f, expli);
            excess_e.push_back (ex_e);
        }

        for (size_t i = 0; i < q.size(); i++) {
            size_t rule;
            if (not tight)
                rule = 4 + floor (q[i]/.04);
            else
                rule = 5 + floor (q[i]/.03);
            while ((available_table(rule) == 0) and (rule < 65))      // Maximum rule is 65, rarely use up to this number though
                rule++;
            vector<double> I_SAXS;
            if (not decomp)
                I_SAXS = cal_I (excess_e, 1, pdb_coord, list, q[i], rule, flex, anom_f, expli);
            else
                I_SAXS = cal_I (dx, 0, pdb_coord, list, q[i], rule, flex, anom_f, expli);

            // Output results
            OUTPUT << setw(12) << setprecision(8) << setiosflags(ios::fixed) << q[i];
            for (size_t j = 0; j < I_SAXS.size(); j++)
                OUTPUT << setw(16) << setprecision(8) << setiosflags(ios::fixed) << scientific << I_SAXS[j];
            OUTPUT << endl;
        }
        OUTPUT.close();
    } else {
        cout << "Unable to write to file " << outfile << endl;
        exit (0);
    }
    return 0;
}
/////////////////////////////////////    END MAIN      /////////////////////////////////////////

static void usage () {
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    cout << "                           A program for computing Small Angle X-ray Scattering intensity from 3D-RISM\n";
    cout << "                                      Author - Hung Nguyen    tienhung@rutgers.edu\n";
    cout << "                                                   Casegroup 2013\n\n";
    cout << "Usage:   SAXS    -g   --grid_dir     folder where all the rism3d output found (expecting guv.* files there)\n";
    cout << "                 -s   --solute       pdb file of the solute\n";
    cout << "                 -m   --conc_ion     bulk concentration of salt [M]\n";
    cout << "                 -w   --conc_wat     water concentration [default 55.34M]\n";
    cout << "                 -q   --qcut         momentum transfer q cutoff [default 0.5 A^-1]\n";
    cout << "                 -i   --dq           q spacing [default 0.01 A^-1]\n";
    cout << "                 -c   --cutoff       distance cutoff [default 20 A]\n";
    cout << "                 -a   --anom_f       f' of atomic scattering factor, used for ASAXS calculation,\n";
    cout << "                                     currently only applied to Rb+, Sr2+ and Br- [default 0: off-edge]\n";
    cout << "                 -e   --expli        flag for accounting for explicit H atoms in pdb file\n";
    cout << "                 -d   --decomp       flag for decomposing SAXS intensity into site contributions (lead to 2-5x computational time)\n";
    cout << "                 -t   --tight        flag for using tighter convergence criteria for Lebedev quadrature (expect more time)\n";
    cout << "                 -f   --off_cutoff   flag for turning off cutoff, using all grid points for the calculation\n";
    cout << "                 -b   --bfactor      using B-factor in the PDB file to account for solute flexibility\n";
    cout << "                 -o   --output       output file\n";
#ifdef _OPENMP
    cout << "                 -n   --ncpus        number of cpus used [default: 0 - using all available cpus]\n";
#endif
}

/////////////////////////
// Read all the dx files
/////////////////////////
static void read_dx (const string &dx_dir,
                     vector<dx_type> &dx) {
    DIR *dir = NULL;
    struct dirent *file = NULL;
    if ((dir = opendir (dx_dir.c_str())) != NULL) {
        size_t type = 0;            // monovalent ion -> 1     ;       divalent ion -> 2
        while ((file = readdir (dir)) != NULL) {
            string filename = file->d_name;
            if (filename.substr(0,3) == "guv") {
                string path = dx_dir + "/" + filename;
                dx_type newdx;

                size_t pos = filename.find ("O");
                if (pos != std::string::npos) {
                    newdx.conc = conc_wat;
                    newdx.type = "Ow";
                } else {
                    pos = filename.find ("H1");
                    if (pos != std::string::npos) {
                        newdx.conc = 2*conc_wat;
                        newdx.type = "Hw";
                    } else {
                        pos = filename.find ("F-");
                        if (pos != std::string::npos)
                            newdx.type = "F-";
                        else {
                            pos = filename.find ("Cl-");
                            if (pos != std::string::npos)
                                newdx.type = "Cl-";
                            else {
                                pos = filename.find ("Br-");
                                if (pos != std::string::npos)
                                    newdx.type == "Br-";
                                else {
                                    pos = filename.find ("I-");
                                    if (pos != std::string::npos)
                                        newdx.type == "I-";
                                    else {
                                        pos = filename.find ("Na+");
                                        if (pos != std::string::npos) {
                                            newdx.type = "Na+";
                                            type = 1;
                                        } else {
                                            pos = filename.find ("K+");
                                            if (pos != std::string::npos) {
                                                newdx.type = "K+";
                                                type = 1;
                                            } else {
                                                pos = filename.find ("Rb+");
                                                if (pos != std::string::npos) {
                                                    newdx.type = "Rb+";
                                                    type = 1;
                                                } else {
                                                    pos = filename.find ("Cs+");
                                                    if (pos != std::string::npos) {
                                                        newdx.type = "Cs+";
                                                        type = 1;
                                                    } else {
                                                        pos = filename.find ("Li+");
                                                        if (pos != std::string::npos) {
                                                            newdx.type = "Li+";
                                                            type = 1;
                                                        } else {
                                                            pos = filename.find ("Mg2+");
                                                            if (pos != std::string::npos) {
                                                                newdx.type = "Mg2+";
                                                                type = 2;
                                                            } else {
                                                                pos = filename.find ("Sr2+");
                                                                if (pos != std::string::npos) {
                                                                    newdx.type = "Sr2+";
                                                                    type = 2;
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Read dx file
                ifstream DXFILE (path.c_str());
                if (DXFILE.is_open()) {
                    string line;
                    while (getline(DXFILE, line)) {
                        istringstream iss(line);
                        double tmp1;
                        while (iss >> tmp1)
                            newdx.value.push_back(tmp1);
                        if (line.find("object 1 class gridpositions counts") != string::npos) {
                            line.replace(0,35," ");
                            istringstream iss(line);
                            newdx.ngrid.resize(3);
                            iss >> newdx.ngrid[0] >> newdx.ngrid[1] >> newdx.ngrid[2];
                        } else if (line.find("origin") != string::npos) {
                            line.replace(0,6," ");
                            istringstream iss(line);
                            newdx.origin.resize(3);
                            iss >> newdx.origin[0] >> newdx.origin[1] >> newdx.origin[2];
                        } else if (line.find("delta") != string::npos) {
                            line.replace(0,5," ");
                            istringstream iss(line);
                            newdx.delta.resize(3);
                            double tmp2, tmp3;
                            iss >> tmp1 >> tmp2 >> tmp3;
                            if (tmp1 != 0)
                                newdx.delta[0] = tmp1;
                            else if (tmp2 != 0)
                                newdx.delta[1] = tmp2;
                            else if (tmp3 != 0)
                                newdx.delta[2] = tmp3;
                        }
                    }
                    DXFILE.close();
                    if (newdx.ngrid[0] * newdx.ngrid[1] * newdx.ngrid[2] != newdx.value.size()) {
                        cout << "Number of grid points not matched in " << file->d_name << endl;
                        exit (0);
                    } else dx.push_back(newdx);
                } else {
                    cerr << "Unable to open file " << file->d_name << endl;
                    exit (0);
                }
            }
        }
        closedir (dir);

        // Swap Ow and Hw to the beginning, for output reading convenience
        for (size_t i = 0; i < dx.size(); i++) {
            if ((dx[i].type == "Ow") and (i > 0)) {
                dx_type tmp = dx[0];
                dx[0] = dx[i];
                dx[i] = tmp;
            }
            if ((dx[i].type == "Hw") and (i != 1)) {
                dx_type tmp = dx[1];
                dx[1] = dx[i];
                dx[i] = tmp;
            }
        }

        // Set concentration depending on type
        for (size_t i = 2; i < dx.size(); i++)
            if (type == 1)
                dx[i].conc = conc_salt;
            else if (type == 2)
                if ((dx[i].type == "Sr2+") or (dx[i].type == "Mg2+"))
                    dx[i].conc = conc_salt;
                else if ((dx[i].type == "Cl-") or (dx[i].type == "F-") or (dx[i].type == "Br-") or (dx[i].type == "I-"))
                    dx[i].conc = 2*conc_salt;
    } else {
        cout << "Not able to open " << dx_dir << endl;
        exit (0);
    }
}

/////////////////////////////////////////////
// Check to see dx grids are different or not
/////////////////////////////////////////////
static bool check_dx (const vector<dx_type> &dx) {
    bool check = 0;
    for (size_t i = 1; i < dx.size(); i++) {
        if (dx[i].ngrid != dx[0].ngrid) {
            check = 1;
            break;
        }
        if (dx[i].origin != dx[0].origin) {
            check = 1;
            break;
        }
        if (dx[i].delta != dx[0].delta) {
            check = 1;
            break;
        }
    }
    return check;
}

/////////////////
// Read pdb file
//////////////////
static void read_pdb (const string &pdb_file,
                      vector<coordinate> &pdb_coord) {
    ifstream PDBFILE (pdb_file.c_str());
    if (PDBFILE.is_open()) {
        string line;
        while (getline(PDBFILE, line))
            if ((line.find("ATOM") == 0) or (line.find("HETATM") == 0)) {
                coordinate coord;
                coord.x = atof (line.substr(30,8).c_str());
                coord.y = atof (line.substr(38,8).c_str());
                coord.z = atof (line.substr(46,8).c_str());
                coord.B_factor = atof (line.substr(60,6).c_str());

                string atm = line.substr(12,4);
                if (atm.substr(0,1) == " ")
                    atm = atm.substr(1,3);
                string type = atm.substr(0,1);
                if ((atm.substr(0,2) == "FE") or (atm.substr(0,2) == "Fe"))
                    type = "Fe";

                if (type == "H")
                    coord.r = 1.2;
                else if (type == "C")
                    coord.r = 1.7;
                else if (type == "N")
                    coord.r = 1.55;
                else if (type == "O")
                    coord.r = 1.52;
                else if ((type == "S") or (type == "P"))
                    coord.r = 1.8;

                coord.type = type;
                coord.nHyd = 0;
                pdb_coord.push_back (coord);
            }
        PDBFILE.close();
    } else {
        cerr << "Unable to open file " << pdb_file << endl;
        exit (0);
    }
}

//////////////////////////////////////////////////////
// Merge H atoms into heavier atoms by distance-based
//////////////////////////////////////////////////////
static void mergeH (vector<coordinate> &pdb_coord) {
    #pragma omp parallel for schedule (dynamic)
    for (size_t i = 0; i < pdb_coord.size(); i++)
        if (pdb_coord[i].type == "H") {
            double min_squared = 100;
            size_t index;
            for (size_t j = 0; j < pdb_coord.size(); j++)
                if (pdb_coord[j].type != "H") {
                    double x = pdb_coord[i].x - pdb_coord[j].x;
                    double y = pdb_coord[i].y - pdb_coord[j].y;
                    double z = pdb_coord[i].z - pdb_coord[j].z;
                    double distsq = x*x + y*y + z*z;
                    if (distsq < min_squared) {
                        min_squared = distsq;
                        index = j;
                    }
                }
            #pragma omp critical
                pdb_coord[index].nHyd++;
        }
}

///////////////////////////////////////
// Find min and max of pdb coordinates
///////////////////////////////////////
static void extreme_pdb (const vector<coordinate> &pdb_coord,
                         coordinate &max,
                         coordinate &min) {
    double min_x = -9e9;    double max_x = 9e9;
    double min_y = -9e9;    double max_y = 9e9;
    double min_z = -9e9;    double max_z = 9e9;
    for (size_t i = 0; i < pdb_coord.size(); i++) {
        if (pdb_coord[i].x < min_x)
            min_x = pdb_coord[i].x;
        else if (pdb_coord[i].x > max_x)
            max_x = pdb_coord[i].x;
        if (pdb_coord[i].y < min_y)
            min_y = pdb_coord[i].y;
        else if (pdb_coord[i].y > max_y)
            max_y = pdb_coord[i].y;
        if (pdb_coord[i].z < min_z)
            min_z = pdb_coord[i].z;
        else if (pdb_coord[i].z > max_z)
            max_z = pdb_coord[i].z;
    }
    max.x = max_x;  max.y = max_y;  max.z = max_z;
    min.x = min_x;  min.y = min_y;  min.z = min_z;
}

/*////////////////////////////////////////////////////////////
// Read RDF of electron around water (for smearing purpose)
///////////////////////////////////////////////////////////
static void read_RDF (const string &filename,
                      double &dr,
                      vector<double> &rdf) {
    ifstream OPENFILE (filename.c_str());
    if (OPENFILE.is_open()) {
        string line;
        double diff_bk = 0;
        double tmp1 = 0; 
        while (getline(OPENFILE, line)) {
            istringstream iss(line);
            // Skip line starting with #
            size_t found = line.find_first_not_of(" \t");
            if (found != std::string::npos)
                if (line[found] == '#')
                    continue;
            double tmp2, value;
            while (iss >> tmp2 >> value)
                rdf.push_back(value);
            double diff = tmp2 - tmp1;
            if ((abs(diff_bk) < 1e-10) and (diff - diff_bk < 1e-10)) {
                cout << "dr is not a constant!! Quit!!!!\n";
                cout << diff_bk << endl << diff << endl;
                exit(0);
            }
            diff_bk = diff;
            tmp1 = tmp2;
        }
        OPENFILE.close();
        dr = diff_bk;
    }
}*/

/////////////////////////////////////////////////////
// Get grid point coordinate from the 1d index
/////////////////////////////////////////////////////
static coordinate get_coord (const dx_type &dx,
                             size_t index) {
    size_t x_index = floor ((double)index/((double)dx.ngrid[1]*dx.ngrid[2]));
    size_t y_index = floor ((double)(index - x_index * dx.ngrid[1] * dx.ngrid[2]) / (double)dx.ngrid[2]);
    size_t z_index = index - x_index*dx.ngrid[1]*dx.ngrid[2] - y_index*dx.ngrid[2];

    coordinate coord;
    coord.x = dx.origin[0] + dx.delta[0] * x_index;
    coord.y = dx.origin[1] + dx.delta[1] * y_index;
    coord.z = dx.origin[2] + dx.delta[2] * z_index;

    return coord;
}

//////////////////////////////////////////////////
// Return the 1D index of the 3D grid coordinate
/////////////////////////////////////////////////
static size_t cal_index (const dx_type &dx,
                         size_t x,
                         size_t y,
                         size_t z) {
    return x*dx.ngrid[1]*dx.ngrid[2] + y*dx.ngrid[2] + z;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Empirically compute the atomic factor of atom in vacuo,
// from "Intensity of diffracted intensities", International Tables for Crystallography, p554
/////////////////////////////////////////////////////////////////////////////////////////////
static double atom_fact (const coeff_f &atom,
                         double q) {
    double q2 = q*q/(16*PI*PI);
    return atom.a1*exp(-atom.b1*q2) + atom.a2*exp(-atom.b2*q2) + atom.a3*exp(-atom.b3*q2) + atom.a4*exp(-atom.b4*q2) + atom.c;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Assign the atomic factor for atoms, summing up the atomic factors of hydrogen atoms into heavy atom
// Using the analytical approximation for scattering factors,
// taken from "Intensity of diffracted intensities", International Tables for Crystallography, p554
/////////////////////////////////////////////////////////////////////////////////////////////////////////
static void assign_val (const string &type,
                        double q,
                        size_t nHyd,
                        bool expli,
                        double &atomic_factor) {
    if (type == "C")
        atomic_factor = atom_fact (C, q) + nHyd*atom_fact (H, q);
    else if (type == "O")
        atomic_factor = atom_fact (O, q) + nHyd*atom_fact (H, q);
    else if (type == "N")
        atomic_factor = atom_fact (N, q) + nHyd*atom_fact (H, q);
    else if (type == "P")
        atomic_factor = atom_fact (P, q) + nHyd*atom_fact (H, q);
    else if (type == "S")
        atomic_factor = atom_fact (S, q) + nHyd*atom_fact (H, q);
    else if (type == "Fe")
        atomic_factor = atom_fact (Fe, q);
    else if (type == "H")
        atomic_factor = atom_fact (H, q);

    // For grid points, use 3D Fourier transform of electron density, thus only need the total number of electron
    else if (type == "Ow")
        if (expli) {
            atomic_factor = 8;
//            if (corr) atomic_factor = 8.8476;
        } else atomic_factor = 10;
    else if (type == "Hw") {
        atomic_factor = 1;
//        if (corr) atomic_factor = 0.5762;
    } else if (type == "F-")
        atomic_factor = 10;
    else if (type == "Cl-")
        atomic_factor = 18;
    else if (type == "Br-")
        atomic_factor = 36;
    else if (type == "I-")
        atomic_factor = 54;
    else if (type == "Li+")
        atomic_factor = 2;
    else if (type == "Na+")
        atomic_factor = 10;
    else if (type == "K+")
        atomic_factor = 18;
    else if (type == "Rb+")
        atomic_factor = 36;
    else if (type == "Cs+")
        atomic_factor = 54;
    else if (type == "Mg2+")
        atomic_factor = 10;
    else if (type == "Sr2+")
        atomic_factor = 36;
    else {
        cerr << "Unable to recognize atom " << type << endl;
        exit (0);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////
// Calculate the list of indexes of grid points within cutoff distance from the solute
/////////////////////////////////////////////////////////////////////////////////////////
static vector<size_t> list_cutoff (const dx_type &dx,
                                   const vector<coordinate> &pdb_coord,
                                   const coordinate &max,
                                   const coordinate &min,
                                   double cutoff) {
    // Form an initial cube
    vector<size_t> prelist;
    #pragma omp parallel for schedule (dynamic) shared (cutoff, prelist)
    for (size_t i = 0; i < dx.value.size(); i++) {
        coordinate grid = get_coord (dx, i);
        double cut = std::max (double(5.), cutoff);       // Keep small distance, at least 5A from the solute to the box edge
        if ((grid.x > min.x - cut) and (grid.x < max.x + cut) and (grid.y > min.y - cut) and (grid.y < max.y + cut) and \
                                                                 (grid.z > min.z - cut) and (grid.z < max.z + cut))
            #pragma omp critical
            prelist.push_back (i);
    }

    // Compute distance of every grid points in the initial cube
    vector<size_t> list;
    #pragma omp parallel for shared (prelist, cutoff, list)
    for (size_t i = 0; i < prelist.size(); i++)
        for (size_t j = 0; j < pdb_coord.size(); j++) {
            coordinate grid = get_coord (dx, prelist[i]);
            double deltax = grid.x - pdb_coord[j].x;
            double deltay = grid.y - pdb_coord[j].y;
            double deltaz = grid.z - pdb_coord[j].z;
            double r_cut;

            if (cutoff >= 0)
                r_cut = cutoff;
            else    // Excluded volume
                r_cut = pdb_coord[j].r + 1.4;

            if ((deltax*deltax + deltay*deltay + deltaz*deltaz <= r_cut*r_cut)) {
                #pragma omp critical
                    list.push_back (prelist[i]);
                break;
            }
        }
    sort (list.begin(), list.end());
    return list;
}

/////////////////////////////////////////////////////////////
// Calculate the atomic form factor of the solute (pdb file)
/////////////////////////////////////////////////////////////
static complex<double> form_factor (const vector<coordinate> &pdb_coord,
                                    const coordinate &q_vector,
                                    double q,
                                    bool flex,
                                    bool expli) {
    complex<double> f (0, 0);
    for (size_t i = 0; i < pdb_coord.size(); i++)
        if ((expli) or (pdb_coord[i].type != "H")) {        // Not account for H atoms for the implicit case
            double atomic_factor;
            if (expli)
                assign_val (pdb_coord[i].type, q, 0, expli, atomic_factor);
            else
                assign_val (pdb_coord[i].type, q, pdb_coord[i].nHyd, expli, atomic_factor);

            if (flex)
                atomic_factor *= exp (-pdb_coord[i].B_factor*q*q/(16*PI*PI));

            double qr = q_vector.x*pdb_coord[i].x + q_vector.y*pdb_coord[i].y + q_vector.z*pdb_coord[i].z;
            f = f + atomic_factor * exp(complex<double> (0,1) * qr);
        }
    return f;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calculate contribution from one type of grid (water, counterions or electron)
// The boolean excess variable is used to specify whether the values are excess relatively to the bulk or not
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
static complex<double> grid_factor (const dx_type &dx,
                                    bool excess,
                                    double atomic_factor,
                                    const vector<size_t> &list_cutoff,
                                    const coordinate &q_vector) {
    // 3D Fourier transform of the electron density in the unit volume
    // f = 8*rho*sin(q_x*a/2)*sin(q_y*a/2)*sin(q_z*a/2)/(q_x*q_y*q_z)
    double sinc_x, sinc_y, sinc_z;
    if (q_vector.x != 0)
        sinc_x = 2*sin(.5*q_vector.x*dx.delta[0]) / q_vector.x;
    else
        sinc_x = dx.delta[0];
    if (q_vector.y != 0)
        sinc_y = 2*sin(.5*q_vector.y*dx.delta[1]) / q_vector.y;
    else
        sinc_y = dx.delta[1];
    if (q_vector.z != 0)
        sinc_z = 2*sin(.5*q_vector.z*dx.delta[2]) / q_vector.z;
    else
        sinc_z = dx.delta[2];

    complex<double> sum (0, 0);

    for (size_t i = 0; i < list_cutoff.size(); i++) {
        coordinate grid = get_coord (dx, list_cutoff[i]);
        double value;

        if (excess)
            value = dx.value[list_cutoff[i]];
        else
            value = dx.value[list_cutoff[i]] - 1;

        double qr = q_vector.x*grid.x + q_vector.y*grid.y + q_vector.z*grid.z;
        sum = sum + value * exp(complex<double> (0, 1) * qr);
    }
    return atomic_factor*dx.conc*avogadro27*sinc_x*sinc_y*sinc_z*sum;
}

//////////////////////////////////////////////
// How much to put into the center grid point
//////////////////////////////////////////////
/*static double increasing_center (const dx_type &dx,
                                 double cutoff,
                                 size_t Z,
                                 double dr,
                                 const vector<double> &rdf) {
    int stepx = floor (cutoff / dx.delta[0]);
    int stepy = floor (cutoff / dx.delta[1]);
    int stepz = floor (cutoff / dx.delta[2]);
    double sum = 0;
    for (int x = -stepx; x <= stepx; x++) {
        for (int y = -stepy; y <= stepy; y++)
            for (int z = -stepz; z <= stepz; z++)
                if ((x != 0) or (y != 0) or (z != 0)) {
                    double dist = distance (dx, x, y, z, 0, 0, 0);
                    size_t rdf_index = floor(dist/dr);
                    // If the rdf_index falls outside the range then assuming there is no electron outside, which is reasonable for a 4A maximum in rdf file
                    if (rdf_index <= rdf.size() - 1)
                        sum += rdf[rdf_index] * dx.delta[0] * dx.delta[1] * dx.delta[2];
                }
    }
    // sum contains all contribution from grid points within the cutoff EXCEPT the center grid
    return Z - sum;
}*/

////////////////////////////////////////////////////////////////////
// Map an excess electron map based on the atomic distribution grid
////////////////////////////////////////////////////////////////////
static dx_type ex_elec (const vector<dx_type> &dx,
                        double anom_f,
                        bool expli) {
    dx_type result;
    result.ngrid = dx[0].ngrid;
    result.delta = dx[0].delta;
    result.origin = dx[0].origin;
    result.conc = 1;    // Concentration for each sites will be precalculated and stored
    result.type = "ex_e";
    result.value.resize(dx[0].value.size(), 0);

    for (size_t type = 0; type < dx.size(); type++)
        if ((expli) or (dx[type].type != "Hw")) {
            double Z;
            if ((dx[type].type != "Ow") or (expli))
                assign_val (dx[type].type, 0, 0, expli, Z);
            else
                assign_val (dx[type].type, 0, 2, expli, Z);
            if ((dx[type].type == "Rb+") or (dx[type].type == "Sr2+") or (dx[type].type == "Br-")) 
                Z += anom_f;

            #pragma omp parallel for shared (result, type, Z)
            for (size_t x = 0; x < dx[type].ngrid[0]; x++)
                for (size_t y = 0; y < dx[type].ngrid[1]; y++)
                    for (size_t z = 0; z < dx[type].ngrid[2]; z++) {
                        size_t index_center = cal_index(dx[type], x, y, z);
                        double grid = (dx[type].value[index_center] - 1) * dx[type].conc * Z;
                        result.value[index_center] += grid;
                    }
        }
    return result;
}

/////////////////////////////////////////////////////////
// Integrate over the sphere using Lebedev quadrature
////////////////////////////////////////////////////////
static vector<double> cal_I (const vector<dx_type> &dx,
                             bool excess,
                             const vector<coordinate> &pdb_coord,
                             const vector< vector<size_t> > &list_cutoff,
                             double q,
                             size_t rule,
                             bool flex,
                             double anom_f,
                             bool expli) {
    vector<double> intensity;
    intensity.resize(dx.size()+2, 0);

    if (q > 0) {
        size_t Npoint = order_table (rule);     // Number of points in the unit sphere, for integration using Lebedev quadrature

        // Generate points on the unit sphere; x,y,z coordinates; w weight
        double *w_leb = new double[Npoint];
        double *x_leb = new double[Npoint];
        double *y_leb = new double[Npoint];
        double *z_leb = new double[Npoint];
        ld_by_order (Npoint, x_leb, y_leb, z_leb, w_leb);

        #pragma omp parallel for schedule (dynamic) shared (excess, q, x_leb, y_leb, z_leb, w_leb, anom_f, expli, flex, intensity)
        for (size_t i = 0; i < Npoint; i++)
            // Only need to compute I for one hemisphere since I(q) = I(-q)
            if (z_leb[i] >= 0.) {
                unsigned weight = 1;
                if (z_leb[i] > 0.)
                    weight = 2;

                // Scaling q_vector
                coordinate q_vector;
                q_vector.x = q*x_leb[i];    q_vector.y = q*y_leb[i];    q_vector.z = q*z_leb[i];

                vector< complex<double> > ampl;
                ampl.resize(dx.size() + 2, 0);
                // Amplitude for the solute in vacuum, index 0
                ampl[0] = form_factor (pdb_coord, q_vector, q, flex, expli);
                ampl[dx.size()+1] = ampl[0];

                // Amplitude for each of dx file, save to next indices
                for (size_t j = 0; j < dx.size(); j++) {
                    double atomic_factor;
                    if (excess)     // Excess electron map
                        atomic_factor = 1;
                    else if ((expli) or (dx[j].type != "Hw")) {    // Not account for Hw grid in the implicit case
                        if ((expli) or (dx[j].type != "Ow"))
                            assign_val (dx[j].type, q, 0, expli, atomic_factor);
                        else
                            assign_val (dx[j].type, q, 2, expli, atomic_factor);

                        if ((dx[j].type == "Rb+") or (dx[j].type == "Sr2+") or (dx[j].type == "Br-"))
                            atomic_factor += anom_f;
                    }
                    if ((expli) or (dx[j].type != "Hw")) {   // Not spend time for Hw grid in the implicit case
                        if (list_cutoff.size() == 1)
                            ampl[j+1] = grid_factor (dx[j], excess, atomic_factor, list_cutoff[0], q_vector);
                        else
                            ampl[j+1] = grid_factor (dx[j], excess, atomic_factor, list_cutoff[j], q_vector);
                        ampl[dx.size()+1] = ampl[dx.size()+1] + ampl[j+1];
                    }
                }
                // Compute intensity (complex)
                for (size_t j = 0; j < ampl.size(); j++)
                    #pragma omp critical
                        intensity[j] += w_leb[i] * weight * norm(ampl[j]);
           }
    } else {        // q = 0 case
        coordinate q_vector;
        q_vector.x = 0;     q_vector.y = 0;     q_vector.z = 0;

        vector< complex<double> > ampl;
        ampl.resize(dx.size() + 2, 0);
        // Amplitude for the solute in vacuum, index 0
        ampl[0] = form_factor (pdb_coord, q_vector, q, flex, expli);
        ampl.back() = ampl[0];

        // Amplitude for each of dx grid
        for (size_t j = 0; j < dx.size(); j++) {
            double atomic_factor;
            if (excess)     // Excess electron map
                atomic_factor = 1;
            else if ((expli) or (dx[j].type != "Hw")) {
                if ((expli) or (dx[j].type != "Ow"))
                    assign_val (dx[j].type, 0, 0, expli, atomic_factor);
                else
                    assign_val (dx[j].type, 0, 2, expli, atomic_factor);

                if ((dx[j].type == "Rb+") or (dx[j].type == "Sr2+"))
                    atomic_factor += anom_f;
            }
            if ((expli) or (dx[j].type != "Hw")) {
                if (list_cutoff.size() == 1)
                    ampl[j+1] = grid_factor (dx[j], excess, atomic_factor, list_cutoff[0], q_vector);
                else
                    ampl[j+1] = grid_factor (dx[j], excess, atomic_factor, list_cutoff[j], q_vector);
                ampl[dx.size()+1] = ampl[dx.size()+1] + ampl[j+1];
            }
        }
        // Compute intensity
        for (size_t j = 0; j < ampl.size(); j++)
            intensity[j] = norm(ampl[j]);
    }
    return intensity;   // Since this is averaging, not integrating, one does not multiply by 4*PI, aka the two 4PI cancel each other
}
