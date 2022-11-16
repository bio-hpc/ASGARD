/* A program to compute SAXS and ASAXS profiles from MD simulation
   More details are in JCP 2009; 130; 134114

   Written by Hung Nguyen, Case's group, 2013 */

#include "sphere_lebedev_rule.hpp"
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

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

double avogadro = 6.02214129e23;
double PI = 3.14159265359;

struct coeff_f {
    double a1, b1, a2, b2, a3, b3, a4, b4, c;
} H, C, O, N, P, S, Fe, CL, BR, NA, K, RB, CS, MG, SR;

struct coordinate {
    string type;        // Atom type
    double x, y, z;
    size_t nHyd;        // Number of hydrogen atoms attached
};

string pdb_solu, pdb_solv, outfile, exper;

double qcut = 1;        // Cutoff q
double anom_f = 0;      // Variable controls whether on-edge or off-edge calculation
double dcutoff = 5;
bool expli = 0;         // Account for explicit H atoms in pdb file
bool tight = 0;         // Use tighter convergence for Lebedev quadrature
bool corr = 0;          // Using corrected atomic factor for water as in J Chem Phys 2000, 113, 9149
double dq = 0.01;
unsigned ncpus;

// Declare functions to use
static void usage ();

static void read_pdb (const string &pdb_file,
                      vector< vector<coordinate> > &pdb_coord,
                      vector<unsigned> &weight);

static void mergeH (vector<coordinate> &model);

static vector<coordinate> solute_coord (vector<coordinate> &model);

static vector<coordinate> strip (vector<coordinate> &model,
                                 vector<coordinate> &solu_coord,
                                 double dcutoff);

static double atom_fact (const coeff_f &atom,
                         double q);

static double f_atm (const string &type,
                     double q,
                     size_t nHyd,
                     bool corr);

static complex<double> form_factor (const vector<coordinate> &model,
                                    double q,
                                    const coordinate &Leb,
                                    double anom_f,
                                    bool corr,
                                    bool expli);

static vector< complex<double> > model_form_factor (const vector< vector<coordinate> > &pdb_coord,
                                                    double q,
                                                    const coordinate &Leb,
                                                    double anom_f,
                                                    bool corr,
                                                    bool expli);

static complex<double> mean_complex (const vector< complex<double> > &v);

static double D11 (const vector< vector<coordinate> > &solu,
                   const vector<unsigned> &weight_solu,
                   const vector< vector<coordinate> > &solv,
                   const vector<unsigned> &weight_solv,
                   double q,
                   const coordinate &Leb,
                   double anom_f,
                   bool corr,
                   bool expli);

static double cal_I (const vector< vector<coordinate> > &solu,
                     const vector<unsigned> &weight_solu,
                     const vector< vector<coordinate> > &solv,
                     const vector<unsigned> &weight_solv,
                     double q,
                     int rule,
                     double anom_f,
                     bool corr,
                     bool expli);

//////////////////////////////////////////////////////////////////////////
int main (int argc,
          char *argv[]) {
    int option_char;
    do {
        static struct option long_options[] =
        {
            {"help",        no_argument,        0,    'h'},
            {"anom_f",      required_argument,  0,    'a'},
            {"system",      required_argument,  0,    'i'},
            {"solvent",     required_argument,  0,    'w'},
            {"qcut",        required_argument,  0,    'q'},
            {"dq",          required_argument,  0,    'd'},
            {"cutoff",      required_argument,  0,    'c'},
            {"expli",       no_argument,        0,    'e'},
//          {"corrected",   no_argument,        0,    'b'},
//          {"exper",       required_argument,  0,    'x'},
            {"output",      required_argument,  0,    'o'},
            {"ncpus",       required_argument,  0,    'n'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        option_char = getopt_long (argc, argv, "ha:i:w:q:d:c:eo:n:", long_options, &option_index);
        // followed by 1 colon - required an argument; 2 colon - not required argument
        if (option_char == -1)
            usage();
        else switch (option_char) {
            case 'h':
            case '?':
                usage();
                exit (0);
            case 'a':
                anom_f = atof (optarg);
                break;
            case 'i':
                pdb_solu = optarg;
                break;
            case 'w':
                pdb_solv = optarg;
                break;
            case 'q':
                qcut = atof (optarg);
                break;
            case 'd':
                dq = atof (optarg);
                break;
            case 'c':
                dcutoff = atof (optarg);
                break;
            case 'e':
                expli = 1;
                break;
/*            case 'b':
                corr = 1;
                break;
            case 'x':
                exper = optarg;
                break;*/
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
    // Cl-
    CL.a1 = 18.2915;        CL.b1 = 0.00660;
    CL.a2 = 7.20840;        CL.b2 = 1.17170;
    CL.a3 = 6.53370;        CL.b3 = 19.5424;
    CL.a4 = 2.33860;        CL.b4 = 60.4486;
    CL.c  = -16.378;
    // Br-
    BR.a1 = 17.1718;        BR.b1 = 2.20590;
    BR.a2 = 6.33380;        BR.b2 = 19.3345;
    BR.a3 = 5.57540;        BR.b3 = 0.28710;
    BR.a4 = 3.72720;        BR.b4 = 58.1535;
    BR.c  = 3.17760;
    // Na+
    NA.a1 = 3.25650;        NA.b1 = 2.66710;
    NA.a2 = 3.93620;        NA.b1 = 6.11530;
    NA.a3 = 1.39980;        NA.b3 = .200100;
    NA.a4 = 1.00320;        NA.b4 = 14.0390;
    NA.c  = .404000;
    // K+
    K.a1  = 7.95780;        K.b1  = 12.6331;
    K.a2  = 7.49170;        K.b2  = .767400;
    K.a3  = 6.35900;        K.b3  = -.00200;
    K.a4  = 1.19150;        K.b4  = 31.9128;
    K.c   = -4.9978;
    // Rb+
    RB.a1 = 17.5816;        RB.b1 = 1.71390;
    RB.a2 = 7.65980;        RB.b2 = 14.7957;
    RB.a3 = 5.89810;        RB.b3 = 0.16030;
    RB.a4 = 2.78170;        RB.b4 = 31.2087;
    RB.c  = 2.07820;
    // Cs+
    CS.a1 = 20.3524;        CS.b1 = 3.55200;
    CS.a2 = 19.1278;        CS.b2 = .308600;
    CS.a3 = 10.2821;        CS.b3 = 23.7128;
    CS.a4 = .961500;        CS.b4 = 59.4565;
    CS.c  = 3.27910;
    //Mg2+
    MG.a1 = 3.49880;        MG.b1 = 2.16760;
    MG.a2 = 3.83780;        MG.b2 = 4.75420;
    MG.a3 = 1.32840;        MG.b3 = .185000;
    MG.a4 = .849700;        MG.b4 = 10.1411;
    MG.c  = .485300;
    // Sr2+
    SR.a1 = 18.0874;        SR.b1 = 1.49070;
    SR.a2 = 8.13730;        SR.b2 = 12.6963;
    SR.a3 = 2.56540;        SR.b3 = 24.5651;
    SR.a4 = -34.193;        SR.b4 = -0.01380;
    SR.c  = 41.4025;

    ofstream OUTPUT (outfile.c_str());
    if (OUTPUT.is_open()) {
        vector< vector<coordinate> > solu_box, solv_box;
        vector<unsigned> weight_solu;
        vector<unsigned> weight_solv;
        cout << "Reading pdb ...\n";
        read_pdb (pdb_solu, solu_box, weight_solu);
        read_pdb (pdb_solv, solv_box, weight_solv);

        /* Currently assume :
           + weight_solu and weight_solv are equal, this is needed for the statistics later
           + sizes of solu and solv must be equal*/

        if (solu_box.size() != solv_box.size()) {
            cout << "!!!!!!!!!!!  Please keep the total snapshots of solute and solvent equal   !!!!!!!!!!!!!\n";
            cout << "Currently        # solute = " << solu_box.size() << ",      # solvent = " << solv_box.size() << ".   QUIT!!!!\n";
            exit (0);
        }
        for (size_t i = 0; i < weight_solu.size(); i++)
            weight_solv[i] = weight_solu[i];

        if (not expli) {
            cout << "Merging H atoms ...\n";

            #pragma omp parallel for shared (solu_box, solv_box)
            for (size_t i = 0; i < solu_box.size(); i++) {
                mergeH (solu_box[i]);
                mergeH (solv_box[i]);
            }
        }
        cout << "Stripping ...\n";
        vector< vector<coordinate> > solu;
        solu.resize(solu_box.size());
        for (size_t i = 0; i < solu.size(); i++)
            solu[i] = solute_coord (solu_box[i]);

        vector< vector<coordinate> > solu_strip, solv_strip;
        solu_strip.resize(solu_box.size());
        solv_strip.resize(solv_box.size());

        #pragma omp parallel for schedule(dynamic) shared (solu_box, solv_box, solu, solu_strip, solv_strip)
        for (size_t i = 0; i < solu_box.size(); i++) {
            solu_strip[i] = strip (solu_box[i], solu[i], dcutoff);
            solv_strip[i] = strip (solv_box[i], solu[i], dcutoff);
        }
        solu_box.clear();   solv_box.clear();

        vector<double> q;
        for (size_t i = 0; i <= floor (qcut/dq); i++)
            q.push_back(i*dq);
        OUTPUT << "# Program options:\n";
        OUTPUT << "#      + Solute          " << pdb_solu << endl;
        OUTPUT << "#      + Solvent         " << pdb_solv << endl;
        OUTPUT << "#      + Anomalous f'    " << anom_f << endl;
        OUTPUT << "#      + Explicit H      ";
        if (expli)
            OUTPUT << "ON\n";
        else OUTPUT << "OFF\n";
        OUTPUT << "#      + Tight conv      ";
        if (tight)
            OUTPUT << "ON\n";
        else OUTPUT << "OFF\n";
        OUTPUT << "\n\n####     q       Intensity\n";

        vector<double> intensity;
        intensity.resize(q.size());
        cout << "Computing SAXS ...\n";
        for (size_t i = 0; i < q.size(); i++) {
            size_t rule;
            if (not tight)
                rule = 4 + floor (q[i]/.04);
            else
                rule = 5 + floor (q[i]/.03);
            while ((available_table(rule) == 0) and (rule < 65))      // Maximum rule is 65, rarely use up to this number though
                rule++;

            intensity[i] = cal_I (solu_strip, weight_solu, solv_strip, weight_solv, q[i], rule, anom_f, corr, expli);
            // Output results
            OUTPUT << setw(12) << setprecision(8) << setiosflags(ios::fixed) << q[i];
            OUTPUT << setw(16) << setprecision(8) << setiosflags(ios::fixed) << scientific << intensity[i];
            OUTPUT << endl;
        }
        OUTPUT.close();
    } else {
        cout << "Unable to write to file " << outfile << endl;
        exit (0);
    }
    return 0;
}

//////////////////////////////    END MAIN      //////////////////////////////////

static void usage () {
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    cout << "                        A program for computing Small angle X-ray scattering curves from MD simulation\n";
    cout << "                                         Author - Hung Nguyen    tienhung@rutgers.edu\n";
    cout << "                                                   Casegroup 2013\n\n";
    cout << "Usage:  SAXS_MD  -i   --system       pdb file of the solute\n";
    cout << "                 -w   --solvent      pdb file of the solvent\n";
    cout << "                 -q   --qcut         momentum transfer q cutoff [default 1.0 A^-1]\n";
    cout << "                 -d   --dq           q spacing [default 0.01 A^-1]\n";
    cout << "                 -c   --cutoff       distance cutoff to the solute, keep only waters and ions within cutoff from the solute [default 5A]\n";
    cout << "                 -t   --tight        use tighter convergence criteria for Lebedev quadrature\n";
    cout << "                 -a   --anom_f       f' for anomalous scattering, used for ASAXS calculation,\n";
    cout << "                                     currently only support Rb+, Sr2+ and Br- [default 0: off-edge]\n";
    cout << "                 -e   --expli        flag for accounting for explicit H atoms in pdb file\n";
//  cout << "                 -b   --corrected    using corrected atomic factor for water\n";
//  cout << "                 -x   --exper        experiment data file for q generating, expect the first column is q (A^-1)\n";
    cout << "                 -o   --output       output file\n";
#ifdef _OPENMP
    cout << "                 -n   --ncpus        number of cpus used [default: 0 - using all available cpus]\n";
#endif
}

///////////////////////////////////////////
// Read pdb file, each model to one vector
///////////////////////////////////////////
static void read_pdb (const string &pdb_file,
                      vector< vector <coordinate> > &pdb_coord,
                      vector<unsigned> &weight) {
    ifstream PDBFILE (pdb_file.c_str());
    if (PDBFILE.is_open()) {
        string line;
        vector<coordinate> model;
        while (getline(PDBFILE, line))
            if (line.find("MODEL") == 0)
                model.clear();          // Remove the previous model, ready to read the new model
            else if ((line.find("ATOM") == 0) or (line.find("HETATM") == 0)) {
                coordinate coord;
                coord.x = atof (line.substr(30,8).c_str());
                coord.y = atof (line.substr(38,8).c_str());
                coord.z = atof (line.substr(46,8).c_str());

                string atm = line.substr(12,4);
                if (atm.substr(0,1) == " ")
                    atm = atm.substr(1,3);

                string type;
                string ion = atm.substr(0,2);
                if ((ion == "Na") or (ion == "NA"))
                    type = "Na+";
                else if (ion == "K ")
                    type = "K+";
                else if ((ion == "Rb") or (ion == "RB"))
                    type = "Rb+";
                else if ((ion == "Cs") or (ion == "CS"))
                    type = "Cs+";
                else if ((ion == "Cl") or (ion == "CL"))
                    type = "Cl-";
                else if ((ion == "Br") or (ion == "BR"))
                    type = "Br-";
                else if ((ion == "Mg") or (ion == "MG"))
                    type = "Mg2+";
                else if ((ion == "Sr") or (ion == "SR"))
                    type = "Sr2+";
                else {
                    type = atm.substr(0,1);
                    if ((line.find("WAT") != string::npos) or (line.find("SPC") != string::npos) or (line.find("T3P") != string::npos) or (line.find("T4E") != string::npos))
                        type += "w";        // Append "w" to differentiate water O and H atoms
                }
                coord.type = type;
                coord.nHyd = 0;
                model.push_back (coord);
               } else if ((line.find("ENDMDL") == 0) or (line.find("END") == 0))
                pdb_coord.push_back (model);
            else if (line.find("WEIGHT") == 0) {
                unsigned w = atoi (line.substr(9,6).c_str());
                // Assign weight to the current model
                weight.resize(pdb_coord.size()+1);
                weight.back() = w;
            }
        // In case there is no MODEL and ENDMDL keyword
        if (pdb_coord.size() == 0)
            pdb_coord.push_back (model);

        PDBFILE.close();
    } else {
        cerr << "Unable to open file " << pdb_file << endl;
        exit (0);
    }
    // Models must have the same size
    for (size_t i = 0; i < pdb_coord.size(); i++)
        if (pdb_coord[i].size() != pdb_coord[0].size()) {
            cerr << "Model " << i << " does not have the same size as the others. Quit!!!!\n";
            exit (0);
        }
    // Fill 1 to all the other models (models that don't have WEIGHT keyword)
    weight.resize(pdb_coord.size());
    for (size_t i = 0; i < weight.size(); i++)
        if (weight[i] <= 1)
            weight[i] = 1;
}

///////////////////////////////////////////////////
// Check if this atom belongs to solute or solvent
///////////////////////////////////////////////////
static bool check_solute (coordinate &coord) {
    if ((coord.type == "Ow") or (coord.type == "Hw") or (coord.type == "Na+") or (coord.type == "K+") or (coord.type == "Rb+") or \
        (coord.type == "Cs+") or (coord.type == "Mg2+") or (coord.type == "Sr2+") or (coord.type == "Cl-") or (coord.type == "Br-"))
        return 0;
    else return 1;
}

//////////////////////////
// Get solute coordinates
//////////////////////////
static vector<coordinate> solute_coord (vector<coordinate> &model) {
    vector<coordinate> solu;
    for (size_t i = 1; i < model.size(); i++)
        if (check_solute (model[i]))
            solu.push_back (model[i]);
    return solu;
}

//////////////
// Stripping
//////////////
static vector<coordinate> strip (vector<coordinate> &model,
                                 vector<coordinate> &solu_coord,
                                 double dcutoff) {
    vector<coordinate> solu_strip;
    for (size_t i = 0; i < model.size(); i++)
        if (not check_solute (model[i]))
            for (size_t j = 0; j < solu_coord.size(); j++) {
                double x = model[i].x - solu_coord[j].x;
                double y = model[i].y - solu_coord[j].y;
                double z = model[i].z - solu_coord[j].z;
                if (x*x + y*y + z*z <= dcutoff*dcutoff) {
                    solu_strip.push_back (model[i]);
                    break;
                }
            }
        else solu_strip.push_back (model[i]);
    return solu_strip;
}

///////////////////////////////////////////////////////
// Merge H atoms into heavier atoms by distance-based
///////////////////////////////////////////////////////
static void mergeH (vector<coordinate> &model) {
    for (size_t i = 0; i < model.size(); i++)
        if (model[i].type == "H") {
            double min_squared = 100;
            size_t index;
            for (size_t j = 0; j < model.size(); j++)
                if (model[j].type != "H") {
                    double x = model[i].x - model[j].x;
                    double y = model[i].y - model[j].y;
                    double z = model[i].z - model[j].z;
                    double distsq = x*x + y*y + z*z;
                    if (distsq < min_squared) {
                        min_squared = distsq;
                        index = j;
                    }
                }
            model[index].nHyd++;
        }
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
static double f_atm (const string &type,
                     double q,
                     size_t nHyd,
                     bool corr) {
    double atomic_factor;
    if (type == "C")
        atomic_factor = atom_fact (C, q) + nHyd*atom_fact (H, q);
    else if ((type == "O") or (type == "Ow")) {
        atomic_factor = atom_fact (O, q) + nHyd*atom_fact (H, q);
        if ((corr) and (type == "Ow") and (nHyd == 0))     // If correction and explicit
            atomic_factor *= 1.10595;
    } else if (type == "N")
        atomic_factor = atom_fact (N, q) + nHyd*atom_fact (H, q);
    else if (type == "P")
        atomic_factor = atom_fact (P, q) + nHyd*atom_fact (H, q);
    else if (type == "S")
        atomic_factor = atom_fact (S, q) + nHyd*atom_fact (H, q);
    else if (type == "Fe")
        atomic_factor = atom_fact (Fe, q);
    else if ((type == "H") or (type == "Hw")) {
        atomic_factor = atom_fact (H, q);
        if ((corr) and (type == "Hw"))
            atomic_factor *= 0.5762;
    } else if (type == "Cl-")
        atomic_factor = atom_fact (CL, q);
    else if (type == "Br-")
        atomic_factor = atom_fact (BR, q);
    else if (type == "Na+")
        atomic_factor = atom_fact (NA, q);
    else if (type == "K+")
        atomic_factor = atom_fact (K, q);
    else if (type == "Rb+")
        atomic_factor = atom_fact (RB, q);
    else if (type == "Cs+")
        atomic_factor = atom_fact (CS, q);
    else if (type == "Mg2+")
        atomic_factor = atom_fact (MG, q);
    else if (type == "Sr2+")
        atomic_factor = atom_fact (SR, q);
    else {
        cerr << "Unable to recognized atom " << type << endl;
        exit (0);
    }
    return atomic_factor;
}
////////////////////////////////////////////////////
// Calculate the form factor of the pdb file A(q)
////////////////////////////////////////////////////
static complex<double> form_factor (const vector<coordinate> &model,
                                    double q,
                                    const coordinate &Leb,
                                    double anom_f,
                                    bool corr,
                                    bool expli) {
    complex<double> f (0, 0);
    for (size_t i = 0; i < model.size(); i++)
        if ((expli) or (model[i].type != "H")) {        // Not account for H atoms for the implicit case
            double atomic_factor;
            if (expli)
                atomic_factor = f_atm (model[i].type, q, 0, corr);
            else
                atomic_factor = f_atm (model[i].type, q, model[i].nHyd, corr);
            if ((model[i].type == "Rb+") or (model[i].type == "Cs+") or (model[i].type == "Br-"))
                atomic_factor += anom_f;

            double qr = q * (Leb.x*model[i].x + Leb.y*model[i].y + Leb.z*model[i].z);
            f = f + atomic_factor * exp(complex<double> (0,1) * qr);
        }
    return f;
}
////////////////////////////////////////////////////////////////////////////////
// Compute the form factors for all models, save into a vector for later access
////////////////////////////////////////////////////////////////////////////////
static vector< complex<double> > model_form_factor (const vector< vector<coordinate> > &pdb_coord,
                                                    double q,
                                                    const coordinate &Leb,
                                                    double anom_f,
                                                    bool corr,
                                                    bool expli) {
    vector< complex<double> > mod_f;
    mod_f.resize (pdb_coord.size());
    for (size_t i = 0; i < pdb_coord.size(); i++)
        mod_f[i] = form_factor (pdb_coord[i], q, Leb, anom_f, corr, expli);
    return mod_f;
}
///////////////////////////////////////////////////////////////////////////////////////
// Compute the mean of a complex vector with weights, v and w must have the same sizes
///////////////////////////////////////////////////////////////////////////////////////
static complex<double> mean_complex (const vector< complex<double> > &v,
                                     const vector<unsigned> &w) {
    complex<double> mean (0, 0);
    unsigned total_weight = 0;
    for (size_t i = 0; i < v.size(); i++) {
        double weight = w[i];
        mean = mean + weight*v[i];
        total_weight += w[i];
    }
    double inv_size = 1.0/total_weight;
    return mean*inv_size;
}
/////////////////////////////////////////////////////////////////////////////////
// Compute D11(q), as described in Makowski et al, J Chem Phys 2009, 130, 134114
/////////////////////////////////////////////////////////////////////////////////
static double D11 (const vector< vector <coordinate> > &solu,
                   const vector<unsigned> &weight_solu,
                   const vector< vector <coordinate> > &solv,
                   const vector<unsigned> &weight_solv,
                   double q,
                   const coordinate &Leb,
                   double anom_f,
                   bool corr,
                   bool expli) {
    vector< complex<double> > solu_f = model_form_factor (solu, q, Leb, anom_f, corr, expli);
    vector< complex<double> > solv_f = model_form_factor (solv, q, Leb, anom_f, corr, expli);
    complex<double> ensbl_solu = mean_complex (solu_f, weight_solu);
    complex<double> ensbl_solv = mean_complex (solv_f, weight_solv);

    unsigned totalw_solu = 0;
    unsigned totalw_solv = 0;
    for (size_t i = 0; i < weight_solu.size(); i++)
        totalw_solu += weight_solu[i];
    for (size_t i = 0; i < weight_solv.size(); i++)
        totalw_solv += weight_solv[i];

    double diff_solu = 0;
    for (size_t i = 0; i < solu.size(); i++)
        diff_solu += weight_solu[i] * norm (solu_f[i] - ensbl_solu);
    diff_solu /= totalw_solu;

    double diff_solv = 0;
    for (size_t i = 0; i < solv.size(); i++)
        diff_solv += weight_solv[i] * norm (solv_f[i] - ensbl_solv);
    if (totalw_solv > 1)
        diff_solv *= (double) (totalw_solv + 1) / (totalw_solv*(totalw_solv - 1));

    return norm(ensbl_solu - ensbl_solv) + diff_solu - diff_solv;
}
/////////////////////////////////////////////////////////
// Integrate over the sphere using Lebedev quadrature
////////////////////////////////////////////////////////
static double cal_I (const vector< vector<coordinate> > &solu,
                     const vector<unsigned> &weight_solu,
                     const vector< vector<coordinate> > &solv,
                     const vector<unsigned> &weight_solv,
                     double q,
                     int rule,
                     double anom_f,
                     bool corr,
                     bool expli) {
    double intensity = 0;
    if (q > 0) {
        int Npoint = order_table (rule);    // Number of points in the unit sphere, for integration using Lebedev quadrature

        // Generate points on the unit sphere; x,y,z coordinates; w weight
        double *leb_weight = new double[Npoint];
        double *x_leb = new double[Npoint];
        double *y_leb = new double[Npoint];
        double *z_leb = new double[Npoint];
        ld_by_order (Npoint, x_leb, y_leb, z_leb, leb_weight);

        #pragma omp parallel for schedule (dynamic) shared (q, rule, anom_f, corr, expli) reduction (+ : intensity)
        for (size_t i = 0; i < Npoint; i++)
            // Only need to compute I for one hemisphere since I(q) = I(-q)
            if (z_leb[i] >= 0.) {
                size_t scale = 1;
                if (z_leb[i] > 0.)
                    scale = 2;
                coordinate Leb;
                Leb.x = x_leb[i];    Leb.y = y_leb[i];    Leb.z = z_leb[i];
                intensity += leb_weight[i] * scale * D11 (solu, weight_solu, solv, weight_solv, q, Leb, anom_f, corr, expli);
            }
    } else {        // q = 0 case
        coordinate Leb;
        Leb.x = 0;     Leb.y = 0;     Leb.z = 0;
        intensity = D11 (solu, weight_solu, solv, weight_solv, q, Leb, anom_f, corr, expli);
    }
    return intensity;   // Since this is averaging, not integrating, one does not multiply by 4*PI, aka the two 4PI cancel each other
}
