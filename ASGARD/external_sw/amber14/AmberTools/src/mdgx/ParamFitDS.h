#ifndef ParamFitDataStructures
#define ParamFitDataStructures

#ifndef PREP_API
#include "Constants.h"

#include "ChargeFitDS.h"
#include "CrdManipDS.h"
#endif

struct EquivAtomGroup {
  int natom;
  char* types;
};
typedef struct EquivAtomGroup eagrp;

struct ExtendedAtomDef {
  int inreport;   // Flag to indicate that this atom type is needed in the
                  //   final report
  int dup;        // Flag to indicate that this atom type is a duplicate,
                  //   branched and renamed from one of the original types
  double mass;    // The atom mass
  double apol;    // The atomic polarizability
  double ljsig;   // The Lennard-Jones Sigma parameter
  double ljeps;   // The Lennard-Jones Sigma parameter
  char atype[8];  // The atom type
  char* comment;  // Comment from the parameter file
};
typedef struct ExtendedAtomDef xatomdef;

struct ExtendedBondDef {
  int fitcol;     // Column of the fitting matrix to which bonds of this sort
                  //   are mapped
  int dup;        // Flag to indicate that this bond type is a duplicate,
                  //   branched and renamed from one of the original types
  int samprow;    // Instances of the bond populate this row of the sampling
                  //   table
  int ninst;      // The number of instances of this bond in the fit
  int inreport;   // Flag to indicate that this bond type is needed in the
                  //   final report
  double K;       // Spring constant
  double l0;      // Equilibrium length
  char atype[8];  // The type of atom A
  char btype[8];  // The type of atom B
  char* comment;  // Comment from the parameter file
};
typedef struct ExtendedBondDef xbonddef;

struct ExtendedHBondDef {
  int fitcol;     // Column of the fitting matrix to which this maps
  int dup;        // Flag to indicate that this H-bond type is a duplicate,
                  //   branched and renamed from one of the original types
  int ninst;      // Instances of this H-bond
  int inreport;   // Flag to indicate that this is needed in the final report
  double Aterm;   // r^12, r^10 related constants 
  double Bterm;   //
  char atype[8];  // The type of atom A
  char btype[8];  // The type of atom B
  char* comment;  // Comment from the parameter file
};
typedef struct ExtendedHBondDef xhb1012def;

struct ExtendedAngleDef {
  int fitcol;     // Column of the fitting matrix to which angles of this sort
                  //   are mapped
  int dup;        // Flag to indicate that this angle type is a duplicate,
                  //   branched and renamed from one of the original types
  int samprow;    // Instances of the angle populate this row of the sampling
                  //   table
  int ninst;      // The number of instances of this angle in the fit
  int inreport;   // Flag to indicate that this bond type is needed in the
                  //   final report
  double K;       // Spring constant
  double th0;     // Equilibrium length
  char atype[8];  // The type of atom A
  char btype[8];  // The type of atom B
  char ctype[8];  // The type of atom C
  char* comment;  // Comment from the parameter file
};
typedef struct ExtendedAngleDef xangldef;

struct TorsionTerm {
  int fitcol;     // Column of the fitting matrix to which torsions of this
                  //   sort are mapped
  int dup;        // Flag to indicate that this torsion type is a duplicate,
                  //   branched and renamed from one of the original types
  int samprow;    // Instances of the torsion populate this row of the sampling
                  //   table
  int ninst;      // The number of instances of this dihedral in the fit
  int inreport;   // Flag to indicate that this bond type is needed in the
                  //   final report
  int impr;       // Flag to indicate that this dihedral is an improper
  double K;       // The basic torsional potential for an angle phi is defined:
  double phase;   //   U = (pk/idivf) * (1 + cos(pn*phi - phase)), K = pk/idivf
  double pn;      //
  double singlet; // Indicates whether the dihedral is singlet or not; +1.0 if
                  //   there is only one term, -1.0 if there are more terms 
  char atype[8];  // 
  char btype[8];  // The A, B, C, and D atom types (i.e. CT, HC, N3) in a
  char ctype[8];  //   torsion arrangement  A--B--C--D or improper A(D)--B--C
  char dtype[8];  // 
  char* comment;  // Comment from the parameter file
};
typedef struct TorsionTerm torterm;

struct BondIndex {
  int a;          // The A and B atoms of the bond, indexed according to the
  int b;          //   system's own topology
  int key;        // The index into the master list of bonds, stored in the
                  //   parameter set and spanning all systems
};
typedef struct BondIndex bidx;

struct AnglIndex {
  int a;          // The A, B, and C atoms of the angle, indexed according to
  int b;          //   the system's own topology
  int c;          //
  int key;        // The index into the master list of angles, stored in the
                  //   parameter set and spanning all systems
};
typedef struct AnglIndex aidx;

struct TorsionIndex {
  int a;          //
  int b;          // The A, B, C, and D atoms of the torsion term, indexed
  int c;          //   according to the system's own topology
  int d;          //
  int key;        // The index into the master list of angles, stored in the
                  //   parameter set and spanning all systems
};
typedef struct TorsionIndex hidx;

struct BondMap {
  int nbond;        // The number of bonds in this system
  bidx* id;         // Numbers of atoms participating in each bond, and indices
                    //   into the master list of bonds
  double* val;      // Values of the underlying bond length in the fitting set
  double* Ukernel;  // The contribution to the fitting matrix
  double* Ucontrib; // The kernel x stiffness constant = contribution to energy
  double* Fkernel;  // Same as Ukernel and Ucontrib above, but for forces.
  double* Fcontrib; //   These are significantly larger arrays, giving x, y,
                    //   and z contributions for each atom in the bond
};
typedef struct BondMap bondmap;

struct AngleMap {
  int nangl;        // The number of angles in this system
  aidx* id;         // Indices into the master list of angles
  double* val;      // Values of the underlying angle in the fitting set
  double* Ukernel;  // The contribution to the fitting matrix
  double* Ucontrib; // The kernel x stiffness constant = contribution to energy
  double* Fkernel;  // Same as Ukernel and Ucontrib above, but for forces.
  double* Fcontrib; //   These are significantly larger arrays, giving x, y,
                    //   and z contributions for each atom in the bond
};
typedef struct AngleMap anglmap;

struct TorsionMap {
  int ntterm;       // The number of torsional terms in this system
  hidx* id;         // Indices into the master list of torsion terms
  double* val;      // Values of the underlying angle in the fitting set
  double* Ukernel;  // The contribution to the fitting matrix
  double* Ucontrib; // The kernel x stiffness constant = contribution to energy
  double* Fkernel;  // Same as Ukernel and Ucontrib above, but for forces.
  double* Fcontrib; //   These are significantly larger arrays, giving x, y,
                    //   and z contributions for each atom in the bond
};
typedef struct TorsionMap tormap;

struct EnergyContributor {
  int fitcol;       // The fitting column that this packet pertains to
  double eave;      // The average energy
  double estd;      // Standard deviation in the energy
};
typedef struct EnergyContributor epacket;

struct AtomTypeSwitch {
  char orig[8];     // The original atom type name
  char pnew[8];      // The new atom type name
};
typedef struct AtomTypeSwitch typeswitch;

struct AtomTypeBranch {
  char instances[MAXLINE];  // The instances in which the atom type of name
                            //   orig is to be recast to new
  char orig[8];             // The original atom type name
  char pnew[8];              // The new atom type name
};
typedef struct AtomTypeBranch typebranch;

struct MMSystem {
  int GroupNum;    // Number of the topology group, computed in GroupSystems 
                   //   to correlate systems with similar topologies
  int PassedEtol;  // Flag to indicate whether this conformation passed the
                   //   check on conformational energy according to the
                   //   original MM force field (0 if not, 1 if the
                   //   conformation passed without problems, 2 if the
                   //   conformation passed after rearrangement)
  prmtop *tp;      // System topology pointer (points to tpencyc)
  coord crd;       // System coordinates
  int* atmap;      // Map of atom types in the system into the master list
  bondmap bmap;    // Bond terms mapping to a list spanning the entire fit
  anglmap amap;    // Angle terms mapping to a list spanning the entire fit
  tormap hmap;     // Torsion terms mapping to a list spanning the entire fit
  dmat excl;       // Matrix of exclusions; 0.0 for total exclusion (1:1
                   //   virtual site anchoring, 1:2 bonded, and 1:3 angle
                   //   interactions fall under this category), 1.0 for no
                   //   exclusion (nonbonded interactions 1:5 and more distal
                   //   have no exclusions), and any other value for a partial
                   //   exclusion (1:4 interactions)
  dmat nbnrg;      // The nonbonded energy matrix, electrostatics above the
                   //   diagonal and van-der Waals interactions below it
  double EEkernel; // The kernel for 1-4 adjustable electrostatic interactions 
  double LJkernel; // The kernel for 1-4 adjustable Lennard-Jones interactions 
  double EEnonfit; // The unfitted, unscaled non-bonded electrostatic energy
  double LJnonfit; // The unfitted, unscaled non-bonded Lennard-Jones energy
  double etrg;     // The target energy of the conformation, derived from
                   //   quantum calculations most likely
  double enorm;    // The "normalized" energy of this conformation--obtained by
                   //   subtracting from this system's etrg the average of etrg
                   //   from all conformations sharing the same topology 
  double eorig;    // The final energy of this conformation according to the
                   //   input Hamiltonian (topology file)
  double efin;     // The final energy of this conformation according to the
                   //   fitted parameters
  double nonfitmm; // The energy of non-adjustable molecular mechanics terms
  double wt;       // The weight that this conformation will get in the fit
  char crdsrc[MAXNAME];  // Source file for the coordinates (inpcrd format)
  char tpsrc[MAXNAME];   // The source file for the topology (prmtop format)
};
typedef struct MMSystem mmsys;

struct InstanceTracker {
  int sysid;       // The system ID number
  int sysno;       // The system number in the master list of conformations
  int order;       // The order of this instance (2, 3, or 4 for bonds, angles,
                   //   and dihedrals)
  int tnum;        // The number of the term within system structs for orders
                   //   2, 3, or 4
  char res[32];    // The residue names of atoms involved in this instance
  char atom[32];   // The atom names of atoms involved in this instance
};
typedef struct InstanceTracker itrack;

struct ParameterFit {
  int nconf;          // The number of conformations in this fitting set
  int natom;          //
  int nbond;          // The number of unique atom, bond, hydrogen bond, angle,
  int nhb1012;        //   and torsion terms in the parameter file and (if 
  int nangl;          //   supplied) in the frcmod file
  int ntor;           // 
  int ndihe;          // The number of proper dihedrals
  int nimpr;          // The number of impropers (ndihe + nimpr = ntor)
  int ncnst;          // The total number of constraints
  int nparm;          // The number of adjustable parameters, which must be no
                      //   greater than nbond + nangl + ntor + nscee + nscnb
  int nunisys;        // The number of unique systems (each topology file with
                      //   a unique name identifies a unique system)
  int nchng;          // The number of changes that have been made to atom
                      //   types in specific cases across all topologies
  int ljbuck;         // This parameter serves as a placeholder for the prmtop
                      //   struct attribute of the same name.  The value found
                      //   in this field will be passed down into every other
                      //   topology file in the conf array
  int FitAllBonds;    // Flags to activate fitting for all identified bond,
  int FitAllAngles;   //   angle, and torsion types
  int FitAllTorsions; //
  int reportall;      // Flag to indicate that all parameters encountered
                      //   should be reported
  int zeroNonfit;     // Flag to indicate that only fitted energy terms should
                      //   contribute to the molecular mechanics energy
  int nbadj;          // The numbers of adjustable bonds, angles, and torsions
  int naadj;          //   specified by the user
  int nhadj;          // 
  int nbvar;          // The number of adjustable bonds, angles, and torsions
  int navar;          //   actually found
  int nhvar;          // 
  int nzerocol;       // The number of zero-values columns in the fitting
                      //   matrix; this is generally not a problem as all
                      //   fitted parameters are restrained
  int ncorrcol;       // The number of highly correlated columns in the fitting
                      //   matrix
  int nrecast;        // The number of atom types to recast
  int ncleave;        // The number of atom types to branch (cleave)
  int fitscnb;        // Flag to activate fitting of the scnb (Lennard-Jones
                      //   1:4 scaling) term
  int fitscee;        // Flag to activate fitting of the scee (electrostatic
                      //   1:4 scaling) term
  int neqgroups;      // Number of equivalent atom groups
  int verbose;        // Display progress for the user (default 1, yes)
  int RemoveOutliers; // Flag to activate removal of conformations whose
                      //   energies are outliers (default 0, no)
  int* zerocol;       // List of zero-valued columns
  int* corrcol;       // List of pairs of highly correlated columns
  int* GroupCount;    // The number of each system, as differentiated by 
                      //   topologies
  int* FirstConf;     // The first conformations of each topology / system
  double lj14fac;     // The Lennard-Jones 1:4 scaling default value
  double elec14fac;   // The electrostatic 1:4 scaling default value
  double cnstB;       // The general constraints on bond, angle, and tosion
  double cnstA;       //   stiffnesses, to keep these values tethered to zero
  double cnstH;       //   and thus small in the result
  double cnst14;      // The 1:4 scaling term constraint factor (if 1:4 scaling
                      //   terms are being fitted)
  double mmtol;       // The molecular mechanics (bonded terms) energy
                      //   tolerance; configurations with energies higher than
                      //   this will be flagged for re-arrangement and reported
                      //   if the energy cannot be reduced
  double esigtol;     // Tolerance for the deviation of the total energy of any
                      //   particular conformation from the mean for all
                      //   conformations of that system
  double wtfloor;     // The minimum weight that any conformation may have;
                      //   default 0.5.  Set to smaller values to increase the
                      //   slant of the training set towards data points with
                      //   favorable energies.  Larger values make the training
                      //   set treat all conformations with equanimity.
  double* corrval;    // Correlations among highly aligned columns
  mmsys* conf;        // The system conformations, each with its own 
                      //   coordinates and topology
  xbonddef* badj;     // The array of adjustable bonds
  xangldef* aadj;     // The array of adjustable angles
  torterm* hadj;      // The array of adjustable torsions
                      //   (the arrays of adjustable bonds, angle, and torsions
                      //    will store the final fitted parameters at the end
                      //    of the fit)
  xatomdef* atoms;    // Atom definitions (masses and type names are included)
  xbonddef* bonds;    // Bond definitions (stiffnesses are fittable parameters
                      //   in a linear least-squares fit)
  xhb1012def* hb1012; // Hydrogen bondind 10-12 parameters (currently only kept
                      //   for the purposes of printing them out in the final
                      //   parameter file)
  xangldef* angls;    // Angle definitions (stiffnesses can be fitted by linear
                      //   least-squares optimization)
  torterm* torsions;  // Torsion terms (amplitudes are fitted by linear
                      //   least-squares optimization)
  dmat LJAcoef;       // Lennard-Jones A and B coefficient matrices
  dmat LJBcoef;       //
  cmat Hydrophilics;  // Hydrophilic atom types
  cmat ChangeLog;     // Log of changes to atom types
  eagrp* eqgroups;    // Equivalent atom groups
  prmtop* tpencyc;    // The encyclopedia of topologies in this fitting run
  typebranch* cleave; // Array of types that are to be branched
  typeswitch* recast; // Array of types whose names are to be recast
  char* ititl;        // Title for output parameter file
  char* icomm;        // Comment for fitted parameters
  char NrgUnits[64];  // Units of the energies for conformations supplied in
                      //   the input file
  char WaterName[8];  // The name of water molecules (for implementing SETTLE),
                      //   serves as a placeholder for the prmtop attribute of
                      //   the same name
  char ep[MAXNAME];   // The name of the EP rules file, a placeholder for the
                      //   prmtop eprulesource attribute; only one file may be
                      //   specified, as it will be passed down to all
                      //   topologies in the conf array
  char sr[MAXNAME];   // The series report file; if specified every torsion
                      //   series solved in the fit will be printed
  char ao[MAXNAME];   // The  accuracy output file; if specified a MatLab
                      //   script will be printed to display the results,
                      //   system by system
};
typedef struct ParameterFit prmset;

#endif
