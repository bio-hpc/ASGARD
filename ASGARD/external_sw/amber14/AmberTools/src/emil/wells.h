#ifndef HAVE_WELLS_H
#define HAVE_WELLS_H

#define WOBBLIUM_MAX_BONDS 28
#define WOBBLIUM_MIN_BONDS  6
#define WOBBLIUM_BOND_CUT   3.5

using namespace std;


class Well {

 public:
  Well            *next;
  double          *R, rCut, req, req2, scale;  
  double           rcut2, rcut2Inv, epsilon, rBond;
  Particle        *theParticle; // particle that is attracted to this well
  Particle        *firstParticle; // of list of particles in range of attraction
  int              npart;      // number of particles in range of attraction
  int              nChains;    // number of identical chains in range of attraction
  int              liquidIndex;
  Cell            *cell;       // cell of radius rcut
  
  //only relevant for wobblium...
  //  sub-classing is getting a bit stretched here.
  Well    **bondedWells;
  double   *bond_eq_l2, *bond_eq_l;
  double    bond_k;
  int       num_bonds;
  Particle *rootP;

  Well(){
    next          =   0;
    cell          =   0;
    firstParticle =   0;
    npart         =   0;
    rBond         =   1.0;
    epsilon       = 100.0;
    num_bonds     =   0;
    rootP         =   NULL;
  }

  void clean();
  
  //only relevant for wobblium...
  //  sub-classing is getting a bit stretched here.
  inline double pairDist2(Well *w, Box *b){
      double dx, dR;
      dR = 0.0;
      for( int j = 0; j < 3; j++){
        dx = w->R[j] - R[j];
        while( fabs(dx) > b->halfx[j] ) 
          dx -= copysign(b->x[j], dx);
        dR += dx * dx;
      }
      return(dR);
  }
  //only relevant for wobblium...
  //  sub-classing is getting a bit stretched here.
  inline double pairDist2_noImg(Well *w){
      double dx, dR;
      dR = 0.0;
      for( int j = 0; j < 3; j++){
        dx = w->R[j] - R[j];
        dR += dx * dx;
      }
      return(dR);
  }
  double computePhi(Box *b);
  double compute_dPhi(double *dPhiDr, Box *b);
  

  //overload these functions to give meaningful results in daughter classes
  virtual inline double computePhiOfDist(double rsq){cerr << "1Warning...called stub!\n";return(0.0);};
  virtual inline double compute_dPhiOfDist(double rsq){cerr << "2Warning...called stub!\n";return(0.0);};
  
  virtual        double reportFreeEnergy(){cerr << "3Warning...called stub!\n";return(0.0);};
  virtual        double reportFreeEnergy(double V){cerr << "4Warning...called stub!\n";return(0.0);};
  virtual        void   testMe(ostream *out){cerr << "5Warning...called stub!\n";};
  virtual        void   wobbliumAddBonds(Well **wells, int num_in, double r, Box *b){cerr << "6Warning...called stub!\n";};
};


//Use C++ inheritance to define separate types of well: conewell has a conical potential.
class coneWell: public Well {
  
  //takes rij^2 as input
  inline double computePhiOfDist(double rsq) {
    if (rsq > rcut2) return 0;
    else return (epsilon * ( sqrt(rsq*rcut2Inv) - 1.0 ) );
  };
  
  inline double compute_dPhiOfDist(double rsq) {
    if (rsq > rcut2) return 0.0;
    else return ( -epsilon );
  };
  
};

//This subclass defines a well which does nothing
class nullWell: public Well {
  
  using Well::reportFreeEnergy;

  //takes rij^2 as input
  inline double computePhiOfDist(double rsq) {
    return(0.0);
  };
  
  inline double compute_dPhiOfDist(double rsq) {
    return(0.0);
  };

  double reportFreeEnergy(double V){
    return(0.0);
  }
  
};
//harmonicWell inherits Well; and has a harmonic potential with a cutoff.
class harmonicWell: public Well {
 
  using Well::reportFreeEnergy;

  void  testMe( ostream *out ){
    double r;
    int    i;
    for( i = 0; i < 1000; i++ ){
      r = sqrt(rcut2) * double(i) / 1000;
     *out << i << " " << r << " " << computePhiOfDist(r * r) << " " << compute_dPhiOfDist( r * r ) << endl;
    }
  };
  
  //takes rij^2 as input
  inline double computePhiOfDist(double rsq) {
    if (rsq > rcut2) return 0.0;
    else return ( epsilon  * ( rsq * rcut2Inv - 1.0) );
  };

  inline double compute_dPhiOfDist(double rsq) {
    if (rsq > rcut2) return 0.0;
    else return ( -2.0 * epsilon * sqrt(rsq) * rcut2Inv);
  };
  
  double reportFreeEnergy(double V){
    double A, rootep, Z;
    double sfac, vfac, rcut3;
    
    //get some useful values
    rcut3  = rcut2 * sqrt(rcut2);
    rootep = sqrt(epsilon);
    
    //component due to vol outside the cutoff
    //vfac  = -log(V - 4.0 * M_PI * rcut3  / 3.0);
    
    //volume-dependent part
    vfac  = ( V - 4.0 * M_PI * rcut3  / 3.0 );
    
    //component due to vol inside cutoff
    //epsilon is already in units of kT.
    sfac  = M_PI * rcut3 * pow(epsilon,-1.5) * ( sqrt(M_PI) * exp(epsilon) * erf(rootep) - 2.0 * rootep );
    
    Z = vfac + sfac;
    
    //A = -kT ln Z.
    A  =  -1.0 * log(Z);
    
    return( A );//free energy in units of kT.
    
  };

};

//harmonicWell inherits Well; and has a harmonic potential.
class harmonicWellNoCut: public Well {
 
  using Well::reportFreeEnergy;

  //takes rij^2 as input... still mentions rcut, confusingly:
  //... this is just to have the zero of the potential consistently defined.
  inline double computePhiOfDist(double rsq) {
    return ( epsilon * ( rsq * rcut2Inv - 1.0) );
  };

  inline double compute_dPhiOfDist(double rsq) {
    return ( -2.0 * epsilon * sqrt(rsq) * rcut2Inv );
  };

  double reportFreeEnergy(){
    double A;
    
    //epsilon is already in units kt

    //harmonic well without cutoff has a manageable partition function:

    //Z = pow(M_PI * rcut2 / epsilon, 1.5) * exp( epsilon );
    //A  = -1.0 * log( Z ); this way is numerically unstable, use this--> 
    A  = -1.5 * ( log( M_PI * rcut2 ) - log( epsilon ) ) - epsilon;

    return( A );//free energy in units of kT.
    
  };
};

//harmonicBondWell inherits Well;
//and has a harmonic potential, relative to some root particle
class harmonicBondWell: public Well {

  using Well::reportFreeEnergy;


  //takes rij^2 FROM WELL MINIMUM as input
  inline double computePhiOfDist(double rsq) {
      return ( epsilon * ( rsq * rcut2Inv - 1.0) );
  };

  //takes rij^2 FROM WELL MINIMUM as input
  inline double compute_dPhiOfDist(double rsq) {
    double r;

    r = sqrt( rsq );
    return ( -2.0 * epsilon * r * rcut2Inv);
  };

  double reportFreeEnergy(){
    double A;

    //epsilon is already in units kt

    //Z = pow(M_PI * rcut2 / epsilon, 1.5) * exp( epsilon );
    //A  = -1.0 * log( Z ); this way is numerically unstable, use this--> 

    A  = -1.5 * ( log(M_PI) + log(rcut2) - log(epsilon) ) - epsilon;

    return( A );//free energy in units of kT.

  };
};
//winged harmonicWell inherits Well; and has a harmonic potential which switches to linear at large r.
class wingWell: public Well {
 
  using Well::reportFreeEnergy;
  public:
  
  double depth, reqInv, reqInv2;
  double wingForce, scalex2;
  
  wingWell(double req_in, double rcut_in, double target_depth, double wingForce_in ){
    
    rCut     = rcut_in;
    rcut2    = rCut * rCut;
    
    req           = req_in;  //r_eq is defined as the position such that U = U_min + kT/2.
    req2          = req * req;
    reqInv        = 1.0 / req;
    reqInv2       = reqInv * reqInv;
    
    wingForce         = -1.0 * wingForce_in;
    depth             = -1.0 * target_depth;
    scale             = wingForce * (rCut - req) - depth; //use this scale factor to modify the harmonic part only.  
    scalex2           = -2.0 * scale;
  }
  
  void report( ostream *out ){
              *out << "#   k1:         " 
                             << -1 * scalex2 * reqInv2 << "\n";
              *out << "#   r1:         " << sqrt(req2)  << "\n";
              *out << "#   wing force: " << -1*wingForce   << "\n";
              *out << "#   r2:         " << sqrt(rcut2) << "\n";
              *out << "#   depth:      " << -1*depth       << "\n"; 
  }

  void  testMe( ostream *out ){
    double r;
    int    i;
   *out << "#k1: "    << scalex2 * reqInv2 << " r1: " << sqrt(req2) 
        << " k2: "    << wingForce         << " r2: " << sqrt(rcut2) 
        << " depth: " << depth << endl;    
    for( i = 0; i < 1000; i++ ){
      r = rCut * double(i) / 1000;
     *out << i << " " << r << " " << computePhiOfDist(r * r) << " " << compute_dPhiOfDist( r * r ) << endl;
    }
  };
  
  //takes rij^2 as input
  inline double computePhiOfDist(double rsq) {
    if (rsq > rcut2) return ( 0.0 );//flat part
    else if( rsq > req2 ) return ( wingForce * ( rCut - sqrt(rsq) ) ); //linear part
    else return ( depth + scale * rsq * reqInv2 );    //harmonic part
  };

  inline double compute_dPhiOfDist(double rsq) {
    if (rsq > rcut2) return 0.0;
    else if( rsq > req2 ) return ( wingForce ); //linear part
    else return ( scalex2 * sqrt(rsq) * reqInv2 );//harmonic part
  };
  
  double reportFreeEnergy(double V){
    
    double Z_3, A, r, sfac_harm, vfac, rcut3;
    double deps, delta, eps_c, req3, root_deps, poly1, poly2;    


    //Z_2 is the partition func of the harmonic part:
    //for this we use the formula for a well with a cut, 
    //then just add an offset later.
    req3      = req * req * req;
    deps      = fabs(scale);  //depth of the harmonic part < total depth.
    root_deps = sqrt(deps);
    sfac_harm = M_PI * req3 * pow(deps,-1.5) * ( sqrt(M_PI) * exp(deps) * erf(root_deps) - 2.0 * root_deps );
    
    //offset of the depth of the harmonic part
    delta     = depth + scale;
    sfac_harm = sfac_harm * exp(delta); //constant offset of the well factors out of the partition function.

    //get some useful values
    r      = sqrt(rcut2);
    rcut3  = rcut2 * r;
    
    //Z_3 is the partition function for the conical part:
    eps_c  = fabs(wingForce);
    poly1  = req2  * eps_c * eps_c + 2.0 * req  * eps_c + 2.0; 
    poly2  = rcut2 * eps_c * eps_c + 2.0 * r    * eps_c + 2.0; 

    Z_3   = -4.0 * M_PI * (exp(-1*eps_c*req)*poly1 - exp(-1*eps_c*r)*poly2) / (eps_c*eps_c*eps_c);

    //volume-dependent part
    vfac = ( V - 4.0 * M_PI * rcut3  / 3.0 );
    
    //partition function
    A = -1.0 * log(  sfac_harm + Z_3 + vfac );   
    
    return( A );//free energy in units of kT.
    
  };

};
//wobblium: the whole system (solvent+solute) is an 
//elastic network
class wobbliumWell: public Well {
  

 public : wobbliumWell(){
    bondedWells  = new Well*[WOBBLIUM_MAX_BONDS];
    bond_eq_l    = new double[WOBBLIUM_MAX_BONDS];
    bond_eq_l2   = new double[WOBBLIUM_MAX_BONDS];
    num_bonds    = 0;
    for( int i =  0 ; i < WOBBLIUM_MAX_BONDS; i++) {
      bondedWells[i] = NULL;
    }
  }
  
  double reportFreeEnergy(){
    return(0.0); //need to write this still.
  }
  double reportFreeEnergy(double V){
    return(0.0); //need to write this still.
  }

  void wobbliumAddBonds(Well **bondedWells_in, int num_in, double bond_k_in, Box *b){
    
    //save the spring constant
    bond_k    = bond_k_in;
    num_bonds = num_in;

    //build up the bond list.
    for(int i = 0; i < num_bonds; i++){
      bondedWells[i] = bondedWells_in[i];
      bond_eq_l2[i]  = pairDist2(bondedWells[i], b);
      //bond_eq_l2[i]  = pairDist2_noImg(bondedWells[i]);
      bond_eq_l[i]   = sqrt(bond_eq_l2[i]);
    }
  }


  inline double computePhiOfDist(double rsq) {
      cerr << "Warning, called stub! computePhiOfDist()" << endl;
      return( 0.0 );
  };
  inline double compute_dPhiOfDist(double rsq) {
      cerr << "Warning, called stub! compute_dPhiOfDist()" << endl;
      return( 0.0 );
  };

};
#endif //HAVE_WELLS_H
