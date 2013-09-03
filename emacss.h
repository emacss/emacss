/**************************************************************/
/*Headers - not all actually needed at present but here anyway*/
#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <cstring>
#include <cmath>

#include <stdlib.h>
#include <stdio.h>

using namespace std;
/**************************************************************/

/**************************************************************/
/*Declares structures and the parameter classes*/

typedef struct{ //Energy information - output, steady state, generation, 
  double zeta, value;
  int source;
} energy;

typedef struct{  //The properties defining the tidal field
  int type;
  double M, R, v;
} tidal_field;

class dynamics;

class node{                       //Binding of cluster parameters at given time
  void initialise();
  void solve_odes(double[],dynamics);
  void convert();
  double E_calc(), trh(), r_jacobi();                      //Calculated factors
  double trc(), rhoc();                                    //Relaxation times
  double G_star, M_star, R_star, T_star;                   //Conversion factors
  double frac, tstep;                                      //Timesteps
  double *nbody[11];                                       //Data arrays 
  int nvar;
 public:
  node();
  int units;
  void input(int, char*[]);
  void zero();
  void evolve(dynamics);
  void output(dynamics);
  tidal_field galaxy;                                      //Galaxy conditions
  double time, t_rh, t_rc, tcc;                            //Times
  energy E;                                                //Energies
  double N, kappa, trhelapsed;                             //Dimensionless
  double mm;                                               //Masses
  double rh, rv, rj, rc;                                   //Radii
  double Rhj, Rch;                                         //Ratios
  double gamma, k0, k1, Rch0;                              //Changeable parameters
};

/*The following is the dynamics module - models pure dynamical effects*/
class dynamics{
  double P();                             //Tidal Factor (AG2012, eq: 25)
  double K();                             //Core Factor (GALB2013, eq: 29)
  double k();
  double y[];
  node *mynode;
  double R1, N1, x, z, xi1, f;            //Tidal characteristics (set at start)
  double delta_1, delta_2;                //Core characteristics
  double N2, N3;                          //Core Scalings
 public:
  dynamics(node*);
  
  void  reset_K_constants(), set_tcc(), set_k0(); //Functions
  double Rchmin();
  double dNdt(), dmmdt(), drdt(), dkdt(), dtrhdt(), drcdt(), dEdt();  
  //Differential equations
  double xi(), gamma_dyn(), mu(), lambda(), delta(), epsilon();    
  //Dimensionless factors
};

/**************************************************************/
/*Non-physical function calls (i.e., program tidying, help file etc.)*/
void help();
void version();
/**************************************************************/
