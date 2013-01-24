/**************************************************************/
/*Headers - not all actually needed at present but here anyway*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <cstring>
#include <cmath>
#include <new>

#include <stdlib.h>
#include <stdio.h>

using namespace std;
/**************************************************************/

/**************************************************************/
/*Declares structuresand the parameter classes*/

typedef struct{  //Time information - output, N-body, real, elapsed relaxation
  double nbody, Myr;
} t;

typedef struct{  //Mass information - N-body, real
  double nbody, Msun;
} mass;

typedef struct{  //Radius information - N-body, real
  double nbody, pc;
} radius;

typedef struct{  //Radius information - N-body, real
  double nbody, kms;
} velocity;

typedef struct{ //Energy information - output, steady state, generation, 
  double epsilon, zeta, nbody, real;
  int source;
} energy;

typedef struct{  //The properties defining the tidal field
  int type;
  mass M; radius R; velocity v;
} tidal_field;

class stellar_evo;
class dynamics;

class node{                       //Binding of cluster paramters at given time
  void zero();
  void initialise();
  void solve_odes(double[],stellar_evo,dynamics);
  void convert();
  double E_calc(), trh(), r_jacobi();                      //Calculated factors
  tidal_field galaxy            ;                          //Galaxy conditions
  double G_star, M_star, R_star, T_star;                   //Conversion factors
  double frac, tstep, R;         //Timesteps, input
  double *nbody[10], *real[10];  //Data arrays
 public:
  node();
  double gamma; 
  int s, units;
  void input(int, char*[]);
  void evolve(stellar_evo,dynamics);
  void output();
  t time, out_time, t_relax;
  energy E;
  double N, kappa, trhelapsed;
  mass M_total, DM_SE;
  radius r, rj;
};

/*The following is the stellar evolution module - models stellar effects*/
class stellar_evo{
  double y[];
  node *mynode;
 public:
  stellar_evo(node*);
  double nu, T_SE, X;              //defining characters (set at intit)
  void set_X();
  double dEdt(), dmsedt();
  double epsilon(), gamma_se(); 
};

/*The following is the dynamics module - models pure dynamical effects*/
class dynamics{
  double P();
  double y[];
  node *mynode;
  stellar_evo *myse;
 public:
  dynamics(node*,stellar_evo*);
  double R1, N1, x, z, xi0, T_DYN, a, F;   //defining characters (set at intit)
  void set_tcc();
  double dNdt(), dMdt(), drdt(), dkdt(), dtrhdt();
  double xi(), gamma_dyn(), mu(), lambda();
};

/**************************************************************/
/*Non-physical function calls (i.e., program tidying, help file etc.)*/
void help();
void version();
/**************************************************************/
