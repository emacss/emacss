/**************************************************************/
/*Headers - not all actually needed at present but here anyway*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <cstring>
#include <cmath>
#include <limits>

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
  double zeta, nbody, real;
  int source;
} energy;

typedef struct{  //The properties defining the tidal field
  int type;
  mass M; radius R; velocity v;
} tidal_field;

class stellar_evo;
class dynamics;

class node{                       //Binding of cluster paramters at given time
  void check_input();
  void initialise();
  void solve_odes(double[],stellar_evo,dynamics);
  void convert();
  double E_calc(), trh(), r_jacobi(), t_crj();             //Calculated factors       ;                          //Galaxy conditions
  double G_star, M_star, R_star, T_star;                   //Conversion factors
  double frac, R;         //Timesteps, input
  double *nbody[13], *real[13];  //Data arrays
 public:
  node();
  
  double gamma, rhrj, tstep, mass_seg; 
  int s, units;
  void input(int, char*[]);
  void zero();
  void evolve(stellar_evo,dynamics);
  void output(stellar_evo,dynamics);
  t time, out_time, t_relax, tcrj;
  tidal_field galaxy;     
  energy E;
  double N, kappa, k_0, trhelapsed, trhp, T_DYN, T_SE;
  mass mm, mm_se;
  radius r, rj;
};

/*The following is the stellar evolution module - models stellar effects*/
class stellar_evo{        
  node *mynode;
  dynamics *mydyn;
  double tidal_escape(), pcc();
 public:
  stellar_evo();
  void load(node*,dynamics*);
  double  nu, MS_1, y, z0, esc_frac;
  double dEdt(), dmsedt(), dmsegdt(); 
  double epsilon(), gamma_se(), chi(), MS_max(), zeta();
};

/*The following is the dynamics module - models pure dynamical effects*/
class dynamics{
  double P();
  double f_delay(), f_max(), td;
  node *mynode;
  stellar_evo *myse;
 public:
  dynamics();
  void load(node*,stellar_evo*);
  double R1, N1, x, z, xi0, F, k_1, n, f; //defining characters (set at intit)
  double dNdt(), dmdyndt(), drdt(), dkdt();
  double dtrhdt(), dtrhpdt(), dmmdt();
  double xi(), xi_ind(), gamma(), gamma_dyn(), mu(), lambda();
};

/**************************************************************/
/*Non-physical function calls (i.e., program tidying, help file etc.)*/
void help();
void version();
/**************************************************************/
