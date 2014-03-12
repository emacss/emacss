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
  int type, f;
  double M, R, R2, v;
} tidal_field;

class stellar_evo;
class dynamics;

class node{                       //Binding of cluster parameters at given time
  void solve_odes(double[],stellar_evo,dynamics);
  void convert(stellar_evo,dynamics);
  double E_calc(), trh(), trhp(), r_jacobi();              //Calculated factors
  double trc(), rhoc();                                    //Relaxation times
  double step_min();                                      //Minimum stepsize
  double G_star, M_star, R_star, T_star, pcMyr;            //Conversion factors
  double tstep;                                            //Timesteps
  double *nbody[13];                                       //Data arrays 
  int nvar;
 public:
  node();
  int units, s;
  void input(int, char*[]);
  void initialise(stellar_evo,dynamics);
  void zero();
  void load(dynamics*,stellar_evo*);
  void evolve(stellar_evo,dynamics);
  void output(stellar_evo,dynamics);
  tidal_field galaxy;                                      //Galaxy conditions
  double time, t_rh, t_rhp, t_rc, tcc, out_time, tdf;      //Times
  energy E;                                                //Energies
  double N, kappa, trhelapsed, trhpelapsed, MS, frac, psi; //Dimensionless
  double mm, mm_se, m_max, m_max0,  m_min;                 //Masses
  double rh, rv, rj, rc;                                   //Radii
  double Rhj, Rch;                                         //Ratios
  double gamma, k0, k1, Rch0;                              //Changeable parameters
};

/*The following is the stellar evolution module - models stellar effects*/
class stellar_evo{     
  double tidal_escape(), p();
  node *mynode;
  dynamics *mydyn;
  double  nu, y, esc_frac;
  double m_ns, m_ref, t_ref, t_se;
  double psi1, psi0;                                        //Ratios 
 public:
  stellar_evo();
  void load(node*,dynamics*);
  void tse(double, double);
  double MS_1;
  //Differential equations
  double dEdt(), dmsedt(), dmsegdt(); 
  //Dimensionless factors
  double epsilon(), gamma_se(), chi(), zeta(), psi(), m_up();
};

/*The following is the dynamics module - models pure dynamical effects*/
class dynamics{
  double P(), K(), k(), f_ind(), F();   //Tidal Factor (AG2012, eq: 25)
  double y[];
  node *mynode;
  stellar_evo *myse;
  double R1, N1, x, z, xi1, f, t_df, Fej;    //Tidal characteristics (set at start)
  double delta_1, delta_2;                       //Core characteristics
  double N2, N3;                                 //Core Scalings
  double Y, b, q;                                //Misc
 public:
  dynamics();
  void load(node*,stellar_evo*);
  void  reset_K_constants(), set_tcc(), set_k0(); //Functions
  void tdf(double, double, double, double);
  double nc;                                     //Needed for energy in myse
  double Rchmin(), m_ej();
  //Differential equations
  double dNdt(), dmmdt(), drdt(), dkdt(), drcdt();  
  double dtrhdt(), dtrhpdt(), dr2dt();
  //Dimensionless factors
  double xi(), xi_e(), xi_i(), gamma(), gamma_dyn();
  double mu(), lambda(), delta(), epsilon();    
};

/**************************************************************/
/*Non-physical function calls (i.e., program tidying, help file etc.)*/
void help();
void version();
/**************************************************************/
