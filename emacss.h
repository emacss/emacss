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
/*Declares and the parameter classes*/

typedef struct{  //The properties defining the tidal field
  double M, R, v;  
} tidal_field;

typedef struct{  //Fixed Parameters to describe the evolution
  double xi0, rhrj1, gamma, N1, x, z, frac;
} parameters;

typedef struct{  //Cluster variables 
  double N0, r0, mm0;                 //Initial properties
  double f_N, f_r, f_mm;              //Pre-core collapse evoltuion
  double zeta, rjrh, tcc;             //Other properties
  double tout;                        //Output time (if wanted)
  int units, s;                       //User unit choice, stellar evoltuiom.
  bool s1,s2,s3;                      //Checks which variables are set by user.
  tidal_field galaxy;
} variables;

class node{
  bool first, virial;              //Checks output options
  bool l1, l2, l3, l4, l5, l6;     //Lazy Evaluation flags
  int interp, s;                   //Interpolation flag
  double y[4], yerr[4];            //yerr not used unless error checking (rk5)
  double dr1[4], dr2[4],dr3[4],dr4[4],dr5[4],dr6[4],trial[4];
  double Nstart, tstep, collapse_time, M_gal, R_gal, trh_old;
  double G_star, T_star, M_star, R_star;
  double P();
  double xi();
  double sigma();	
  double mu();
  double dNdt();
  double drdt();
  double dmmdt();
  double dtrhdt();
  void set_units();
  void rk4();
  void solve_diffs(double[]);
  double exp_get_c(int);
  double exp_get_b(int);
  double poly_get_a(int);
  double poly_get_b(int);
  void pre_collapse();
  void precc_evolve(double[]);
  void check_post_collapse_state(float);
  void end_collapse();
  variables init;
  parameters params;
 public:
  node();
  double t, N, r, mm, zeta, n_relax;
  void initialise(parameters,variables);
  double trh();
  double r_jacobi();
  double energy();
  void prepare_interpolators();
  void core_collapse();
  void evolve();
  void output_time(double);
  void output();
};
/*************************************************************/

/*************************************************************/
/*Input functions declarations*/
void readinput(variables *, parameters *, int, char*[]);
void readfile(char[],variables *);
void zero(variables *);
void set_cluster(variables *);
void set_galaxy(variables *);
void early_evolution(variables *, float);
void galaxy_r(variables *);
void galaxy_rhrj(variables *);
void galaxy_vel(variables *);
void galaxy_all(variables *);
double r_j(variables *);
double trh(variables *);
void get_params(parameters *);
void help();
void version();
/*************************************************************/
