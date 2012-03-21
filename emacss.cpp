/***********************************************************************/
//  Evolve Me A Cluster of StarS (EMACSS)                              //
//                                                                     //
//  By: Poul Alexander (University of Cambridge)                       //
//      Mark Gieles (University of Cambridge)                          //
//                                                                     //
//  EMACSS is a numerical integrator that solves for the differential  //
//  equations:                                                         //
//                                                                     //
//  dN/dt = f(N,r,p) = - xi*N/t_rh                                     // 
//                                                                     //
//  dr/dt = g(N,r,p) = (zeta - 2*xi) * r/t_rh                          //
//                                                                     // 
//  where                                                              // 
//                                                                     // 
//  xi = f(N,r,rt)                                                     //
//  p is an array of parameters/variables, which define the properties //
//  of star clusters at any given time in the life-cycle               //
//                                                                     //
// Equations referenced AG2012 refer to Alexander & Gieles 2012        //
/***********************************************************************/
#include "emacss.h"

int main(int argc, char* argv[]){
  node cluster;
  variables init;
  parameters params;

  readinput(&init,argc,argv);
  get_params(&params);
  cluster.initialise(params,init);
  cluster.core_collapse();
  while (cluster.N > 200){
    cluster.evolve();
    cluster.output();
  }
  return 0;
}

//----------------------------------------------------------------------------//
//Initialisation functions
	  
node::node(){ //Iterator for a cluster at any given time(initialised)
}

void node::initialise(parameters in1, variables in2){
  //Assigns initial properties for a cluster
  params = in1;
  init = in2;

  M_gal = init.galaxy.M;                   //Sets galaxy properties
  R_gal =  init.galaxy.R;
  if (init.units==1) set_units();
  
  l1 = l2 = l3 = l4 = l5 = true;           //Sets flags for lazy evaluation
  
  t = 0;                                   //Assigns first data point
  N = Nstart = y[0] = trial[0] = init.N0;
  r = y[1] =  trial[1] = 1.0; 
  mm = y[2] =  trial[2] = 1.0/init.N0;

  zeta = init.zeta;                        //Other user specifications
  toff = init.tcc;

  n_relax = 0;                             //Other initialisations
  first = true;
  interp = 1;
  output();
}

void node::set_units(){
  G_star = 0.00449857;                    //Grav constant (pc^3M_sun^-1Myr^-2)
  M_star = init.N0*init.mm0;              //Initial Mass of cluster (M_sun)
  R_star = init.r0;                       //Virial radii (parsec / N-body)
  T_star = sqrt(pow(R_star,3)/(M_star*G_star)); //N-body time (Myr)
  
  //Conversion of Galaxy M and R to N-body units
  M_gal = init.galaxy.M/(init.mm0*init.N0);     
  R_gal =  init.galaxy.R/R_star;
}

//---------------------------------------------------------------------------//
//Evolution functions

void node::evolve(){
  /*Does what it says in the title  - evolves by rk4 kernel.*/

  tstep = params.frac*trh();
  rk4();
  N = y[0]; r = y[1]; mm = y[2]; t += tstep;
  n_relax = y[3];
}

void node::solve_diffs(double deriv[]){
  /*Differential equations in usable form for the rk4 routine.*/

  l1 = l2 = l3 = l4 = l5 = true; //Sets flags for lazy evaluation
  deriv[0] = dNdt();
  deriv[1] = drdt();
  deriv[2] = dmmdt();
  deriv[3] = dtrhdt();
}

/*---------------------------------------------------------------------------*/
// Defining differential equations (see AG2012, section 2)
double node::dNdt(){
  return -(xi()*trial[0])/trh();                   //Equation (6) AG2012
}

double node::drdt(){
  return (mu()*trial[1])/trh();                    //Equation (7) AG2012
}

double node::dmmdt(){
  return 0.0;                                      //To be implemented in future
}

double node::dtrhdt(){
  return 1.0/trh();                                //t_rh per unit time
}

/*---------------------------------------------------------------------------*/
//Defining Physical Parameter equations (see AG2012, sections 2 & 2.2)

double node::r_jacobi(){                            //Equation (18) AG2012
  static double rj;
  if (l1) {
    l1 = false;
    rj = pow((trial[0]*trial[2])/(3.0*M_gal),(1.0/3.0))*R_gal;
  }
  return rj;
}

double node::trh(){                                  //Equation (4) AG2012
  static double trh;
  if (l2) {
    l2 = false;
    trh = 0.138*sqrt(trial[0]*pow(trial[1],3)/(trial[2]))/
      (log(params.gamma*trial[0]));
  }
  return trh;
}

/*---------------------------------------------------------------------------*/
//Defining Dimensionless Parameter equations (see AG2012, sections 2 & 3)

double node::xi(){                                   //Equation (26) AG2012
  static double xi;
  if (l4){
    l4 = false;
    double tidal = (3.0/5.0)*zeta;
    double iso = params.xi0;
    xi = tidal*P()+iso*(1.0-P());
  }
  return xi;
}

double node::mu(){                                    //Equation (8) AG2012
  static double mu;
  if (l5){
    l5 = false;
    mu = zeta-2.0*xi();
  }
  return mu;
}

double node::energy(){                                //Equation (1) AG2012
  return -0.25*pow(trial[0]*trial[2],2)/trial[1];
}

double node::P(){                                     //Equation (25) AG2012
  static double P;
  if (l3){
    l3 = false;
    double f = pow((double)(trial[0]*log(params.gamma*params.N1))/
		 (params.N1*log(params.gamma*trial[0])),(double)(1.0-params.x));
    double g = pow((double)((trial[1]/r_jacobi())/params.rhrj1),
		   (double)params.z);
    P = f*g;
  }
  return P;
}

/*---------------------------------------------------------------------------*/
//Numerical integrator (standard rk4 routine, see documentation anywhere)
void node::rk4(){
  
  for (int i=0;i<4;i++) trial[i] = y[i];
  solve_diffs(dr1);
  for (int i=0;i<4;i++) trial[i]=y[i]+(0.5*tstep*dr1[i]);
  solve_diffs(dr2); 
  for (int i=0;i<4;i++) trial[i]=y[i]+(0.5*tstep*dr2[i]);
  solve_diffs(dr3);
  for (int i=0;i<4;i++) trial[i]=y[i]+(tstep*dr3[i]);
  solve_diffs(dr4);
  for (int i=0;i<4;i++) 
    y[i]=y[i]+(tstep/6.0)*(dr1[i]+2.0*dr2[i]+2.0*dr3[i]+dr4[i]); 
}

//Here endeth the code.
