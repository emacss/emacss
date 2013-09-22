/*Cluster module - contains the general functions for the cluster,such 
 * as the initialiser, exoltuion modules, and output functionality*/

#include "../emacss.h"

double node::E_calc(){                                //Equation (1) AG2012
  return -kappa*pow(N*mm,2)/rh;
}

double node::r_jacobi(){                             //Equation (18) AG2012
  double rj = 0;
  if (galaxy.type == 1)
    rj = pow((N*mm)/(3.0*galaxy.M),(1.0/3.0))*galaxy.R;
  else if (galaxy.type == 2)
    rj = pow((N*mm)/(2.0*galaxy.M),(1.0/3.0))*galaxy.R;
  return rj;
}

double node::trh(){                                   //Equation (4) AG2012
  return 0.138*sqrt(N*pow(rh,3)/mm)/(log(gamma*N));
}

double node::trc()                                    //Equation (14) GALB2013
{
  double sigc2 = 8./3.*M_PI*rhoc()*pow(rc,2.0); // Core velocity dispersion
  return pow(sigc2,1.5)/(15.4*mm*rhoc()*log(0.11*N));
}

double node::rhoc() // MGi: compute core density
{
  double rhoc = 0;
  if (E.source == 0)
    rhoc = 0.055*pow(rc,-2.2);

  if (E.source == 1){
    double rhoh = 3.0*N*mm/(8*M_PI*pow(rh,3.0));
    rhoc = pow(rc/rh,-2.0)*rhoh;
  }
  return rhoc;
}

void node::evolve(dynamics dyn){
  /*Does what it says in the title  - evolves by rk4 kernel.*/
    
  static double duplicate_nbody[11]; 
  static double dr1[11],dr2[11],dr3[11],dr4[11]; 
	
  if (E.source == 0) tstep = 1.0/(E.zeta/(10*t_rc) + E.zeta/(0.01*t_rh));
  else tstep = 0.01*t_rh/E.zeta; 

  
  for (int i=0;i<nvar;i++){
    duplicate_nbody[i] = *nbody[i]; //Backup
  }
  
  solve_odes(dr1,dyn);
  convert();
  for (int i=0;i<nvar;i++) 
    *nbody[i]=duplicate_nbody[i]+(0.5*tstep*dr1[i]);
  solve_odes(dr2,dyn);
  convert();

  for (int i=0;i<nvar;i++) 
    *nbody[i]=duplicate_nbody[i]+(0.5*tstep*dr2[i]);
  solve_odes(dr3,dyn);
  convert();

  for (int i=0;i<nvar;i++) 
    *nbody[i]=duplicate_nbody[i]+(tstep*dr3[i]);
  solve_odes(dr4,dyn);
  convert();

  for (int i=0;i<nvar;i++) 
    *nbody[i]=duplicate_nbody[i]+(tstep/6.0)*(dr1[i]+2.0*dr2[i]+2.0*dr3[i]+dr4[i]);
  convert(); 

  // Switch energy source on at core collapse. 
  if  (Rch <= dyn.Rchmin() && E.source == 0){
      tcc = time;                            //Stores core collapse time
      dyn.reset_K_constants();               //Changes parameters
      E.source = 1;                          //Changes energy source
    }
}

void node::solve_odes(double dvdt[],dynamics dyn){
/*Essentially a wrapper function; wraps the solving of all requisite ODES into 
 * one function. Calls to a variety of dynamics module functions.*/	
   
    dvdt[0] = 1;
    dvdt[1] = dyn.dEdt();
    dvdt[2] = dyn.dNdt();
    dvdt[3] = dyn.dmmdt();
    dvdt[4] = dyn.drdt();
    dvdt[5] = 0;
    dvdt[6] = 0;
    dvdt[7] = dyn.dkdt();
    dvdt[8] = dyn.dtrhdt();
    dvdt[9] = dyn.drcdt(); 
    dvdt[10] = 0;
}

void node::convert(){
 /*Housekeeping function. Calls functions to determine the quantities that need 
  * to be redefined*/
  
    t_rh = trh(); t_rc = trc(); rj = r_jacobi(); rv = rh/(4.0*kappa);
    Rhj = rh/rj; Rch = rc/rh;
}
