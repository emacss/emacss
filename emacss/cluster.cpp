/*Cluster module - contains the general functions for the cluster,such 
 * as the initialiser, exoltuion modules, and output functionality*/

#include "../emacss.h"

double node::E_calc(){                                //Equation (1) AG2012
  return -0.25*pow(N*mm.nbody,2)/r.nbody;
}

double node::r_jacobi(){                            //Equation (18) AG2012
  double rj = 0;
  if (galaxy.type == 0)
	rj = pow((N*mm.nbody)/(3.0*galaxy.M.nbody),(1.0/3.0))*galaxy.R.nbody;
  else if (galaxy.type == 1)
	rj = pow((N*mm.nbody)/(2.0*galaxy.M.nbody),(1.0/3.0))*galaxy.R.nbody;
  return rj;
}

double node::trh(){                                  //Equation (4) AG2012
return 0.138*sqrt(N*pow(r.nbody,3)/(mm.nbody))/(log(gamma*N));
}

void node::evolve(stellar_evo se,dynamics dyn){
  /*Does what it says in the title  - evolves by rk4 kernel.*/
  static double duplicate_nbody[10]; 
  static double dr1[10],dr2[10],dr3[10],dr4[10]; 

  tstep = 1.0/(1.0/(frac*t_relax.nbody)+1.0/(frac*time.nbody/se.nu));
  if (tstep ==0 ) tstep = 0.1;       //incase t = 0
  
  for (int i=0;i<9;i++) duplicate_nbody[i] = *nbody[i]; //Backup
  
  solve_odes(dr1,se,dyn);
  convert();
  for (int i=0;i<10;i++) *nbody[i]=duplicate_nbody[i]+(0.5*tstep*dr1[i]);
  solve_odes(dr2,se,dyn);
  convert();
  for (int i=0;i<10;i++) *nbody[i]=duplicate_nbody[i]+(0.5*tstep*dr2[i]);
  solve_odes(dr3,se,dyn);
  convert();
  for (int i=0;i<10;i++) *nbody[i]=duplicate_nbody[i]+(tstep*dr3[i]);
  solve_odes(dr4,se,dyn);
  convert();
  for (int i=0;i<10;i++) 
  *nbody[i]=duplicate_nbody[i]+(tstep/6.0)*(dr1[i]+2.0*dr2[i]+2.0*dr3[i]+dr4[i]);
  convert(); 
  
  
}

void node::solve_odes(double dvdt[],stellar_evo se,dynamics dyn){
/*Essentially a wrapper function; wraps the solving of all requsite differential
 equations into one function. Calls to a variety of dynamics and se functions.*/	
   
    dvdt[0] = 1;
    dvdt[1] = 0;
    dvdt[2] = dyn.dNdt();
    dvdt[3] = dyn.dmmdt();
    dvdt[4] = dyn.drdt();
    dvdt[5] = 0;
    dvdt[6] = 0;
    dvdt[7] = se.dmsedt();
    dvdt[8] = dyn.dkdt();
    dvdt[9] = dyn.dtrhdt();
}

void node::convert(){
 /*Housekeeping function. Keeps n-body and real units aligned when called.
  Also calls functions to determine the quantities that need to be redefined*/
  
    E.nbody = E_calc(); E.real = E.nbody*pow(M_star,2)/R_star;
    t_relax.nbody = trh(); t_relax.Myr = t_relax.nbody*T_star;
    rj.nbody = r_jacobi(); rj.pc = rj.nbody*R_star;
    
    *real[0] = *nbody[0]*T_star;
    *real[2] = *nbody[2];
    *real[3] = *nbody[3]*M_star;
    *real[4] = *nbody[4]*R_star;
    *real[7] = *nbody[7]*M_star;
    *real[8] = *nbody[8];
    *real[9] = *nbody[9];
}
