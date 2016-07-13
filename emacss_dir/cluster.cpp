/*Cluster module - contains the general functions for the cluster,such 
 * as the initialiser, exoltuion modules, and output functionality*/

#include "../emacss.h"

static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
        b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
	b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
	b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
	b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
	c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
	dc5 = -277.00/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
	 dc4=c4-13525.0/55296.0,dc6=c6-0.25;

double node::E_calc(){                              //Equation (1) AG2012
  double E;
  
  E = -kappa*pow(N*mm,2)/rh;
  
  return E;
}

double node::r_jacobi(){                             //Equation (18) AG2012
  double rj = 0;
  
  if (s == 0)
    rj = pow((N*mm)/(3.0*galaxy.M),(1.0/3.0))*galaxy.R;
  else if (s == 1)
    rj = pow((N*mm)/(2.0*galaxy.M),(1.0/3.0))*galaxy.R;
  
  return rj;
}

double node::trh(){                                   //Equation (4) AG2012
 double trh = 0;	
 trh = 0.138*sqrt(N*pow(rh,3)/mm)/(log(gamma*N));
// if (time < 2) cerr << time*T_star << ' ' << trh*T_star << endl;
 
 return trh;
}

double node::trhp(){                                  //Equation (4) AG2012
double trhp;

trhp =  t_rh/psi;

return trhp;
}

double node::trc(){                                    //Equation (14) GALB2013
double trc = 0, sigc2 = 0;
sigc2 = 8./3.*M_PI*rhoc()*pow(rc,2.0); // Core velocity dispersion
trc = pow(sigc2,1.5)/(15.4*mm*rhoc()*log(0.11*N));
return trc;
}

double node::rhoc(){ // MGi: compute core density
  double rhoc = 0;
  if (E.source == 0)
    rhoc = 0.055*pow(rc,-2.2);

  if (E.source == 1){
    double rhoh = 3.0*N*mm/(8*M_PI*pow(rh,3.0));
    rhoc = pow(rc/rh,-2.0)*rhoh;
  }
  return rhoc;
}

double node::step_min(){                           //Ensures progress!
   double min = 0;
   if (s == 0) min = 1.0/(1.0/(frac*t_rc)+1.0/(frac*t_rh));
   else min = 1.0/(1.0/(frac*t_rhp)+1.0/(frac*time));
   
   return min;
}

void node::evolve(stellar_evo se,dynamics dyn){
  /*Does what it says in the title  - evolves by rk4 kernel.*/
    
  static double duplicate_nbody[15]; 
  static double dr1[13],dr2[13],dr3[13],dr4[13],dr5[13],dr6[13],errs[13];
  static double safety = 0.9, shrink = -0.25, grow = -0.2;
  double error, step_test;
  
  if (time == 0 && s == 0) tstep = step_min();         //Initial stepsize
  if (time == 0 && s == 1) tstep = 0.1;
  
  if (time+tstep > out_time)     //incase Completing
	tstep = 1.001*(out_time-time);
  
  for (int i=0;i<13;i++) duplicate_nbody[i] = *nbody[i]; //Backup 
  convert(se, dyn);

  for (;;){
    for (int i=0;i<13;i++) *nbody[i] = duplicate_nbody[i];
    solve_odes(dr1,se,dyn);
    convert(se, dyn);
    
    for (int i=0;i<13;i++) 
      *nbody[i]=duplicate_nbody[i]+tstep*(b21*dr1[i]);
    solve_odes(dr2,se,dyn);
    convert(se, dyn); 

    for (int i=0;i<13;i++) 
      *nbody[i]=duplicate_nbody[i]+tstep*(b31*dr1[i]+b32*dr2[i]);
    solve_odes(dr3,se,dyn);
    convert(se, dyn);

    for (int i=0;i<13;i++) 
      *nbody[i]=duplicate_nbody[i]+tstep*(b41*dr1[i]+b42*dr2[i]+b43*dr3[i]);
    solve_odes(dr4,se,dyn);
    convert(se, dyn);
    
    for (int i=0;i<13;i++) 
      *nbody[i]=duplicate_nbody[i]+tstep*(b51*dr1[i]+b52*dr2[i]+b53*dr3[i]+b54*dr4[i]);
    solve_odes(dr5,se,dyn);
    convert(se, dyn);

    for (int i=0;i<13;i++) 
      *nbody[i]=duplicate_nbody[i]+tstep*(b61*dr1[i]+b62*dr2[i]+b63*dr3[i]+b64*dr4[i]+b65*dr5[i]);
    solve_odes(dr6,se,dyn);
    convert(se, dyn);

    for (int i=0;i<13;i++){
      *nbody[i]=duplicate_nbody[i]+tstep*(c1*dr1[i]+c3*dr3[i]+c4*dr4[i]+c6*dr6[i]);
      errs[i] = tstep*(dc1*dr1[i]+dc3*dr3[i]+dc4*dr4[i]+dc5*dr5[i]+dc6*dr6[i]);
    }
    convert(se, dyn);
    
    if (time+tstep > out_time) break;
    
    error = 0;
    for (int i=0;i<13;i++)                   //Checks internal error calculated
       error = fmax(error, fabs(errs[i] / *nbody[i])*1e6); 
    if (error<= 1.0) break;                 //If error < 0.01% accept step
    if (tstep < 1.01*step_min()) break;          //Minimum stepsize - must proceed!
    step_test = safety*tstep*pow(error,shrink);
    tstep=(tstep >= 0.0 ? fmax(step_test,0.1*tstep) : fmin(step_test,0.1*tstep)); 	
  }
  
  if (tstep < step_min()) tstep = step_min();
  else if (error >  1.89e-4) tstep = safety*tstep*pow(error,grow);
  else tstep = 5.0*tstep;
    
 // cerr << Rch << ' ' << dyn.Rchmin() << endl;
  // Switch energy source on at core collapse. 
  if  (Rch <= dyn.Rchmin() && E.source == 0){
      tcc = time;                            //Stores core collapse time
      dyn.reset_K_constants();               //Changes parameters
    }
}

void node::solve_odes(double dvdt[],stellar_evo se,dynamics dyn){
/*Essentially a wrapper function; wraps the solving of all requisite ODES into 
 * one function. Calls to a variety of dynamics module functions.*/	
   
    dvdt[0] = 1;
    dvdt[1] = 0;
    dvdt[2] = dyn.dNdt();
    dvdt[3] = dyn.dmmdt();
    dvdt[4] = se.dmsedt();
    dvdt[5] = dyn.drdt();
    dvdt[6] = dyn.drcdt();
    dvdt[7] = 0;	   
    dvdt[8] = dyn.dtrhdt();
    dvdt[9] = dyn.dtrhpdt();
    dvdt[10] = dyn.dkdt();
    dvdt[11] = se.dmsegdt(); 
    dvdt[12] = dyn.dr2dt();
}

void node::convert(stellar_evo se,dynamics dyn){
 /*Housekeeping function. Calls functions to determine the quantities that need 
  * to be redefined*/
  
    //General functions
    E.value = E_calc(); t_rh = trh(); t_rc = trc();
    psi = se.psi(); t_rhp = trhp(); m_max = se.m_up();
    
    //RV (if used)
    rv = rh/(4.0*kappa);
    
    //updates RG if using dynamical friction
    galaxy.R = sqrt(galaxy.R2);
    
    //Radii calculations, as needed
    rj = r_jacobi(); Rhj = rh/rj; Rch = rc/rh;
    
    //Dynamical Friction Functions
    if (galaxy.f == 1) dyn.tdf(M_star, R_star, pcMyr, T_star);
}
