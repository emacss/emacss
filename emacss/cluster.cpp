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

double node::E_calc(){                                //Equation (1) AG2012
  return -kappa*pow(N*mm.nbody,2)/r.nbody;
}

double node::E_orb(){       //Orbital energy in log potential
  double r2 = 0, v2 = 0, energy = 0;

    r2 = pow(pos[0],2.0)+pow(pos[1],2.0)+pow(pos[2],2.0); 
    v2 = pow(vel[0],2.0)+pow(vel[1],2.0)+pow(vel[2],2.0);

  energy = 0.5*v2 + pow(galaxy.vc.nbody,2)*log(sqrt(r2));
  
  return energy;
}

double node::r_jacobi(){                            //Equation (18) AG2012
  double rj = 0;
  if (galaxy.type == 1)
	rj = pow((N*mm.nbody)/(3.0*pow(galaxy.v.nbody,2.0)),1.0/3.0)*\
		pow(galaxy.R.nbody,2.0/3.0);
  else if (galaxy.type == 2)
      	rj = pow((N*mm.nbody)/(2.0*pow(galaxy.v.nbody,2.0)),1.0/3.0)*\
		pow(galaxy.R.nbody,2.0/3.0);
  return rj;
}

double node::trh(){                                  //Equation (4) AG2012
  return 0.138*sqrt(N*pow(r.nbody,3)/(mm.nbody))/(log(gamma*N));
}

void node::get_acc(double pos[], double acc[]){    //Orbital acceleration
  double r2 = 0, fmp = 0;
    
  r2 = pow(pos[0],2.0)+pow(pos[1],2.0)+pow(pos[2],2.0);
  fmp = pow(galaxy.vc.nbody,2)/r2;
 
  for (int i = 0; i < 3; i++) acc[i] = -pos[i]*fmp;

 
}

double node::step_min(){                           //Ensures progress!
   return 1.0/(1.0/(frac*t_rhp.nbody)+1.0/(frac*time.nbody));
}

void node::evolve(stellar_evo se,dynamics dyn){
  /*Does what it says in the title  - evolves by rk4 kernel.*/
  static double duplicate_nbody[12]; 
  static double dr1[12],dr2[12],dr3[12],dr4[12],dr5[12],dr6[12],errs[12];
  static double dp1[3],dp2[3],dp3[3],dp4[3],dp5[3],dp6[3],pos_dup[3];
  static double dv1[3],dv2[3],dv3[3],dv4[3],dv5[3],dv6[3],vel_dup[3];
  static double safety = 0.9, shrink = -0.25, grow = -0.2, etest;
  double error, step_test;
  
  convert(se);                              //Checks all quantities calculated

  if (time.nbody == 0) tstep = 0.1;

  if (time.nbody+tstep > out_time.nbody)     //incase Completing
	tstep = 1.001*(out_time.nbody-time.nbody);
  
  for (int i=0;i<12;i++) duplicate_nbody[i] = *nbody[i]; //Backup
  for (int i=0;i<3;i++) pos_dup[i] = pos[i];
  for (int i=0;i<3;i++) vel_dup[i] = vel[i];
  
  for (;;){
    for (int i=0;i<12;i++) *nbody[i] = duplicate_nbody[i]; //Restores backup
    for (int i=0;i<3;i++) pos[i] = pos_dup[i];
    for (int i=0;i<3;i++) vel[i] = vel_dup[i];
	    
    for (int i=0;i<3;i++) dp1[i] = vel_dup[i];                 //Sets orbital int
    get_acc(pos,dv1);
    solve_odes(dr1,se,dyn);
    convert(se);
  
    for (int i=0;i<12;i++) 
      *nbody[i]=duplicate_nbody[i]+tstep*(b21*dr1[i]);
    for (int i=0;i<3;i++){
      pos[i] = pos_dup[i] + tstep*(b21*dp1[i]);
      dp2[i] = vel[i] = vel_dup[i] + tstep*(b21*dv1[i]);	
    }
    get_acc(pos,dv2);
    solve_odes(dr2,se,dyn);
    convert(se);
   
    for (int i=0;i<12;i++) 
      *nbody[i]=duplicate_nbody[i]+tstep*(b31*dr1[i]+b32*dr2[i]);
    for (int i=0;i<3;i++){
      pos[i] = pos_dup[i] + tstep*(b31*dp1[i]+b32*dp2[i]);
      dp3[i] = vel[i] = vel_dup[i] + tstep*(b31*dv1[i]+b32*dv2[i]);
    }
    get_acc(pos,dv3);
    solve_odes(dr3,se,dyn);
    convert(se);
  
    for (int i=0;i<12;i++) 
      *nbody[i]=duplicate_nbody[i]+tstep*(b41*dr1[i]+b42*dr2[i]+b43*dr3[i]);
    for (int i=0;i<3;i++){
      pos[i] = pos_dup[i] + tstep*(b41*dp1[i]+b42*dp2[i]+b43*dp3[i]);
      dp4[i] = vel[i] = vel_dup[i] + tstep*(b41*dv1[i]+b42*dv2[i]+b43*dv3[i]);
    }
    get_acc(pos,dv4);
    solve_odes(dr4,se,dyn);
    convert(se);
  
    for (int i=0;i<12;i++) 
      *nbody[i]=duplicate_nbody[i]+tstep*(b51*dr1[i]+b52*dr2[i]+b53*dr3[i]+b54*dr4[i]);
    for (int i=0;i<3;i++){
      pos[i] = pos_dup[i] + tstep*(b51*dp1[i]+b52*dp2[i]+b53*dp3[i]+b54*dp4[i]);
      dp5[i] = vel[i] = vel_dup[i] + tstep*(b51*dv1[i]+b52*dv2[i]+b53*dv3[i]+b54*dv4[i]);
    }
    get_acc(pos,dv5);
    solve_odes(dr5,se,dyn);
    convert(se);
  
    for (int i=0;i<12;i++) 
      *nbody[i]=duplicate_nbody[i]+tstep*(b61*dr1[i]+b62*dr2[i]+b63*dr3[i]+b64*dr4[i]+b65*dr5[i]);
    for (int i=0;i<3;i++){
      pos[i] = pos_dup[i] + tstep*(b61*dp1[i]+b62*dp2[i]+b63*dp3[i]+b64*dp4[i]+b65*dp5[i]);
      dp6[i] = vel[i] = vel_dup[i] + tstep*(b61*dv1[i]+b62*dv2[i]+b63*dv3[i]+b64*dv4[i]+b65*dv5[i]);
    }
    get_acc(pos,dv6);
    solve_odes(dr6,se,dyn);     
    convert(se); 
  
    for (int i=0;i<12;i++){
      *nbody[i]=duplicate_nbody[i]+tstep*(c1*dr1[i]+c3*dr3[i]+c4*dr4[i]+c6*dr6[i]);
      errs[i] = tstep*(dc1*dr1[i]+dc3*dr3[i]+dc4*dr4[i]+dc5*dr5[i]+dc6*dr6[i]);
    }
    for (int i=0;i<3;i++){
      pos[i] = pos_dup[i]+tstep*(c1*dp1[i]+c3*dp3[i]+c4*dp4[i]+c6*dp6[i]);
      vel[i] = vel_dup[i]+tstep*(c1*dv1[i]+c3*dv3[i]+c4*dv4[i]+c6*dv6[i]);
     }
    convert(se); 
    
   if (time.nbody+tstep > out_time.nbody) break;
   error = 0;
   for (int i=0;i<12;i++)                   //Checks internal error calculated
       error = fmax(error, fabs(errs[i] / *nbody[i])*100); 
   etest = E_orb();
   error = fmax(fabs(((E_orbital-etest) / E_orbital)*1e9),error);
      
   if (error<= 1.0) break;                 //If error < 1% accept step
   if (tstep < step_min()) break;          //Minimum stepsize - must proceed!
   step_test = safety*tstep*pow(error,shrink);
   tstep=(tstep >= 0.0 ? fmax(step_test,0.1*tstep) : fmin(step_test,0.1*tstep)); 	
   
  }
  if (tstep > step_min()){
    if (error >  1.89e-4) tstep = safety*tstep*pow(error,grow);
    else tstep = 5.0*tstep;
  } 
 E_orbital = etest;
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
    dvdt[7] = dyn.dkdt();
    dvdt[8] = dyn.dtrhdt();
    dvdt[9] = dyn.dtrhpdt();
    dvdt[10] = se.dmsegdt();
    dvdt[11] = se.dmsedt();
}

void node::convert(stellar_evo se){
 /*Housekeeping function. Keeps n-body and real units aligned when called.
  Also calls functions to determine the quantities that need to be redefined*/
  
    E.nbody = E_calc(); E.real = E.nbody*pow(M_star,2)/R_star;
    t_relax.nbody = trh(); t_relax.Myr = t_relax.nbody*T_star;
    t_rhp.nbody = t_relax.nbody/se.zeta(); t_rhp.Myr = t_rhp.nbody*T_star;
    galaxy.R.nbody = sqrt(pow(pos[0],2)+pow(pos[1],2)+pow(pos[2],2));
    galaxy.v.nbody = sqrt(pow(vel[0],2)+pow(vel[1],2)+pow(vel[2],2));
    galaxy.R.pc = galaxy.R.nbody*R_star; 
    galaxy.v.kms = galaxy.v.nbody*V_star;   
    rj.nbody = r_jacobi(); rj.pc = rj.nbody*R_star;
    rhrj = r.nbody/rj.nbody;
       
    *real[0] = *nbody[0]*T_star;
    *real[2] = *nbody[2];
    *real[3] = *nbody[3]*M_star;
    *real[4] = *nbody[4]*R_star;
    *real[7] = *nbody[7];
    *real[8] = *nbody[8];
    *real[9] = *nbody[9];
    *real[10] = *nbody[10];
    *real[11] = *nbody[11]*M_star;
}
