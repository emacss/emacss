/***********************************************************************/
//  Evolve Me A Cluster of StarS (EMACSS)                              //
//                                                                     //
//  By: Poul Alexander (University of Cambridge)                       //
//      Mark Gieles (University of Cambridge, University of Surrey)    //
//      Henny Lamers (University of Amsterdam)                         //
//      Holger Baumgardt (Univesity of Brisbane)                       //
//                                                                     //
//  EMACSS is a numerical integrator that solves for the differential  //
//  equations:                                                         //
//
//  dE/dr = e(N,r,p) = -epsilon/t_rh
//                                                                     //
//  dN/dt = f(N,r,p) = - xi*N/t_rh                                     // 
//                                                                     //
//  dr/dt = g(N,r,p) = (zeta - 2*xi) * r/t_rh                          //
//                                                                     // 
//  where                                                              // 
//                                                                     // 
//  epsilon = f(N,r,rt,t)                                              //
//  xi = f(N,r,rt)                                                     //
//  p is an array of parameters/variables, which define the properties //
//  of star clusters at any given time in the life-cycle               //
//                                                                     //
// Equations referenced AG2012 refer to Alexander & Gieles 2012        //
// Equations referenced AGLB2013 refer to Alexander et al. 2013        //
/***********************************************************************/

#include "emacss.h"

void node::input(double in[]){

  zero();
  
  time.Myr = in[0];
  tin.Myr = in[0];
  out_time.Myr = in[1];
  N = in[2];
  mm.Msun = in[3];
  mm_se.Msun = in[4];
  r.pc = in[5];
  galaxy.R.pc = in[6]*1e3;  //in kpc
  galaxy.v.kms = in[7];
  coll = in[8];
  trhp = in[9];
  mass_seg = in[10];
  M0.Msun = in[11];
  r0.pc = in[12];
  
  galaxy.type = 2;
  s = 1;
 
  check_input();
  initialise();      
}

void node::evaluate(stellar_evo se_module, dynamics dynamics_module, double out[]){ 

  while (time.Myr < out_time.Myr && N > 200)
    evolve(se_module,dynamics_module);
   
  if (N < 200 || N != N){
    N = 0;
    r.pc = 0;
    mm.Msun = 0;
  }

  out[0] = time.Myr;
  out[1] = N;
  out[2] = mm.Msun*N;
  out[3] = mm.Msun;
  out[4] = mm_se.Msun;
  out[5] = r.pc;
  out[6] = rj.pc;
  out[7] = dynamics_module.dNdt()/T_star;
  out[8] = dynamics_module.dmmdt()*(M_star/T_star);   
  out[9] = dynamics_module.dmmdt()*(M_star/T_star)*	\
    N+dynamics_module.dNdt()/T_star*mm.Msun;   
  out[10] = dynamics_module.drdt()*(R_star/T_star); 
  out[11] = coll;
  out[12] = trhp;
  out[13] = mass_seg;

}


