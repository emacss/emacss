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

class emacss{
    node cluster;
    stellar_evo se_module;
    dynamics dynamics_module; 
    
public:
    
  emacss(){
    cluster.zero();
    se_module.load(&cluster,&dynamics_module);
    dynamics_module.load(&cluster,&se_module);
   }
     
  void input(double in[]){

    cluster.zero();
   
    cluster.time.Myr = in[0];
    cluster.tin.Myr = in[0];
    cluster.out_time.Myr = in[1];
    cluster.N = in[2];
    cluster.r.pc = in[3];
    cluster.mm.Msun = in[4];
    cluster.mm_se.Msun = in[4];
    cluster.galaxy.R.pc = in[5]*1e3;  //in kpc
    cluster.galaxy.v.kms = in[6];
    cluster.tcc.Myr= in[7];
    cluster.trhp = in[8];
    
    cluster.galaxy.type = 2;
    cluster.s = 1;
      
    cluster.check_input();
    cluster.initialise();      
  }
    
  void evaluate(double out[]){ 
      
    while (cluster.time.Myr < cluster.out_time.Myr-0.5 && cluster.N > 200)
       cluster.evolve(se_module,dynamics_module);

    if (cluster.N < 200 || cluster.N != cluster.N){
	cluster.N = 0;
	cluster.r.pc = 0;
	cluster.mm.Msun = 0;
    }
    out[0] = cluster.time.Myr;
    out[1] = cluster.N;
    out[2] = cluster.mm.Msun*cluster.N;
    out[3] = cluster.r.pc;
    out[4] = cluster.rj.pc;
    out[5] = dynamics_module.dNdt()/cluster.T_star;
    out[6] = dynamics_module.dmmdt()*(cluster.M_star/cluster.T_star);   
    out[7] = dynamics_module.dmmdt()*(cluster.M_star/cluster.T_star)*\
        cluster.N+dynamics_module.dNdt()/cluster.T_star*cluster.mm.Msun;   
    out[8] = dynamics_module.drdt()*(cluster.R_star/cluster.T_star); 
    out[9] = cluster.tcc.Myr;
    out[10] = cluster.trhp;
  }
};

extern "C" {
    emacss* emacss_new(){ return new emacss(); }
    double in[10], out[11];
    void input(emacss* em, double in[]){ em->input(in);}
    void evaluate(emacss* em, double out[]){ em->evaluate(out); }
}
