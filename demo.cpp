#include "emacss.h"

int main(){

  //These must be included only once before any clusters are called/evolved

  node cluster;                            //Definition in "emacss.h". Stores cluster properties
  stellar_evo se_module;                   //Module that runs stellar evolution
  dynamics dynamics_module;                //Module than runs dyanmical evolution

  se_module.load(&cluster,&dynamics_module);//Sets up cross module communication
  dynamics_module.load(&cluster,&se_module);//Only needs to be called once 

  double in[13];    //Input array
  double out[14];  //Outpul array

  double time[30000], M[30000], r[30000], dMdt[30000], drdt[30000];      //Storage arrays

  double N0 =  131072, r0 =  11.44548, mm0 = 0.54349, RG = 8.5, vG = 220; //Initial conditions

  //Properties set to 0 are not needed for the first timestep, or have an initial value of 0.
  out[0] = 0; out[1] = N0; out[2] = mm0*N0; out[3] = mm0; out[4] = mm0;
  out[5] = r0; out[6] = 0; out[7] = 0; out[8] = 0; out[9] = 0; out[10] = 0;
  out[11] = 0; out[12] = 0;  out[13] = 1.0; 

  for (int t = 0; t < 30000; t++){

    in[0] = out[0]; in[1] = out[0]+1; in[2] = out[1]; in[3] = out[3];
    in[4] = out[4]; in[5] = out[5]; in[6] = RG; in[7] = vG;
    in[8] = out[11]; in[9] = out[12]; in[10] = out[13]; in[11] = N0*mm0;
    in[12] = r0;

  //These two commands are needed each time a cluster is evolved from time t1 to t2. First call loads the data, second call evolves the cluster and returns the properties at t2.
    cluster.input(in);
    cluster.evaluate(se_module,dynamics_module,out);

    time[t] = out[0]; M[t] = out[2]; r[t] = out[5];
    dMdt[t] = out[9]; drdt[t] = out[10];
  }

  for (int i = 0; i < 30000; i++) cout << time[i] << ' ' << M[i] << ' ' << r[i]\
				       << ' ' << dMdt[i] << ' ' << drdt[i]\
				       << endl;
   	      
  return 0;
}
