#include "../emacss.h"

void node::input(int argc, char* argv[]){
  /*Gets conditions of N-body run - uses standard C getopt function for 
    command line input. If values are not supplied, defaults are used.*/

  string check;
  int c;
  char *value = NULL;
  
  if (argc > 1) { 
    check = argv[1];
    if (check == "-help" || check == "-h") help();
    if (check == "-version") version();
  }
    
  zero();  //Zero's the code (i.e., sets a default set of conditions) 
  
  while (( c = getopt(argc,argv, "N:r:m:t:M:d:v:z:o:s:R:g:")) != -1) 
    switch (c){
    case ('N'):                              //N bodies
      value = optarg;
      N = atof(value);
      break;
    case ('r'):                              //r at time=0
      value = optarg;
      r.pc = atof(value);
      break;
    case ('R'):                              //Filling Factor
      value = optarg;
      rhrj = atof(value);
      break;
    case ('m'):                              //Mean mass at time=0
      value = optarg;
      mm.Msun = atof(value);
      break;
    case ('t'):                              //output time
      value = optarg;
      out_time.Myr = atof(value);
      break;  
    case ('g'):                              //galaxy halo type flag
      value = optarg;
      galaxy.type = atoi(value);
      break;  
    case ('M'):                              //Galaxy Mass (contained)
      value = optarg;
      galaxy.M.Msun = atof(value);
      break;
    case ('d'):                              //distance from galaxy CoM (/1000)
      value = optarg;
      galaxy.R.pc = atof(value)*1000;        //1000 - input in kpc
      break;  
    case ('v'):                              //Cluster orbital velocity (kms^-2)
      value = optarg;
      galaxy.v.kms = atof(value);
      break;
    case ('o'):                              //Output setting (Unit output)
      value = optarg;
      units = atoi(value);
      break;
    case ('s'):                              //Stellar Evolution on/off
      value = optarg;
      s = atoi(value);
      break; 
    case ('z'):                              //adjustable zeta (mass function)
      value = optarg;
      E.zeta = atof(value);
      break;   
    }
    initialise();
}

void node::initialise(){
  /*Initialises all undefined values - various sections serve various purposes,
    as described.*/
  
  //Cluster - checks sufficient values are entered. 
  if (N == 0){                           //Check number of stars supplied.
    cerr << "No number of stars specified. Assuming N=100000" << endl;
    N = 100000;
  }

  if (mm.Msun == 0){                //Check mass of cluster supplied.
    cerr << "No initial mean mass specified. Assuming mm = 0.547M_sun" << endl;
    mm.Msun = 0.547;  
  } 

  if (r.pc == 0 && rhrj == 0){      //Check radius of cluster supplied.
    cerr << "No initial radius specified. Assuming r = 1pc" << endl;
    r.pc = 1;    
  }

  if (E.zeta == 0){                //Check energy conduction defined
    cerr << "No zeta specified. Assuming mass spectrum (zeta = 0.15)" << endl;
    E.zeta = 0.15;     
  }
  
  if (out_time.Myr == 0){       //Check output criterion
    cerr << "No output time specified; evolving cluster until dissolution." << endl;
    out_time.Myr = numeric_limits<double>::max();
  }
  
  if (units == 0 && s == 1){       //Check output criterion
    cerr << "Fatal error: N-body units not allowed with stellar evolution" << endl;
    exit(1);
  }

 //Galaxy - if unset, uses solar neighbourhood of MW-like galaxy.
  double G = 4.31572e-3;    //pc M_sun^_1 (km^2)(s^-2)
  if (galaxy.M.Msun == 0 && galaxy.R.pc != 0 && galaxy.v.kms != 0) 
    galaxy.M.Msun = (pow(galaxy.v.kms,2)*galaxy.R.pc)/G;
  else if (galaxy.M.Msun != 0 && galaxy.R.pc == 0 && galaxy.v.kms != 0) 
    galaxy.R.pc = (G*galaxy.M.Msun)/pow(galaxy.v.kms,2);
  else if (galaxy.M.Msun != 0 && galaxy.R.pc != 0 && galaxy.v.kms == 0) 
    galaxy.v.kms = pow((G*galaxy.M.Msun)/galaxy.R.pc,1.0/2.0);
  else{
    cerr << "Insufficient galaxy data supplied. Assuming MW-like" << endl;
    galaxy.v.kms = 220;
//    galaxy.R.pc = 8500;
    galaxy.M.Msun = (pow(galaxy.v.kms,2)*galaxy.R.pc)/G;
  }
  
  //Sets radius if not previously defined
  if (r.pc == 0){
	if (galaxy.type == 0)
	   r.pc = rhrj*pow((G*N*mm.Msun)/(3*pow(galaxy.v.kms,2))*pow(galaxy.R.pc,2),
		1.0/3.0);
	else if (galaxy.type == 1)
           r.pc = rhrj*pow((G*N*mm.Msun)/(2*pow(galaxy.v.kms,2))*pow(galaxy.R.pc,2),
		1.0/3.0);
  }
  
    //Sets unit conversions
  G_star = 0.00449857;                    //Grav constant (pc^3M_sun^-1Myr^-2)
  M_star = mm.Msun*N;                     //Initial Mass of cluster (M_sun)
  R_star = r.pc;                          //Virial radii (parsec / N-body)
  T_star = sqrt(pow(R_star,3)/(M_star*G_star)); //N-body time (Myr)
  
  //Converts galaxy to n-body units
  galaxy.M.nbody = galaxy.M.Msun/M_star;
  galaxy.R.nbody = galaxy.R.pc/R_star;

  //Sets final factors needed - masses, E, trh and rj
  mm.nbody = 1.0/N; DM_SE.Msun = mm.Msun*N;
  E.nbody = E_calc(); E.real = E.nbody*pow(M_star,2)/R_star;
  t_relax.nbody = trh(); t_relax.Myr = t_relax.nbody*T_star;
  rj.nbody = r_jacobi(); rj.pc = rj.nbody*R_star; trhelapsed = 0;
  out_time.nbody = out_time.Myr/T_star;
  
  //Sets pointer array for required parameters
  nbody[0] = &time.nbody; nbody[1] = &E.nbody; nbody[2] = &N;
  nbody[3] = &mm.nbody; nbody[4] = &r.nbody; nbody[5] = &rj.nbody;
  nbody[6] = &t_relax.nbody; nbody[7] = &DM_SE.nbody; nbody[8] = &kappa;
  nbody[9] = &trhelapsed;

  real[0] = &time.Myr; real[1] = &E.real; real[2] = &N;
  real[3] = &mm.Msun; real[4] = &r.pc; real[5] = &rj.pc;
  real[6] = &t_relax.Myr; real[7] = &DM_SE.Msun; real[8] = &kappa;
  real[9] = &trhelapsed;
  
  cerr << endl;
}

void node::zero(){
  /*Initialises a cluster to unphysical values. These are default flags used
    by the code to detect which properties are user defined.*/
  cerr << endl; 
  cerr << "Evolve Me A Cluster of Stars v2.0" << endl;
  cerr << "-----------------------------------------"<<endl; 
  
  galaxy.M.Msun = galaxy.R.pc = galaxy.v.kms = galaxy.type = 0;
  time.nbody = time.Myr = out_time.Myr = 0; //Arbitary time to evolve to
  N = 0; DM_SE.nbody = 1;
  mm.nbody = 0; mm.Msun = 0;
  R = 0; r.nbody = 1; r.pc = 0; rhrj = 0;
  kappa = 0.25;
  E.zeta = 0;
  units = 1; s = 1; //Default = stellar evo on, output in Real units
}
