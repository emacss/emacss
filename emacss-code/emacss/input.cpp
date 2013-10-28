#include "../emacss.h"

void node::direct_input(double in[]){
  /*Gets conditions of N-body run - uses standard C getopt function for 
    command line input. If values are not supplied, defaults are used.*/

    zero();
        
    out_time.Myr = in[0];
    N = in[1];
    r.pc = in[2];
    rhrj = in[3];
    mm.Msun = in[4];
    mm_se.Msun = in[4];
    galaxy.M.Msun = in[5];  //in kpc    
    galaxy.R.pc = in[6]*1e3;  //in kpc
    galaxy.v.kms = in[7];
    E.zeta = in[8];
    
    galaxy.type = 2;
    s = 1;
   
    check_input();
    initialise();
}

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
  
  while (( c = getopt(argc,argv, "N:r:m:t:M:d:v:z:o:s:R:g:k:")) != -1) 
    switch (c){
    case ('N'):                              //N bodies
      value = optarg;
      if (atof(value) < 100 || atof(value) > 1e6){
	  cerr << "N out of valid range (100-1000000)" << endl;
	  exit(1);
      }
      N = atof(value);
      break;
    case ('r'):                              //r at time=0
      value = optarg;
      if (atof(value) < 0 || atof(value) > 25){
	  cerr << "r out of valid range (0-25pc)" << endl;
	  exit(1);
      }
      r.pc = atof(value);
      break;
    case ('R'):                              //Filling Factor
      value = optarg;
      if (atof(value) < 0 || atof(value) > 1){	
	  cerr << "Filling factor not in range 0-1" << endl;
	  exit(1);
      }
      rhrj = atof(value);
      break;
    case ('m'):                              //Mean mass at time=0
      value = optarg;
      if (atof(value) < 0 || atof(value) > 2){	
	  cerr << "Stellar Mass not in range 0-2" << endl;
	  exit(1);
      }
      mm.Msun = atof(value);
      mm_se.Msun = atof(value);  
      break;
    case ('t'):                              //time before stellar evolution (Myr)
      value = optarg;
      T_SE = atof(value);
      break;  
    case ('g'):                              //galaxy halo type flag
      value = optarg;
      if (atoi(value) < 0 || atoi(value) > 2){
	  cerr << "Galaxy Type Invalid" << endl;
	  exit(1);
      }
      galaxy.type = atoi(value);
      break;  
    case ('M'):                              //Galaxy Mass (contained)
      value = optarg;
      if (atof(value) < 0){	
	  cerr << "Galaxy mass negative" << endl;
	  exit(1);
      }
      galaxy.M.Msun = atof(value);
      break;
    case ('d'):                              //distance from galaxy CoM (/1000)
      value = optarg;
      galaxy.R.pc = fabs(atof(value)*1000);  //1000 - input in kpc
      break;  
    case ('v'):                              //Cluster orbital velocity (kms^-2)
      value = optarg;
      if (fabs(atof(value)) > 2000){	
	  cerr << "Orbital velocity high..." << endl;
	  exit(1);
      }
      galaxy.v.kms = fabs(atof(value));
      break;
    case ('o'):                              //Output setting (Unit output)
      value = optarg;
      if (atoi(value) < 0 || atoi(value) > 2){
	  cerr << "Output flag: 0 = N-body, 1 = Real, 2 = Both" << endl;
	  exit(1);
      }
      units = atoi(value);
      break;
    case ('s'):                              //Stellar Evolution on/off
      value = optarg;
      if (atoi(value) < 0 || atoi(value) > 1){
	  cerr << "Stellar Evolution flag: 1 = Real units, 0 = N-body" << endl;
	  exit(1);
      }
      s = atoi(value);
      break; 
    case ('z'):                              //adjustable zeta (mass function)
      value = optarg;
      if (atof(value) < 0 || atof(value) > 1){
	  cerr << "zeta negative or greater than 1" << endl;
	  exit(1);
      }
      E.zeta = atof(value);
      break;   
    case ('k'):                              //Core radius
      value = optarg;
      if (fabs(atof(value)) < 0.3 || fabs(atof(value)) < 0.1){
	  cerr << "Error - energy form factor out of expected range" << endl;
	  exit(1);
      }
      kappa = atoi(value);
      k_0 = atoi(value);
      break; 
    }
    check_input();
    initialise();
}

void node::check_input(){
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
    mm_se.Msun = 0.547; 
  } 

  if (r.pc == 0 && rhrj == 0){      //Check radius of cluster supplied.
    cerr << "No initial radius specified. Assuming r = 1pc" << endl;
    r.pc = 1;    
  }

  if (E.zeta == 0){                //Check energy conduction defined
    E.zeta = 0.10;     
  }
  
  if (out_time.Myr == 0){       //Check output criterion
    cerr << "No output time specified; evolving cluster until dissolution." << endl;
    out_time.Myr = numeric_limits<double>::max();
  }
  
  if (units == 0 && s == 1){       //Check output criterion
    cerr << "Fatal error: N-body units not allowed with stellar evolution" << endl;
    exit(1);
  }
  
  if (s == 0){                    //Sets gamma if s = 0
    gamma = 0.11;
  }

 //Galaxy - if unset, uses solar neighbourhood of MW-like galaxy.
  double G = 4.31572e-3;    //pc M_sun^_1 (km^2)(s^-2)
  if (galaxy.R.pc != 0 && galaxy.v.kms != 0) 
    galaxy.M.Msun = (pow(galaxy.v.kms,2)*galaxy.R.pc)/G;
  else if (galaxy.M.Msun != 0 && galaxy.R.pc == 0 && galaxy.v.kms != 0) 
    galaxy.R.pc = (G*galaxy.M.Msun)/pow(galaxy.v.kms,2);
  else if (galaxy.M.Msun != 0 && galaxy.R.pc != 0 && galaxy.v.kms == 0) 
    galaxy.v.kms = pow((G*galaxy.M.Msun)/galaxy.R.pc,1.0/2.0);
  else if (rhrj != 0 && r.pc != 0 && galaxy.R.pc != 0){
    galaxy.v.kms = pow((G*N*mm.Msun*pow(galaxy.R.pc,2))/3.0*pow(rhrj/r.pc,3),
	    1.0/2.0);
    galaxy.M.Msun = (pow(galaxy.v.kms,2)*galaxy.R.pc)/G;
  }
  else if (galaxy.type == 0){
    cerr << "Isolated cluster assigned" << endl;
    }
  else {
    cerr << "Insufficient galaxy data supplied. Assuming MW-like (solar neighbourhood)" << endl;
    galaxy.v.kms = 220;
    galaxy.R.pc = 8500;
    galaxy.M.Msun = (pow(galaxy.v.kms,2)*galaxy.R.pc)/G;
  }
  
  //Sets radius if not previously defined
  if (r.pc == 0){
      if (galaxy.type == 0){
          cerr << "Isolated galaxy and no initial radius. Set r0=1pc" << endl;
	  r.pc = 1.0;
      }
      else if (galaxy.type == 1)
	  r.pc = rhrj*pow((G*N*mm.Msun)/(3*pow(galaxy.v.kms,2))*pow(galaxy.R.pc,2),
			1.0/3.0);
      else if (galaxy.type == 2)
          r.pc = rhrj*pow((G*N*mm.Msun)/(2*pow(galaxy.v.kms,2))*pow(galaxy.R.pc,2),
			1.0/3.0);
  }
}

void node::initialise(){
    //Sets unit conversions
  G_star = 0.00449857;                    //Grav constant (pc^3M_sun^-1Myr^-2)
  M_star = mm.Msun*N;                     //Initial Mass of cluster (M_sun)
  R_star = r.pc;                          //Virial radii (parsec / N-body)
  T_star = sqrt(pow(R_star,3)/(M_star*G_star)); //N-body time (Myr)
  
  //Converts galaxy to n-body units
  galaxy.M.nbody = galaxy.M.Msun/M_star;
  galaxy.R.nbody = galaxy.R.pc/R_star;

  //Sets final factors needed - masses, E, trh and rj
  time.nbody = time.Myr/T_star;
  mm.nbody = 1.0/N; mm_se.nbody = 1.0/N;
  E.nbody = E_calc(); E.real = E.nbody*pow(M_star,2)/R_star;
  t_relax.nbody = trh(); t_relax.Myr = t_relax.nbody*T_star;
  rj.nbody = r_jacobi(); rj.pc = rj.nbody*R_star;
  trhelapsed = 0; trhp = 0; mass_seg = 1.0;
  out_time.nbody = out_time.Myr/T_star;
  tcc.nbody = tcc.Myr/T_star;
  tin.nbody = tin.Myr/T_star;

  if (time.nbody > tcc.nbody || trhp > T_DYN) coll = 1;
  
  //Sets pointer array for required parameters
  nbody[0] = &time.nbody; nbody[1] = &E.nbody; nbody[2] = &N;
  nbody[3] = &mm.nbody; nbody[4] = &r.nbody; nbody[5] = &rj.nbody;
  nbody[6] = &t_relax.nbody; nbody[7] = &kappa; nbody[8] = &trhelapsed;
  nbody[9] = &trhp; nbody[10] = &mass_seg; nbody[11] = &mm_se.nbody;

  real[0] = &time.Myr; real[1] = &E.real; real[2] = &N;
  real[3] = &mm.Msun; real[4] = &r.pc; real[5] = &rj.pc;
  real[6] = &t_relax.Myr; real[7] = &kappa; real[8] = &trhelapsed;
  real[9] = &trhp; real[10] = &mass_seg; real[11] = &mm_se.Msun;
	    
//  cerr << endl;
}

void node::zero(){
  /*Initialises a cluster to unphysical values. These are default flags used
    by the code to detect which properties are user defined.*/
/*  cerr << endl; 
  cerr << "Evolve Me A Cluster of Stars v2.0" << endl;
  cerr << "-----------------------------------------"<<endl; */
  
  tin.nbody = 0; tin.Myr = 0;
  galaxy.M.Msun = galaxy.R.pc = galaxy.v.kms = 0; galaxy.type = 1;
  time.nbody = time.Myr = out_time.Myr = 0; //Arbitary time to evolve to
  N = 0; mm.nbody = 0; mm.Msun = 0;
  R = 0; r.nbody = 1.0; r.pc = 0; rhrj = 0;
  kappa = 0.2; k_0 = 0.2;
  E.zeta = 0; coll= 0;
  units = 1; s = 1; //Default = stellar evo on, output in Real units
}

