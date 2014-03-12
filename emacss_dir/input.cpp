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
  
  while (( c = getopt(argc,argv, "N:r:m:t:M:d:v:z:o:s:R:g:l:u:f:")) != -1) 
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
      rh = atof(value);
      break;
    case ('R'):                              //Filling Factor
      value = optarg;
      if (atof(value) < 0 || atof(value) > 0.3){	
	  cerr << "Filling factor not in range 0-0.3" << endl;
	  exit(1);
      }
      Rhj = atof(value);
      break;
    case ('m'):                              //Mean mass at time=0
      value = optarg;
      if (atof(value) < 0 || atof(value) > 2){	
	  cerr << "Stellar Mass not in range 0-2" << endl;
	  exit(1);
      }
      mm = atof(value);
      mm_se = atof(value);
      break;
    case ('t'):                              //End time (if wanted)
      value = optarg;
      if (atof(value) < 0 ){	
	  cerr << "Time interval for evolution must be positive" << endl;
	  exit(1);
      }
      out_time = atof(value);
      break;
    case ('u'):                              //Maximum of mass function
      value = optarg;
      if (atof(value) < 1 || atof(value) > 100){	
	  cerr << "Maximum stellar mass not in range 1-100" << endl;
	  exit(1);
      }
      m_max = atof(value);
      break;  
    case ('l'):                              //Minimum of mass function
      value = optarg;
      if (atof(value) < 0 || atof(value) > 1){
	  cerr << "Minimum stellar mass not in range 0-1" << endl;
	  exit(1);  
      }
      m_min = atof(value);
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
      galaxy.M = atof(value);
      break;
    case ('d'):                              //distance from galaxy CoM (/1000)
      value = optarg;
      galaxy.R = fabs(atof(value)*1000);     //1000 - input in kpc
      break;  
    case ('v'):                              //Cluster orbital velocity (kms^-2)
      value = optarg;
      if (fabs(atof(value)) > 2000){	
	  cerr << "Orbital velocity high, aborting..." << endl;
	  exit(1);
      }
      galaxy.v = fabs(atof(value));
      break;
    case ('o'):                              //Output setting (Unit output)
      value = optarg;
      if (atoi(value) < 0 || atoi(value) > 1){
	  cerr << "Output flag: 0 = N-body, 1 = Real" << endl;
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
    case ('f'):                              //Dynamical Friction on/off
      value = optarg;
      if (atoi(value) < 0 || atoi(value) > 1){
	  cerr << "Dynamical friction flag: 1 = On, 0 = Off" << endl;
	  exit(1);
      }
      galaxy.f = atoi(value);
      break;
    case ('z'):                              //adjustable zeta (mass function)
      value = optarg;
      if (atof(value) < 0 || atof(value) > 1){
	  cerr << "zeta negative or greater than 1" << endl;
	  exit(1);
      }
      E.zeta = atof(value);
      break;   
    }
}

void node::initialise(stellar_evo se,dynamics dyn){
  /*Initialises all undefined values - various sections serve various purposes,
    as described.*/
  double G = 1;                  //Gravitational constant (N-body input)
  if (units > 0) G = 4.31572e-3; //pc M_sun^_1 (km^2)(s^-2) (recast)
  
  if (units == 0 && s == 1){                 //Checks for valid input
    cerr << "Fatal error: N-body units not allowed with stellar evolution" << endl;
    exit(1);
  } 
   
  if (out_time == 0){                        //Check output criterion
    cerr << "No output time specified; evolving cluster until dissolution." << endl;
    out_time = numeric_limits<double>::max();
  }
   
  //Cluster - checks sufficient values are entered (general)
  if (N == 0){                           //Check number of stars supplied.
    cerr << "No number of stars specified. Assuming N=64k" << endl;
    N = 65536;
  }
     
  if (E.zeta == 0){                //Check energy conduction defined
    cerr << "No zeta specified. Assuming default (zeta = 0.1)" << endl;   
    E.zeta = 0.1;
  }
  
  if (mm == 0){      
    cerr << "No initial mean mass specified. Assuming mm = 0.64M_sun" << endl;
    if (units == 1) mm = 0.64; 
    else mm = 1.0/N;
  }
  
  if (rh== 0){
      if (units == 0) rh = 0.78;
      else rh = 1.0;
  }
         
 //Galaxy set up - if unset, uses isolated.
  if (galaxy.M == 0 && galaxy.R != 0 && galaxy.v != 0){ 
    galaxy.M = (pow(galaxy.v,2)*galaxy.R)/G;
    if (galaxy.type == 0) galaxy.type = 2;
    cerr << "Galaxy set by velocity and radius" << endl;
  }
  else if (galaxy.M != 0 && galaxy.R == 0 && galaxy.v != 0){
    galaxy.R = (G*galaxy.M)/pow(galaxy.v,2);
    if (galaxy.type == 0) galaxy.type = 2;
    cerr << "Galaxy set by radius and mass" << galaxy.R << endl;
  }
  else if (galaxy.M != 0 && galaxy.R != 0 && galaxy.v == 0){
    galaxy.v = pow((G*galaxy.M)/galaxy.R,1.0/2.0);
    if (galaxy.type == 0) galaxy.type = 2;
    cerr << "Galaxy set by mass and radius" << endl;
  }
  else if (Rhj != 0){
    cerr << "Galaxy set by rh/rj" << endl;
    if (galaxy.M == 0) galaxy.M = 1e10; //Assumption - makes it work
    if (galaxy.type == 1)
      galaxy.R = (rh/Rhj)*pow((3.0*galaxy.M)/(N*mm),1.0/3.0);
    else{
      galaxy.R = (rh/Rhj)*pow((2.0*galaxy.M)/(N*mm),1.0/3.0);
      galaxy.type = 2;
    }
    galaxy.v = pow((G*galaxy.M)/galaxy.R,1.0/2.0);
  }
  else{
    cerr << "Insufficient galaxy data supplied. Assuming isolated..." << endl;
    galaxy.type = 0;
  }
 
  
  if (s == 0 && galaxy.type == 0){                 //Checks for valid input
    cerr << "Fatal error: cannot handle isolated cluster without stellar evolution." << endl;
    exit(1);
  } 
   
  //Checks enough data for real output units
  if (units == 1){            
    if (rh == 0){   
      if (galaxy.type == 0){
        cerr << "Isolated galaxy and no initial radius. Set r0=1pc" << endl;
        rh = 1.0;
      }
      else if (galaxy.type == 1)
        rh = Rhj*pow((G*N*mm)/(3*pow(galaxy.v,2))*pow(galaxy.R,2),1.0/3.0);
      else if (galaxy.type == 2)
        rh = Rhj*pow((G*N*mm)/(2*pow(galaxy.v,2))*pow(galaxy.R,2),1.0/3.0); 
    }
    
  //Sets unit conversions - Assumes conversion to Plummer sphere
    pcMyr = 0.977813106 ;                   //Converts km/s to pc/Myr
    G_star = 0.00449857;                    //Grav constant (pc^3M_sun^-1Myr^-2)
    M_star = mm*N;                          //Initial Mass of cluster (M_sun)
    R_star = rh/0.78;                       //Virial radii (parsec / N-body)
    T_star = sqrt(pow(R_star,3)/(M_star*G_star)); //N-body time (Myr)

  //If needed, converts galaxy to N-body units 
    galaxy.M = galaxy.M/M_star;
    galaxy.R = galaxy.R/R_star;
    galaxy.v = (galaxy.v*pcMyr)/(R_star/T_star); 
  }
  
  //Sets up RG^2 for dynamical friction 
  galaxy.R2 = pow(galaxy.R,2);
  
  //Sets Coulomb logarithm
  if (s == 0) gamma = 0.11;
  else gamma = 0.02;
   
  //Sets final factors needed - kappa, mass, E, trh , trc, and rj
  mm = 1.0/N; mm_se = mm; rv = 1.0; rh = 0.78; psi = se.psi();
  t_rh = trh(); t_rc = trc(); t_rhp = trhp();
  rj = r_jacobi(); Rhj = rh/rj; Rch = rc/rh;
  trhelapsed = 0; trhpelapsed = 0; MS = 3.0; out_time = out_time/T_star;
  m_min = m_min/M_star; m_max = m_max0 = m_max/M_star; 
  E.value = E_calc();
  
  if (s == 0) kappa = k0 = rh/(4.0*rv);
  else kappa = 0.2;
  
  if (galaxy.type != 0 && Rhj > 0.5){             //Extra check of ratio Rhj
    cerr << Rhj << endl;     
    cerr << "Error: rh/rj > 0.5 (out of acceptable range). Exiting..." << endl;
    exit(10);
  }
  
  //Sets pointer array for required parameters
  nbody[0] = &time; 
  nbody[1] = &E.value; 
  nbody[2] = &N;
  nbody[3] = &mm; 
  nbody[4] = &mm_se; 
  nbody[5] = &rh; 
  nbody[6] = &rc;
  nbody[7] = &rj;
  nbody[8] = &trhelapsed;
  nbody[9] = &trhpelapsed;
  nbody[10] = &kappa;
  nbody[11] = &MS;
  nbody[12] = &galaxy.R2;
  
  cerr << endl;
}

void node::zero(){
  /*Initialises a cluster to unphysical values. These are default flags used
    by the code to detect which properties are user defined.*/
  cerr << endl; 
  cerr << "Evolve Me A Cluster of Stars v3.10" << endl;
  cerr << "-----------------------------------------"<<endl; 
  
  //General Properties
  N = 0; rh = 0; rv = 1.0; rc = 0.32; Rhj = 0;
  
  //Masses
  mm = mm_se = 0; m_min = 0.1; m_max = 100; 
  
  //Normalisation
  G_star = M_star = R_star = T_star = 1;
  
  //Galaxy
  galaxy.M = galaxy.R = galaxy.v = 0; galaxy.type = 1, galaxy.f = 0;
  
  //Times: out time = termination time, t_ref = reference time.
  time = 0; out_time = 0;
  
  //Form factors
  kappa = k0 = 0; psi = 0;
  
  E.zeta = 0; E.source = 0;
  
  units = 1; s = 1; //Default = stellar evo on, output in Real units
}
