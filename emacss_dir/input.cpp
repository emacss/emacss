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
  
  while (( c = getopt(argc,argv, "N:r:R:m:g:M:d:v:o:z:")) != -1) 
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
    case ('z'):                              //adjustable zeta (mass function)
      value = optarg;
      if (atof(value) < 0 || atof(value) > 1){
	  cerr << "zeta negative or greater than 1" << endl;
	  exit(1);
      }
      E.zeta = atof(value);
      break;   
    }
    initialise();
}

void node::initialise(){
  /*Initialises all undefined values - various sections serve various purposes,
    as described.*/
   double G = 1;                //Gravitational constant (N-body input) 
   
  //Cluster - checks sufficient values are entered (general)
  if (N == 0){                           //Check number of stars supplied.
    cerr << "No number of stars specified. Assuming N=64k" << endl;
    N = 65536;
  }
  
  if (E.zeta == 0){                //Check energy conduction defined
    cerr << "No zeta specified. Assuming default (zeta = 0.1)" << endl;   
    E.zeta = 0.1;
  }
   
  if (rh != 0.78 && units == 0){       //Check radius
    cerr << "Half mass radius entered but N-body units. Ignoring entry" << endl;   
  }
   
  if (Rhj == 0 && galaxy.M == 0 && galaxy.R == 0 && galaxy.v == 0){//Check gal
    cerr << "No tidal field specified - assuming rh/rj = 0.1" << endl;   
    Rhj = 0.1;
  }
  //Cluster - checks sufficient values are entered (real output units)
  if (units > 0){
    if (mm == 0 && units > 0){                //Check mass of cluster supplied.
      cerr << "No initial mean mass specified. Assuming mm = 0.5M_sun" << endl;
      mm = 0.5;  
    } 

    if (rh == 0.78){      //Check radius of cluster supplied.
      cerr << "No initial radius specified. Assuming r = 0.78pc" << endl;
      rh = 0.78;    
    }
    
  //Sets unit conversions - Assumes conversion to Plummer sphere
    pcMyr = 0.977813106 ;                   //Converts km/s to pc/Myr
    G_star = 0.00449857;                    //Grav constant (pc^3M_sun^-1Myr^-2)
    M_star = mm*N;                          //Initial Mass of cluster (M_sun)

    R_star = rh/0.78;                       // (parsec / N-body)
    T_star = sqrt(pow(R_star,3)/(M_star*G_star)); //N-body time (Myr)
  
  //To check galaxy conditions
    G = 4.31572e-3;    //pc M_sun^_1 (km^2)(s^-2), for galaxy setup
  }
   
 //Sets final factors needed - kappa, mass, E, trh , trc, and rj
 if (units < 1){
   mm = 1.0/N; rh = 0.78; kappa = rh/(4.0*rv);
   t_rh = trh(); t_rc = trc(); E.value = E_calc(); trhelapsed = 0; 
 }

 //Galaxy set up - if unset, uses solar neighbourhood of MW-like galaxy.
  if (galaxy.M == 0 && galaxy.R != 0 && galaxy.v != 0){ 
    galaxy.M = (pow(galaxy.v,2)*galaxy.R)/G;
    cerr << "Galaxy set by mass and radius" << endl;
  }
  else if (galaxy.M != 0 && galaxy.R == 0 && galaxy.v != 0){
    galaxy.R = (G*galaxy.M)/pow(galaxy.v,2);
    cerr << "Galaxy set by radius and velocity" << galaxy.R << endl;
  }
  else if (galaxy.M != 0 && galaxy.R != 0 && galaxy.v == 0){
    galaxy.v = pow((G*galaxy.M)/galaxy.R,1.0/2.0);
    cerr << "Galaxy set by mass and velocity" << endl;
  }
  else if (Rhj != 0){
    cerr << "Galaxy set by rh/rj" << endl;
    if (galaxy.M == 0) galaxy.M = 1e10; //Assumption - makes it work
    if (galaxy.type == 1)
      galaxy.R = (rh/Rhj)*pow((3.0*galaxy.M)/(N*mm),1.0/3.0);
    else if (galaxy.type == 2)
      galaxy.R = (rh/Rhj)*pow((2.0*galaxy.M)/(N*mm),1.0/3.0);
    galaxy.v = pow((G*galaxy.M)/galaxy.R,1.0/2.0);
  }
  else{
    cerr << "Insufficient galaxy data supplied. Assuming isolated..." << endl;
    galaxy.type = 0;
  }
  
  //If needed, converts galaxy to N-body units
  if (units > 0){ 
    galaxy.M = galaxy.M/M_star;
    galaxy.R = galaxy.R/R_star;
    galaxy.v = (galaxy.v*pcMyr)/(R_star/T_star); 
    mm = 1.0/N; rh = 0.78; kappa = rh/(4.0*rv);
    t_rh = trh(); t_rc = trc(); E.value = E_calc(); trhelapsed = 0; 
  }
 
  //A few more factors set (that need the galaxy conditions)
  rj = r_jacobi(); Rhj = rh/rj; Rch = rc/rh;
  if (Rhj > 0.33){
    cerr << "Error: rh/rj out of acceptable range. Exiting..." << endl;
    exit(10);
  }

  //Sets pointer array for required parameters
  nbody[0] = &time; 
  nbody[1] = &E.value; 
  nbody[2] = &N;
  nbody[3] = &mm; 
  nbody[4] = &rh; 
  nbody[5] = &rj;
  nbody[6] = &t_rh; 
  nbody[7] = &kappa;
  nbody[8] = &trhelapsed;
  nbody[9] = &rc;
  nbody[10] = &t_rc;
  
  cerr << endl;
}

void node::zero(){
  /*Initialises a cluster to unphysical values. These are default flags used
    by the code to detect which properties are user defined.*/
  cerr << endl; 
  cerr << "Evolve Me A Cluster of Stars v2.03" << endl;
  cerr << "-----------------------------------------"<<endl; 
  
  galaxy.M = galaxy.R = galaxy.v = 0; galaxy.type = 1;
  G_star = M_star = R_star = T_star = 1;
  time = 0;                       //Out time = termination time.
  N = 0;  mm = 0; mm = 0; Rhj = 0; 

  rh = 0.78; rv = 1.0; rc = 0.32; Rhj = 0;  
  kappa = E.zeta = E.source = 0;
  
  units = 1; //Default = output in Real units
}
