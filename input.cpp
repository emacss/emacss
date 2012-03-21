#include "emacss.h"

/* Functions that either accept input from the command line or use inbuilt
   prescriptions to describe a cluster's evolution prior to the core collapse.

   Essential inputs: N, rh0 & rhrj (or galaxy R and M)
   Optional inputs: pre-core collapse expansion (f_r), mass loss (f_N)*/

void readinput(variables *init, int argc, char* argv[]){
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
    
  zero(init); //Zero's the code (i.e., sets a default set of conditions)
  while (( c = getopt 
	   (argc,argv, "u:N:r:m:k:t:R:M:d:v:l:x:e:a:z:c:")) != -1)
    switch (c){
    case ('u'):                              //User unit selection (for input)
      value = optarg;
      init->units = atoi(value);
      break;  
    case ('N'):                              //N bodies
      value = optarg;
      init-> N0 = atof(value);
	  if (init->N0 < 256) {
		cerr << "Error - N to low. Aborting..." << endl;
		exit(1);
	  }
	  break;
    case ('r'):                              //r at time=0
      value = optarg;
      init->r0 = atof(value);
      break;
      //Note that output is in terms of r (virial radius), where r = r_h/0.77
      //based on the assumption of a plummer sphere.
    case ('m'):                              //mm at time=0
      value = optarg;
      init->mm0 = atof(value);
      break;
    case ('t'):                              //tcc (20trh0 in AG2012)
      value = optarg;
      init->tcc = atof(value);
      break;
    case ('R'):                              //r_h/r_j at time=0
      value = optarg;
      init->rjrh = 1.0/(atof(value));        //Handled as rj/rh for ease.
      break;   
    case ('M'):                              //mass of central galaxy (M_sun)
      value = optarg;
      init->galaxy.M = atof(value);
      break;
    case ('d'):                              //distance from galaxy CoM (/1000)
      value = optarg;
      init->galaxy.R = atof(value)*1000;     //1000 - input in kpc
      break;  
    case ('v'):                              //Cluster orbital velocity (kms^-2)
      value = optarg;
      init->galaxy.v=atof(value);
      break;
    case ('l'):                              //core collapse mass-loss 
      value = optarg;                        //fraction retained.
      init->f_N = atof(value);
      init->s1 = false;
      break;  
    case ('x'):                              //core collapse expansion
      value = optarg;
      init->f_r = atof(value);
      init->s2 = false;
      break;  
    case ('e'):                              //core collapse change in mean mass
      value = optarg;
      init->f_mm = atof(value);
      init->s3 = false;
      break;  
    case ('z'):                              //adjustable zeta
      value = optarg;
      init->zeta = atof(value);
      break;   
     }
  set_cluster(init);
  set_galaxy(init);
  early_evolution(init);

  init->galaxy.M =  init->galaxy.M*init->mm0; //Converts to mean masses
}

void set_cluster(variables *init){
  /*Checks if user has specified initial N and mm. If these are unspecified, a 
    default N=100000, mm=0.5M_sun cluster is used. This  are declared to the
    cerr output.*/
  if ( init->N0 == 0 ) {
    cerr << "No N0 supplied. Using N=100000." << endl;
    init->N0 = 1e5;
  }
  if ( init->units == 0 ){
    init->mm0 = 1.0/init->N0;                 //Has to occur for N-body units
   }
  else if ( init->units == 1 && init->mm0 == 0 ){
      cerr << "No mean mass supplied. Using m=0.5M_sun." << endl;   
      init->mm0 = 0.5;                        //0.5M_sun if real units.
  }
}

void set_galaxy(variables *init){
  /*Checks the user specified tidal field properties. If specified, uses R_G 
    and V_G by preference. If no R_G and/or MG, checks for velocity, then rjrh
    ratio, using arbitatary properties if necessary.*/
                                        
  if (init->galaxy.v != 0) galaxy_vel(init); //If user specifies velocity
  if (init->r0 == 0 || init->units == 0) 
    galaxy_r(init);                          //If user doesn't specify rh
  if (init->rjrh != 0) galaxy_rhrj(init);    //If user specifies rhrj
  if (init->galaxy.M == 0 || init->galaxy.R == 0)
    galaxy_all(init); 

  if (r_j(init)/init->r0 < 0.769){            //r/rj = 0.77rh/rj so > 1/10
    cerr << endl;
    cerr << "Invalid tidal field strength (too strong for prediction): R = ";
    cerr <<setprecision(3) << init->r0/r_j(init) << endl;
    exit(1);
  }
}
  
  void galaxy_vel(variables *in){
  /* Using specified orbital velocity and any known galaxy properties, 
     gives remaining galaxy parameters required for the code. 
     USES CIRCULAR ORBITS ONLY!!!!*/

  double G = 1;
  if (in->units == 1) G = 4.30172e-3   ;         //(parsec*(km/s)^2)
  
  if ( in->galaxy.R == 0 && in->galaxy.M != 0)
    in->galaxy.R = (2.0*G*in->galaxy.M)/pow(in->galaxy.v,2);
  else if ( in->galaxy.M == 0 && in->galaxy.R != 0)
    in->galaxy.M = (pow(in->galaxy.v,2)*in->galaxy.R)/(2.0*G);
  else if ( in->galaxy.M == 0 && in->galaxy.R == 0) {
    cerr << "No galaxy mass or galactocentric radius supplied -";
    cerr << "Using MG=1e10M_sun" << endl;
    in->galaxy.M = 1e10;
    in->galaxy.R = (2.0*G*in->galaxy.M)/pow(in->galaxy.v,2);
  }
  else
    cerr << "Velocity, M_G and R_G specified. Ignoring velocity" << endl;
 }

void galaxy_r(variables *in){
  /* Sets initial cluster radius, based on input criteria (if radius is not
     user defined.*/
  if (in->units == 0 ) {
    if (in->r0 > 0)
      cerr << "Option r ignored (working in N-body units); for N-body units ";
    else cerr << "N-body units; setting " ;
    cerr << "rh = 0.77rv." << endl;
    in->r0=1.0;
  }
  else if ( in->rjrh != 0 && in->galaxy.M != 0 && in->galaxy.R != 0){
    in->r0 = r_j(in)/in->rjrh;
  }
  else{
    cerr <<"No initial radius supplied. Using r=3pc." << endl;
    in->r0 = 3.0;
  }      
}

void galaxy_rhrj(variables *in){
  /* Using the ratio of rhrj and any known galaxy properties, gives remaining 
     galaxy parameters required for the code. The rhrj is here defined for a 
     PLUMMER MODEL, ie rh/rj = 0.77*r/rj
     Therefore an rh/rj of 1/100 =/= R of 0.1*/
  if ( in->galaxy.M == 0 ){
    cerr << "Tidal field defined by r/rj." << endl;
    if ( in->galaxy.R == 0 ){
      cerr << "No galaxy mass or radius supplied. Using MG = 1e10" << endl;
      in->galaxy.M = 1e10;
      in->galaxy.R = pow((3.0*in->galaxy.M)/(in->N0),(1.0/3.0))*
	(in->rjrh*in->r0);
    }
    else{
      in->galaxy.M = ((in->N0)/3.0)*
	pow(in->galaxy.R/(in->rjrh*in->r0),3);
    }
  }
  else if ( in->galaxy.R == 0 ){
    cerr << "Tidal field defined by r/rj." << endl;;
    in->galaxy.R = pow((3.0*in->galaxy.M)/(in->N0),(1.0/3.0))*
      (in->rjrh*in->r0);
  }
  else if (in->galaxy.v != 0)
    cerr << "Galaxy properties defined by velocity. Ignoring rj/rh" << endl;
  else 
    cerr << "MG and RG already specified. Ignoring rj/rh." << endl;
}

void galaxy_all(variables *in){
  /* Final check - if only a mass or radius defined, uses default values. If one
     is specified, uses that, and uses an arbitrary value for the other.*/

  if (in->galaxy.M != 0){
     cerr << "No galaxy radius supplied. Setting RG = 8.5kpc" << endl; 
     in->galaxy.R = 8500;
  }
  else if (in->galaxy.R != 0){
    cerr << "No galaxy mass supplied. Setting MG = 1e10M_sun" << endl;  
    in->galaxy.M = 1e10;
  }
  else{
    cerr << "No galaxy data supplied.";
    cerr << "Setting RG = 8.5kpc, v=220 kms^-2" << endl;
    in->galaxy.M = 4.8e10;
    in->galaxy.R = 8500;
  }
  in->rjrh = in->r0/r_j(in);
}

double r_j(variables *in){
  /*Jacobi radius function (ie, find r_j*/
  return pow((in->N0)/(3.0*in->galaxy.M),(1.0/3.0))*in->galaxy.R;
}

void early_evolution(variables *init){
  /* Basic prescriptions for pre-core collapse behaviour. This is based on 
     figure 7 in AG2012, and is roughly accurate for 20 < R < infinity. This
     can be extended to R = 10, although accuracy suffers at this point.*/
  double R = init->r0/(r_j(init)/0.77);   //0.77 on account of half-mass radius
  if (init->s1) init->f_N = 1.95-exp(80.0*pow(R,2.4));
  if (init->s2) init->f_r = 1.75*exp(-2.2*pow(R,0.55)); 
  /*Additional factor of 0.77 present on account of conversion from half mass 
    to virial radius during the pre-collapse expansion. Thus, an initially 
    quoted half mass radius is converted to virial radius by this factor.*/
  if (init->s3) init->f_mm = 1.0;
}

void zero(variables *in){
  /*Initialises a cluster to unphysical values. These are default flags used
    by the code to detect which properties are user defined.*/
  cerr << endl; 
  cerr << "Evolve Me A Cluster of Stars v1.03" << endl;
  cerr << "-----------------------------------------"<<endl; 

  in->N0 = 0;
  in->r0 = 0;              
  in->mm0 = 0;              
  in->galaxy.M = 0;         
  in->galaxy.R = 0;         
  in->tcc = 20.0*sqrt(pow(0.77,3)); //20*trh, where trh0 is defined by rh=0.77rv
  in->zeta = 0.111;                 //Default (for equal mass cluster)
  in->s1 = in->s2 = in->s3 = true;  //Sets quantities to be defined
  in->units = 1;                    //Physical units
  in->rjrh = 0;
  in->galaxy.v = 0;

  /*Note - the core collapse time, 20*trh0, uses the exact value of rh = 0.77rv
    (Plummer model) for the scaling of relaxation time (trh0). This occurs the 
    core collapse time for a Plummer model is given as 17.6 trh0 (Takahashi 
    1995), where rh is defined as 0.77rv. Hence, we use 20*(0.77)^(3/2) to 
    account for this scaling.*/
}

void help(){
  cout << "\t Evolve Me A Cluster of StarS (EMACSS) - version 1.03";
  cout << '\n';
  cout << "\t By: Poul Alexander (University of Cambridge)\n ";
  cout << "\t     Mark Gieles (University of Cambridge)\n";
  cout << '\n';
  cout << "\t EMACSS is a numerical integrator that solves for the\n";
  cout << "\t differential equations:\n";
  cout << '\n';
  cout << "\t  dN/dt = f(N,r,p) = - xi*N/t_rh\n";
  cout << '\n';
  cout << "\t  dr/dt = g(N,r,p) = mu*r/t_rh\n";
  cout << '\n';
  cout << "\t where \n";
  cout << "\t xi = f(N,r,rt)\n";
  cout << "\t mu = zeta - 2*xi\n";
  cout << "\t p is an array of variables defining the conditions at any given time.\n";
  cout << "\t zeta defines the fractional flow of energy per half-mass\n";
  cout << "\t relaxation time.\n";
  cout << '\n';
  cout << "\t Full documentation is available in Alexander & Gieles 2012\n";
  cout << '\n';
  cout << "\t The options to specify a cluster (initial number of stars,\n";
  cout << "\t half mass radius) can be specified with\n";
  cout << "\t appropriate options (see below) \n";
  cout << "\t The output is printed directly onto the command line; we \n";
  cout << "\t therefore suggest piping this output into a file or plotting\n";
  cout << "\t package.\n";
  cout << '\n';
  cout << "\t emacss [options] > [output]\n";
  cout << '\n';
  cout << "\t Options:\n";
  cout << "\t -help: Displays this file.\n";
  cout << "\t -version: Displays the current version number and version log.\n";
  cout << "\t -u: (Input) N-body or Real units. \n";
  cout << "\t\t 0=N-body (see Heggie & Mathieu 1986)\n";
  cout << "\t\t 1=Real [parsec, solar masses, Myr] \n";
  cout << "\t -N: initial number of stars in cluster.\n";
  cout << "\t -r: initial half-mass radius of cluster [pc if u = 1\n";
  cout << "\t     ignored if u = 0 ]\n";
  cout << "\t -m: mean mass of stars [solar masses].\n";
  cout << "\t -t: Start of balanced evolution (default is 20trh0).\n";
  cout << "\t -R: initial ratio of virial to Jacobi radius.\n";
  cout << "\t -M: mass of point mass galaxy around which the cluster is\n";
  cout << "\t     in orbit.\n";
  cout << "\t -d: distance to the point mass galaxy around which the cluster\n";
  cout << "\t     is in orbit [kpc if u = 1, krv (virial radius * 1000)\n";
  cout << "\t     otherwise].\n";
  cout << "\t -v: Orbital velocity of cluster around galaxy. [km/s]\n";
  cout << "\t -l: mass loss prior to core collapse.\n";
  cout << "\t -x: expansion prior to core collapse.\n";
  cout << "\t -e: change in mean mass prior to core collapse.\n";
  cout << "\t -z: value for zeta (otherwise default, 0.111 (equal mass \n";
  cout << "\t     clusters) is used).\n";
  cout << '\n';
  cout << "\t Tidal field is specified first by MG and RG, if these are\n";
  cout << "\t provided. If one or both are missing, the code will first \n";
  cout << "\t attempt to use orbital velocity, then filling factor rhrj.\n";
  cout << '\n';
  cout << "\t When options are not defined a default cluster is used with:\n";
  cout << "\t N=1e5\n";
  cout << "\t mm=0.5 M_sun\n";
  cout << "\t r=3pc\n";
  cout << "\t RG=8.5kpc\n";
  cout << "\t vG=220kms^-2\n";
  cout << '\n';
  cout << "\t If used to represent clusters with mass functions, we\n";
  cout << "\t recommend increasing zeta to ~0.2.\n";
  cout << '\n';
  cout << "\t Early evolution is not properly recovered by this \n";
  cout << "\t prescription- evolution prior to t_cc is found via exponential\n";
  cout << "\t extrapolation, and is marked at output by a 1 in the final \n";
  cout << "\t column.\n";
  cout << '\n';

  exit(1);
}

void version(){
  cout << "\t Evolve Me A Cluster of StarS (EMACSS) - version 1.03\n";
  cout << '\n';
  cout << "\t Version log:\n";
  cout << "\t 0.1 - Basic Outline\n";
  cout << "\t 0.11 - Functions added to define required quantities.\n";
  cout << "\t 0.12 - Rearranged evolve routine using pointers.\n";
  cout << "\t 0.13 - Input data added (hard wired), output to file.\n";
  cout << "\t 0.14 - Structures défient types.\n";
  cout << "\t 0.15 - Adjustment to use rk4 routine. Correction to t_rh.\n";
  cout << "\t 1.00 - Completed initial version (01/11)\n";
  cout << "\t 1.01 - Parameter read migrated to external file, option for\n";
  cout << "\t        input from command line supported.\n";
  cout << "\t 1.02 - Inclusion mass loss in isolated regime (this is xi0.\n";
  cout << "\t        equation (3), Baumguardt, Hut, Heggie).\n ";
  cout << "\t 1.03 - Version published with AG2012 on ArXiv (20/03/2012)\n";
  cout << '\n';

  exit(1);
}
