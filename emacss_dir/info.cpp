#include "../emacss.h"

/*Readme and help files these are not essential for program operation, but
 contain a guide for use, version information, and sufficient documentation
 for basic usage. For full instruction, see Alexander and Gieles 2012 and
 Alexander, Gieles and Lamers 2013.*/

void help(){
  cout << "\t Evolve Me A Cluster of StarS (EMACSS) - version 3.10";
  cout << '\n';
  cout << "\t By: Poul Alexander (University of Cambridge)\n ";
  cout << "\t     Mark Gieles (University of Surrey)\n";
  cout << "\t     Henny Lamers (University of Amsterdam)\n";
  cout << "\t     Holger Baumgardt (University of Brisbane)\n";
  cout << '\n';
  cout << "\t EMACSS is a numerical integrator that solves for the\n";
  cout << "\t differential equations:\n";
  cout << '\n';
  cout << "\t  dE/dt = f(N,r,p) = - epsilon*E/t_rh'\n";
  cout << '\n';
  cout << "\t  dN/dt = g(N,r,p) = - xi*N/t_rh'\n";
  cout << '\n';
  cout << "\t  drh/dt = h(N,r,p) = mu*r/t_rh'\n";
  cout << '\n';
  cout << "\t  drc/dt = h(N,r,p) = delta*r/t_rh'\n";
  cout << '\n';
  cout << "\t  dmm/dt = i(N,r,p) = gamma*mm/t_rh'\n";
  cout << '\n';
  cout << "\t  dk/dt = h(N,r,p) = lambda*r/t_rh'\n";
  cout << '\n';
  cout << "\t where \n";
  cout << "\t epsilon = -lambda+mu+2xi \n";
  cout << "\t xi = f(N,r,rt)\n";
  cout << "\t mu = f(rh,rt,rc,N,kappa)\n";
  cout << "\t delta = f(t_rh,t_rc,N,xi)\n";
  cout << "\t lambda = f(rc,rh,kappa)\n";
  cout << "\t gamma = f(N,r,rt,t)\n";
  cout << "\t t_rh' = t_rh/psi(t)\n"; 
  cout << '\n';	
  cout << "\t p is an array of variables defining the conditions at any given time.\n";
  cout << "\t psi is a decreasing function of time.\n";
  cout << "\t zeta defines the fractional flow of energy per half-mass\n";
  cout << "\t relaxation time and stellar evolution.\n";
  cout << '\n';
  cout << "\t Full documentation is available in Alexander & Gieles 2012, \n";
  cout << "\t Gieles, Alexander, Lamers and Baumgardt 2013, and \n";
  cout << "\t Alexander, Gieles, Lamers and Baumgardt 2014.\n";
  cout << '\n';
  cout << "\t The options to specify a cluster (initial number of stars,\n";
  cout << "\t half mass radius) can be specified with appropriate options \n";
  cout << "\t (see below) \n";
  cout << "\t The output is printed directly onto the command line; we \n";
  cout << "\t therefore suggest piping this output into a file or plotting\n";
  cout << "\t package.\n";
  cout << '\n';
  cout << "\t emacss [options] > [output]\n";
  cout << '\n';
  cout << "\t Although note that the output is formatted (see 'head' for the \n";
  cout << "\t format of output files. \n"; 
  cout << '\n';
  cout << "\t Options:\n";
  cout << "\t -help: Displays this file.\n";
  cout << "\t -version: Displays the current version number and version log.\n";
  cout << "\t -o: (Output) N-body or Real units. \n";
  cout << "\t\t 0=N-body (see Heggie & Mathieu 1986)\n";
  cout << "\t\t 1=Real [parsec, solar masses, Myr] \n";
  cout << "\t\t 2=Real and N-body alternating \n";
  cout << "\t -N: initial number of stars in cluster.\n";
  cout << "\t -r: initial half-mass radius of cluster. \n";
  cout << "\t -m: mean mass of stars [solar masses].\n";
  cout << "\t -R: initial ratio of half-mass to Jacobi radius.\n";
  cout << "\t -g: Galaxy type: 0 - isolated, 1 = in tidal field. Tidal field\n";
  cout << "\t	  assumed to be point mass if s = 0, isothermal otherwise.\n"; 
  cout << "\t -M: mass of point mass galaxy around which the cluster is\n";
  cout << "\t     in orbit [units set by -o].\n";
  cout << "\t -d: distance to the point mass galaxy around which the cluster\n";
  cout << "\t     is in orbit [units set by -o].\n";
  cout << "\t -v: Orbital velocity of cluster around galaxy. [units set by \n";
  cout << "\t     -o]\n";
  cout << "\t -z: value for zeta. Defaults to 0.111 for equal mass clusters.\n";
  cout << "\t -s: Stellar Evolution flag (0 = off, equal mass, 1 = on, IMF).\n";
  cout << "\t -u: Upper limit of IMF (M_sun)\n";
  cout << "\t -l: Lower limit of IMF (M_sun)\n";
  cout << "\t -f: (Experimental) Dynamical Friction flag (0 = off, 1 = on).\n";
  cout << '\n';
  cout << "\t Tidal field is specified first by MG and RG, if these are\n";
  cout << "\t provided. If one or both are missing, the code will first \n";
  cout << "\t attempt to use orbital velocity, then filling factor rhrj.\n";
  cout << '\n';
  cout << "\t When options are not defined a 64k isolated cluster is used.\n";
  cout << '\n';
  cout << "\t If used to represent clusters of equal mass stars, the code\n";
  cout << "\t is only accurate from the start (i.e., throughout core \n";
  cout << "\t collapse.\n";
  cout << '\n';
  cout << "\t The default output is N-body units. If option 2 is used\n";
  cout << "\t N-body and real units alternate. \n";
  cout << '\n';
  cout << "\t Happy evolving!\n";
  cout << '\n';

  exit(1);
}

void version(){
  cout << "\t Evolve Me A Cluster of StarS (EMACSS) - version 2.02\n";
  cout << '\n';
  cout << "\t Version log:\n";
  cout << "\t 0.10 - Basic Outline\n";
  cout << "\t 0.11 - Functions added to define required quantities.\n";
  cout << "\t 0.12 - Rearranged evolve routine using pointers.\n";
  cout << "\t 0.13 - Input data added (hard wired), output to file.\n";
  cout << "\t 0.14 - Structures define types.\n";
  cout << "\t 0.15 - Adjustment to use rk4 routine. Correction to t_rh.\n";
  cout << "\t 1.00 - Completed initial version (01/11)\n";
  cout << "\t 1.01 - Parameter read migrated to external file, option for\n";
  cout << "\t        input from command line supported.\n";
  cout << "\t 1.02 - Inclusion mass loss in isolated regime (this is xi0.\n";
  cout << "\t        equation (3), Baumguardt, Hut, Heggie).\n ";
  cout << "\t 1.03 - Version published with AG2012 on ArXiv (20/03/2012)\n";
  cout << "\t 1.04 - Minor update, with addition of specific time output\n";
  cout << "\t        for finding cluster properties after a specific time \n";
  cout << "\t        (20/03/2012)\n";
  cout << "\t 2.00 - Initial introduction of unbalanced (pre-core collapse)\n";
  cout << "\t        evolution (change in core radius). \n";
  cout << "\t 2.01 - Inclusion of varying kappa, and energetic effects of  \n";
  cout << "\t        unbalanced evolution.\n";
  cout << "\t 2.02 - Completed second version (07/13). \n";
  cout << '\n';

  exit(1);
}
