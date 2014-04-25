#include "../emacss.h"

/*Readme and help files these are not essential for program operation, but
 contain a guide for use, version information, and sufficient documentation
 for basic usage. For full instruction, see Alexander and Gieles 2012 and
 Alexander, Gieles and Lamers 2013.*/

void help(){
  cout << "\t Evolve Me A Cluster of StarS (EMACSS) - version 2.10";
  cout << '\n';
  cout << "\t By: Poul Alexander (University of Cambridge)\n ";
  cout << "\t     Mark Gieles (University of Surrey)\n";
  cout << "\t     Henny Lamers (University of Amsterdam)\n";
  cout << "\t     Holger Baumgardt (University of Brisbane)\n";
  cout << '\n';
  cout << "\t EMACSS is a numerical integrator that solves for the\n";
  cout << "\t differential equations:\n";
  cout << '\n';
  cout << "\t  dE/dt = f(N,r,p) = - epsilon*E/t_rh\n";
  cout << '\n';
  cout << "\t  dN/dt = g(N,r,p) = - xi*N/t_rh\n";
  cout << '\n';
  cout << "\t  dr/dt = h(N,r,p) = mu*r/t_rh\n";
  cout << '\n';
  cout << "\t  dmm/dt = i(N,r,p) = (gamma_se+gamma_dyn)*mm/t_rh\n";
  cout << '\n';
  cout << "\t where \n";
  cout << "\t epsilon = max((2 + X)*gamma_se,zeta) \n";
  cout << "\t xi = f(N,r,rt)\n";
  cout << "\t mu = epsilon - 2*xi - 2*gamma_se - 2*gamma_dyn\n";
  cout << "\t gamma_se = nu*trh/t";
  cout << "\t gamma_dyn = A*xi/(1+(1-DM_SE)/(F*N*mm))";	  
  cout << '\n';	
  cout << "\t p is an array of variables defining the conditions at any given time.\n";
  cout << "\t zeta defines the fractional flow of energy per half-mass\n";
  cout << "\t relaxation time and stellar evolution.\n";
  cout << '\n';
  cout << "\t DM_SE is the mass lost due to stellar evoltuion only\n";
  cout << '\n';
  cout << "\t X, A and F are constants calibrated against N-body runs\n"; 
  cout << '\n';
  cout << "\t Full documentation is available in Alexander & Gieles 2012 and \n";
  cout << "\t in Alexander, Gieles, Lamers and Baumgardt 2013 (in prep) and \n";
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
  cout << "\t -r: initial half-mass radius of cluster [ignored if u = 0 ]\n";
  cout << "\t -m: mean mass of stars [solar masses].\n";
  cout << "\t -t: Output time (i.e., for single output after set time).\n";
  cout << "\t -R: initial ratio of virial to Jacobi radius.\n";
  cout << "\t -M: mass of point mass galaxy around which the cluster is\n";
  cout << "\t     in orbit [M_sun].\n";
  cout << "\t -d: distance to the point mass galaxy around which the cluster\n";
  cout << "\t     is in orbit [kpc].\n";
  cout << "\t -v: Orbital velocity of cluster around galaxy. [km/s]\n";
  cout << "\t -z: value for zeta. Defaults:\n";
  cout << "\t\t 0.111-equal mass clusters, no stellar evolution.\n";
  cout << "\t\t 0.15-stellar evolution on, clusters assumed with Krouper \n";
  cout << "\t\t stellar IMF.\n";
  cout << "\t -s: Stellar Evolution flag (0 = off, equal mass, 1 = on, IMF.\n";
  cout << '\n';
  cout << "\t Tidal field is specified first by MG and RG, if these are\n";
  cout << "\t provided. If one or both are missing, the code will first \n";
  cout << "\t attempt to use orbital velocity, then filling factor rhrj.\n";
  cout << '\n';
  cout << "\t When options are not defined a default cluster is used with:\n";
  cout << "\t N=1e5\n";
  cout << "\t mm=0.547 M_sun\n";
  cout << "\t r=1pc\n";
  cout << "\t RG=8.5kpc\n";
  cout << "\t vG=220kms^-2\n";
  cout << '\n';
  cout << "\t If used to represent clusters of equal mass stars, the code\n";
  cout << "\t is only accurate once balanced evolution has begun (no \n";
  cout << "\t core-collapse).\n";
  cout << '\n';
  cout << "\t Happy evolving!\n";
  cout << '\n';

  exit(1);
}

void version(){
  cout << "\t Evolve Me A Cluster of StarS (EMACSS) - version 2.10\n";
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
  cout << "\t 2.00 - Initial introduction of mass loss through stellar \n";
  cout << "\t        evolution (change in mean mass) \n";
  cout << "\t 2.01 - Additional change in mean mass of stars through \n";
  cout << "\t        the preferential loss of low mass stars\n";
  cout << "\t 2.02 - Energy production (zeta) allowed to vary prior to \n";
  cout << "\t        balanced evolution.\n";
  cout << "\t 2.10 - Code restructured, test suite written. (28/01/2013)\n";
  cout << '\n';

  exit(1);
}
