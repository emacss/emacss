#include "../emacss.h"

/*Readme and help files these are not essential for program operation, but
 contain a guide for use, version information, and sufficient documentation
 for basic usage. For full instruction, see Alexander and Gieles 2012 and
 Alexander, Gieles and Lamers 2013.*/

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
