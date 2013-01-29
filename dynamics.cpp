#include "../emacss.h"

void dynamics::set_tcc(){
  T_DYN = T_DYN*mynode->t_relax.nbody;
}

/*---------------------------------------------------------------------------*/
// Purely dynamical differential equations (see AG2012, s2 & AGLB2013, s?)
double dynamics::dNdt(){
  double dNdt = 0;
  if (mynode->time.nbody > T_DYN) 
      dNdt -= (xi()*mynode->N)/mynode->t_relax.nbody;  //Equation (6) AG2012
  return dNdt;
}

double dynamics::dmmdt(){
  double dMdt = 0;
  dMdt -= ((gamma_dyn()+myse->gamma_se())*mynode->mm.nbody)/
	mynode->t_relax.nbody;
  return dMdt;
}

double dynamics::drdt(){
  double drdt = 0;
  drdt += (mu()*mynode->r.nbody)/mynode->t_relax.nbody;//Equation (7) AG2012
  return drdt;
}

double dynamics::dkdt(){                               //This is not used
  double dkdt = 0;
  dkdt = lambda()*mynode->kappa/mynode->t_relax.nbody;
  return dkdt;
}

double dynamics::dtrhdt(){
  return 1.0/mynode->t_relax.nbody;                          //t_rh per unit time
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
//Dynamical Dimensionless Parameter Equations (see AG2012, s2.2 & AGLB2013, s?)

double dynamics::xi(){                           //Equation (26) AG2012	
  double xi = 0, tidal = 0, iso = 0;
  if (mynode->time.nbody > T_DYN) {
      tidal = (3.0/5.0)*mynode->E.zeta;
      iso = xi0;
  }
  return tidal*P()+iso*(1.0-P());
}

double dynamics::gamma_dyn(){                   //Equation (??) AGL2013
  double gamma = 0;
  if (mynode->time.nbody > T_DYN) 
      gamma -= a*xi()/(1.0+((1.0-mynode->DM_SE.nbody)/
		(F*mynode->N*mynode->mm.nbody)));
  return gamma;
}

double dynamics::lambda(){                      //Not yet in use
  double lambda =0;
  return lambda;
}

double dynamics::mu(){                //Equation (8) AG2012
  double mu;
  mu = myse->epsilon() - 2.0*xi() - 2.0*gamma_dyn() - 2.0*myse->gamma_se() - lambda();
  return mu;
}

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
//Extra required parameter equations (see AG2012, s2 & AGLB2013, s?)

double dynamics::P(){                                     //Equation (25) AG2012
  double f = pow((mynode->N*log(mynode->gamma*N1))/\
	(N1*log(mynode->gamma*mynode->N)),(1.0-x));
  double g = pow((mynode->r.nbody/mynode->rj.nbody)/R1,z);
  return f*g;
}
