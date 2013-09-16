#include "../emacss.h"


/*---------------------------------------------------------------------------*/
// Purely dynamical differential equations (see AG2012, s2 & AGLB2013, s?)
double dynamics::dNdt(){
  double dNdt = 0;
  dNdt = -((xi()+xi_ind())*mynode->N)/mynode->t_relax.nbody;  //Equation (6) AG2012
  return dNdt;
}

double dynamics::dmmdt(){
  double dmdt = 0;
  dmdt = -gamma()*mynode->mm.nbody/mynode->t_relax.nbody;  
  return dmdt;
}

double dynamics::drdt(){
  double drdt = 0;
  drdt = (mu()*mynode->r.nbody)/mynode->t_relax.nbody;//Equation (7) AG2012
  return drdt;
}

double dynamics::dkdt(){                               //This is not used
  double dkdt = 0;
  dkdt = -lambda()*mynode->kappa/mynode->t_relax.nbody;
  return dkdt;
}

double dynamics::dtrhdt(){                      //t_rh per unit time
  return 1.0/mynode->t_relax.nbody;       
}

double dynamics::dtrhpdt(){                            //t_rh' per unit time
  return myse->zeta()/(mynode->t_relax.nbody);   
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
//Dynamical Dimensionless Parameter Equations (see AG2012, s2.2 & AGLB2013, s?)

/*
double dynamics::xi(){                           //Equation (26) AG2012	
  double xi_nb = 0, xi_bal = 0, F = prec;
  xi_nb = (3.0/5.0)*myse->zeta()*P();
  xi_bal = xi_nb+xi0*(1.0-P());
  if (mynode->trhp > mynode->T_DYN) return xi_bal;
  else if (mynode->trhp > 0.5*mynode->T_DYN)
      F = (prec-1.0)*(1.0-mynode->trhp/mynode->T_DYN)+1;
  return xi_nb/F;
}
*/

double dynamics::xi(){                           //Equation (26) AG2012	
 
  double xi = 0, F = 1, f_ind = 0;
  if (mynode->trhp < mynode->T_DYN){
      F = pow(mynode->trhp/mynode->T_DYN,3);
  }
  xi += F*xi0*(1.0-P())+(f+(1-f)*F)*(3.0/5.0)*myse->zeta()*P();
  return xi;
}

double dynamics::xi_ind(){
  double ind = 0;
    if (mynode->galaxy.type > 0 && mynode->trhp < mynode->T_DYN) 
	ind = f_delay()*f_max();
  return ind;
}


double dynamics::gamma(){
  double gamma = 0;
  gamma = myse->gamma_se()-gamma_dyn();
  return gamma;  
}

double dynamics::gamma_dyn(){                   //Equation (??) AGL2013
  double gamma = 0;
  if (mynode->galaxy.type > 0) gamma += F*(xi()+xi_ind());
  return gamma;
}

double dynamics::lambda(){                      //Not yet in use
  double lambda = 0;
  lambda = (mynode->k_0-k_1)/mynode->T_DYN;
  if (fabs(mynode->kappa-k_1) < 0.01) lambda = 0;
  
  return lambda;
}

double dynamics::mu(){                          //Equation (8) AG2012
  double mu;
  mu = myse->epsilon() - 2.0*xi() - 2.0*xi_ind() - lambda()- 2.0*gamma();
//  cerr << mynode->time.Myr << ' ' << myse->epsilon() << ' ' << xi() << ' ' << gamma() << ' ' << mu << endl;
  return mu;
}

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
//Extra required parameter equations (see AG2012, s2 & AGLB2013, s?)

double dynamics::P(){                                     //Equation (25) AG2012
  double f = 0, g = 0;
  
  if (mynode->galaxy.type > 0){
      f += pow((mynode->N*log(mynode->gamma*N1))/\
	   (N1*log(mynode->gamma*mynode->N)),(1.0-x));  
      g += pow((mynode->r.nbody/mynode->rj.nbody)/R1,z);
  }
  return f*g;
}

double dynamics::f_delay(){
  double f_delay = 0, t_delay = 0;
  
  t_delay = n*mynode->tcrj.nbody;
  f_delay += 1.0-exp(-mynode->time.nbody/t_delay);
  
  return f_delay;
}

double dynamics::f_max(){
  double f_max = 0;
  
  if (mynode->rhrj > R1) 
      f_max = 2.5*(mynode->rhrj-R1);
  
  return f_max;
}
