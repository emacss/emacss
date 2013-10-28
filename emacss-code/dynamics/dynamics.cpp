#include "../emacss.h"


/*---------------------------------------------------------------------------*/
// Purely dynamical differential equations (see AG2012, s2 & AGLB2013, s?)
double dynamics::dNdt(){
  double dNdt = 0;
  dNdt = -(xi()*mynode->N)/mynode->t_rhp.nbody;  //Equation (6) AG2012
  return dNdt;
}

double dynamics::dmmdt(){
  double dmdt = 0;
  dmdt = -gamma()*mynode->mm.nbody/mynode->t_rhp.nbody;  
  return dmdt;
}

double dynamics::drdt(){
  double drdt = 0;
  drdt = (mu()*mynode->r.nbody)/mynode->t_rhp.nbody;//Equation (7) AG2012
  return drdt;
}

double dynamics::dkdt(){                               //This is not used
  double dkdt = 0;
  dkdt = -lambda()*mynode->kappa/mynode->t_rhp.nbody;
  return dkdt;
}

double dynamics::dtrhdt(){                      //t_rh per unit time
  return 1.0/mynode->t_relax.nbody;       
}

double dynamics::dtrhpdt(){                            //t_rh' per unit time
  return 1.0/(mynode->t_rhp.nbody);   
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
//Dynamical Dimensionless Parameter Equations (see AG2012, s2.2 & AGLB2013, s?)


double dynamics::xi(){                           //Equation (26) AG2012	
  double xi = 0;
    xi = xi_i()+xi_e();
  return xi;
}

double dynamics::xi_e(){                           //Equation (26) AG2012	
 
  double xi = 0, F = 1, f_ind = 0;
  if (mynode->trhp > mynode->T_DYN || mynode->time.nbody < mynode->tcc.nbody ){
      F = pow(mynode->trhp/mynode->T_DYN,3);
  }
  xi += F*xi0*(1.0-P())+(f+(1-f)*F)*(3.0/5.0)*myse->zeta()*P();
  return xi/myse->zeta();
}

double dynamics::xi_i(){
  double xi = 0;
    if (mynode->galaxy.type > 0 && mynode->trhp < mynode->T_DYN) 
	xi += f_ind()*myse->gamma_se();
  return xi;
}


double dynamics::gamma(){
  double gamma = 0;
  gamma = myse->gamma_se()-gamma_dyn();
  return gamma;  
}

double dynamics::gamma_dyn(){                   //Equation (??) AGL2013
  double gamma = 0;
  if (mynode->galaxy.type > 0)
    gamma += F*xi();
  return gamma;
}

double dynamics::lambda(){                      //Not yet in use
  double lambda, F1;
  F1 = (k_1-mynode->k_0)/mynode->T_DYN;
  lambda = F1*pow(mynode->trhp/mynode->T_DYN,3);
//  lambda = (mynode->k_0-k_1)/mynode->T_DYN;
  if (fabs(mynode->kappa-k_1) < 1e-2) lambda = 0;
  
  return -lambda/myse->zeta();
}

double dynamics::mu(){                          //Equation (8) AG2012
  double mu;
  mu = myse->epsilon() - 2.0*xi() - lambda()- 2.0*gamma();
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


double dynamics::f_ind(){
  double f= 0;
  
 // if (mynode->rhrj > 0.8*R1) 
//      f = 2.5*(mynode->rhrj);
  if (mynode->rhrj > R1) 
      f = (exp(pow(mynode->rhrj/(1.4*R1),5))-1.0);
  
  return f;
}
