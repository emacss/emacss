#include "../emacss.h"

void stellar_evo::set_X(){
  X = 1.0 + 0.72*pow(log10(mynode->t_relax.Myr)-3.5,2);
}

/*---------------------------------------------------------------------------*/
// Purely SE differential equations (AGLB2013, s?)

double stellar_evo::dmsedt(){                            //Equation (?) AGLB2013
  double dmsedt = 0;
  dmsedt = -(gamma_se()*mynode->DM_SE.nbody)/mynode->t_relax.nbody;            
  return dmsedt;
}

double stellar_evo::dEdt(){                       
    double dEdt = 0;
    dEdt = epsilon()*fabs(mynode->E.nbody)/mynode->t_relax.nbody;
  return dEdt;
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
//Dynamical SE Parameter Equations (AGLB2013, s?)

double stellar_evo::epsilon(){                           //Equation (?) AGLB2013
  double epsilon = 0;
  mynode->E.source = 0;        //No energy generation
  epsilon = (2.0+X)*gamma_se();
  if (epsilon < mynode->E.zeta && mynode->time.Myr > T_SE){
     mynode->E.source = 2;
     epsilon = mynode->E.zeta;
  }
  return epsilon;
}

double stellar_evo::gamma_se(){                          //Equation (?) AGLB2013
  double gamma = 0;
  if (mynode->s != 0 && mynode->time.Myr > T_SE){
      gamma += (nu*mynode->t_relax.nbody)/mynode->time.nbody;
      mynode->E.source=1;
  }
  return gamma;
}

/*---------------------------------------------------------------------------*/
