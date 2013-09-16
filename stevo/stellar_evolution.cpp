#include "../emacss.h"

/*---------------------------------------------------------------------------*/
// Purely SE differential equations (AGLB2013, s?)

double stellar_evo::dmsedt(){                            //Equation (?) AGLB2013
  double dmsedt = 0;
  dmsedt -= (gamma_se()*mynode->mm.nbody)/mynode->t_relax.nbody;            
  return dmsedt;
}

double stellar_evo::dmsegdt(){                            //Equation (?) AGLB2013
  double dmsegdt = 0;
  dmsegdt += (chi()*mynode->mass_seg*zeta())/mynode->t_relax.nbody;            
  return dmsegdt;
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
  double epsilon, adi_factor = 0;
  
  mynode->E.source = 0;					  //No energy generation
  adi_factor = (2.0+mynode->mass_seg);
  
  if (mynode->time.Myr > mynode->T_SE){
      mynode->E.source = 1;	
      epsilon = adi_factor*gamma_se() + tidal_escape();
  }

  if (mynode->trhp > mynode->T_DYN){
      mynode->E.source = 2; 
      epsilon = zeta();
  }
  
  return epsilon;
}

double stellar_evo::gamma_se(){                          //Equation (?) AGLB2013
  double gamma = 0;
  if (mynode->s != 0 && mynode->time.Myr > mynode->T_SE)
      gamma += (nu*mynode->t_relax.nbody)/mynode->time.nbody\
	      *(mynode->mm_se.nbody/mynode->mm.nbody);
  return gamma;
}

double stellar_evo::chi(){
    double chi = 0;
    chi += (MS_max()-mynode->mass_seg);
    return chi;

}

double stellar_evo::MS_max(){
    double MS = MS_1;
    if (mynode->time.Myr > 0)
    	MS = (MS_1);//*pow(mynode->time.Myr/mynode->T_SE,y)+1.0);
    return MS;
}


/*---------------------------------------------------------------------------*/

double stellar_evo::zeta(){
    double zeta = z0;
    
    if (mynode->time.Myr > mynode->T_SE){
      zeta =(z0-mynode->E.zeta)*pow(mynode->time.Myr/mynode->T_SE,y)\
		+mynode->E.zeta;
    }
    return zeta;
}

double stellar_evo::tidal_escape(){
    double epsilon = 0, A = 0;
    
    if (mynode->galaxy.type > 0){
       A = (1.0/mynode->kappa)*(mynode->rhrj);
    }
    epsilon = A*mydyn->xi();  
    return epsilon;
}

