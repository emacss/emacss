#include "../emacss.h"

void stellar_evo::tse(double M_star, double T_star){
 if (mynode->s == 0){
     m_ns = mynode->mm; m_ref = mynode->mm; t_ref = mynode->mm;
     t_se = numeric_limits<double>::max();
  }
 else{
   m_ns = m_ns/M_star; m_ref = m_ref/M_star; t_ref = t_ref/T_star;
   t_se = t_ref*pow(mynode->m_max/m_ref,\
		-1.7+0.82*log10(mynode->m_max*M_star));
   psi1 = psi1*pow(t_se/t_ref,y); 
 }
}


/*---------------------------------------------------------------------------*/
// Purely SE differential equations (AGLB2013, s?)

double stellar_evo::dmsedt(){                            //Equation (?) AGLB2013
  double dmsedt = 0;
  dmsedt -= gamma_se()*mynode->mm/mynode->t_rhp;            
  return dmsedt;
}

double stellar_evo::dmsegdt(){                            //Equation (?) AGLB2013
  double dmsegdt = 0;
  dmsegdt += chi()*mynode->MS/mynode->t_rhp;            
  return dmsegdt;
}

double stellar_evo::dEdt(){                       
    double dEdt = 0;
    dEdt = epsilon()*fabs(mynode->E.value)/mynode->t_rhp;
  return dEdt;
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
//Dynamical SE Parameter Equations (AGLB2013, s?)

double stellar_evo::epsilon(){                           //Equation (?) AGLB2013
  double epsilon = 0;
   
  if (mynode->s == 0) epsilon = -mydyn->lambda() + mydyn->mu() + \
	2.0*mydyn->xi();
  
  else{
    mynode->E.source = 0;                                //Needed due to adaptive timestep
    if (mynode->trhpelapsed > mydyn->nc){
      mynode->E.source = 2; 
      epsilon = mynode->E.zeta;
    }
    else if (mynode->time > t_se){
      mynode->E.source = 1;
      epsilon = mynode->MS*gamma_se() + tidal_escape();	 
    }
    else epsilon = tidal_escape();
  }
 return epsilon;
}

double stellar_evo::gamma_se(){                          //Equation (?) AGLB2013
  double gamma = 0;
  if (mynode->s == 1 && mynode->E.source > 0)
      
    gamma += (nu*mynode->t_rhp)/mynode->time*(mynode->mm_se/mynode->mm);

  return gamma;
}

double stellar_evo::chi(){
    double chi = 0;
    
    if (mynode->s == 1) chi += (MS_1-mynode->MS);

    return chi;
}

/*---------------------------------------------------------------------------*/
  
  double stellar_evo::psi(){
  double psi = 1;

  if (mynode->s == 1){
    psi = psi1;
    if (mynode->E.source > 0)
      psi =(psi1-psi0)*pow(mynode->time/t_se,y)+psi0;
  }

  return psi;
}

double stellar_evo::tidal_escape(){
    double epsilon = 0, A = 0;
    
    if (mynode->galaxy.type > 0){
       A = (1.0/mynode->kappa)*(mynode->Rhj);
    }
    epsilon = A*mydyn->xi()*(mydyn->m_ej()/mynode->mm);
    
    return epsilon;
}

double stellar_evo::m_up(){
    
  double m_up = mynode->m_max;
  
  if (mynode->E.source > 0)
    m_up = mynode->m_max0*pow(mynode->m_max0/m_ns,\
	    pow(t_se/mynode->time,0.37) - 1.0); 
  
  m_up = sqrt(pow(m_up,2)+pow(m_ns,2));
  
  return m_up;
}
