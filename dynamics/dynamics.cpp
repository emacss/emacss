#include "../emacss.h"

void dynamics::set_k0(){ 
  // Called at start to set the constant k0 needed in the function K()
  mynode->k0 = mynode->rh/4.0;
}

double dynamics::Rchmin(){          
  // Definition of the minimum ratio r_c/r_h that defines core collapse
  double Rchmin = 0;
  Rchmin += pow(N2/mynode->N + N2/N3,2.0/3.0);
  return Rchmin;
}

void dynamics::reset_K_constants(){          
  //Changes parameters to post-collapse values.
  mynode->Rch0 = 0.2;
  mynode->k0 = 0.19;
  mynode->k1 = 0.27;
}

/*---------------------------------------------------------------------------*/
// Purely dynamical differential equations (see AG2012, s2 & GALB2013, s2)

double dynamics::dEdt(){			//Equation (4) GALB2013 
  double dEdt = 0;
  dEdt -= epsilon()*mynode->E.value/mynode->t_rh;
  return dEdt;
}

double dynamics::dNdt(){
  double dNdt = 0;
  dNdt -= xi()*mynode->N/mynode->t_rh;		//Equation (6) AG2012
  return dNdt;
}

double dynamics::dmmdt(){			//Not used
  double dMdt = 0;
  return dMdt;
}

double dynamics::drdt(){			//Equation (7) AG2012
  double drdt = 0;
  drdt += mu()*mynode->rh/mynode->t_rh;
  return drdt;
}

double dynamics::dkdt(){			//Equation (5) GALB2013
  double dkdt = 0;
  dkdt += lambda()*mynode->kappa/mynode->t_rh;
  return dkdt;
}

double dynamics::dtrhdt(){			//Monitors t_rh per unit time
  double dtrhdt = 0;
  dtrhdt += 1.0/mynode->t_rh;
  return dtrhdt;			
}

double dynamics::drcdt(){			//Equation (10) GALB2013
  double drcdt = 0;
  drcdt += delta()*mynode->rc/mynode->t_rh;
  return drcdt;                          
}

/*---------------------------------------------------------------------------*/
//Dynamical Dimensionless Parameter Equations (see AG2012, s2 & GALB2013, s2,3)

double dynamics::xi(){                          //Equation (26) AG2012
  double xi = 0, tidal = 0, iso = 0;	  	//Equation (22) GALB2013

  tidal = (3.0/5.0)*mynode->E.zeta;
  iso = xi1*mynode->E.zeta/0.1; 

  if (mynode->E.source == 0){
    double F = 0.0;   

    F = Rchmin()/mynode->Rch;

    xi = (f+(1-f)*F)*tidal*P();
    xi += F*iso*(1-P());
  }
  
  if (mynode->E.source > 0){ 
    xi = tidal*P() + iso*(1.0 - P());
  }
  
  return xi;
}

double dynamics::lambda(){			//Equation (12) GALB2013
  double lambda = 0;
  
  lambda += K()*(delta() - mu()); 
  if (mynode->E.source==1) lambda += (k() - mynode->kappa)/k();
 
   return lambda;
}

double dynamics::epsilon(){			//Equation (2) GALB2013
  double epsilon = 0.0;
  
  epsilon = -lambda() + mu() + 2.0*xi();
  
  return epsilon;
}

double dynamics::mu(){				//Equation (19) GALB2013
  double mu = 0.0, rhrj = 0;
  
  if (mynode->galaxy.type > 0) rhrj = mynode->Rhj;
  
  if (mynode->E.source == 0){
      mu = K()*delta();                       //Increase due to shrinking core
      mu += (rhrj/mynode->kappa - 2.0)*xi();//Shrink due to removal of stars
      mu /= (1.0 + K());// Correction because rh appears in expression for kappa
     }

   if (mynode->E.source == 1){
       mu = mynode->E.zeta - 2.0*xi();         // Make the balance (AG2012)
       mu += 2.0/3.0*K()*xi()/(1.0+mynode->N/N3);// Extra terms
     }
  return mu;
 }

 double dynamics::delta(){			//Equation (17) GALB2013
   double delta = 0;
   
   if (mynode->E.source == 0)
     delta += delta_1 + delta_2*mynode->t_rh/mynode->t_rc;
   if (mynode->E.source == 1)
     delta += mu() + 2.0/3.0*xi()/(1.0 + mynode->N/N3);

   return delta;
 }

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
//Extra required parameter equations (see AG2012, s2 & AGLB2013, s2,3)

double dynamics::P(){				 //Equation (25) AG2012
  double f = 0, g = 0;
  
  f += pow((mynode->N*log(mynode->gamma*N1))/  \
		 (N1*log(mynode->gamma*mynode->N)),(1.0-x));
  if (mynode->galaxy.type > 0)
    g += pow((mynode->rv/mynode->rj)/R1,z);
  
  return f*g;
}

double dynamics::K(){				 //Equation (29) GALB2013          
  double K = 0;
  
  K = 2.0*mynode->Rch*(mynode->k0-mynode->k1)*\
	exp(-pow(mynode->Rch/mynode->Rch0,2.0))/(sqrt(M_PI)*mynode->Rch0*k());
  return K;
}

double dynamics::k(){				 //Breakdown, GALB2013           
  double k = 0;
  
  k = (mynode->k1+(mynode->k0-mynode->k1)*erf(mynode->Rch/mynode->Rch0));
  return k;
}
