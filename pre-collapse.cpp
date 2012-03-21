#include "emacss.h"

/*The pcc (pre-core collapse) module contains fitting functions used to interpolate evolution between the initial conditions of the cluster and the state of the cluster when balanced evolution is achieved. This functions in 2 ways.

1. Exponential function (preferred) 
Works by applying exponential functions to N and r between t=0 and the user specified "start of balanced evolution", with continuous derivatives.

    y = a*exp(b*x^c)
    dy/dx = abcx^(c-1)exp(b*x^c)

    so a = y_0
       c = (xdy/dx)/(log(y_1/y_0)*y_1)
       b = (dy/dx)/(c*y_1*x_1^(c-1))

Which we hence solve at intervals throughout pre-core collapse.

2. Cubic function (if y_1 < y_0 but dy_1/dx > 0 (or vice-versa).
Works by applying an cubic function with smoothly varying values and derivatives at the end points. This is used when the exponential function would fail, and does not accurately reproduce the behaviour of a cluster. However, the essential properties - smoothly variant value and derivate - are maintained.

   y = ax^3 + bx^2 + d
   dy/dx = 3ax^2 + 2bx 

   so d = y_0
      a = (2.0/x**3)*(y_0+(dydx*x)/2.0-y_1)
      b = (dydx-3.0*a*tcc**2)/(2*tcc)

which are alternatively used. */

void node::core_collapse(){

  double trh_old = trh();                     //temporary relaxation time store

  /* array 'y' stores initial values (y_0), 'dr1' stores final values (y_1),
     dr2 stores final derivatives (dy_1/dx_1), dr3 stores b and dr4 stores c. */
   
  collapse_time = toff*trh();                //Works out post collapse "targets"

  trial[0] = dr1[0] = N*init.f_N;
  trial[1] = dr1[1] = r*init.f_r;
  trial[2] = dr1[2] = mm*init.f_mm;
  l2 = true;  
  solve_diffs(dr2);                           //Gets post collapse derivatives
  for (int i = 0 ; i < 3 ; i++){              //Decides if polynomial
    if ((dr1[i] < 1.15*y[i] && dr2[i] > 0 ) || (dr1[i] > y[i]/1.15 && dr2[i] < 0 ))
    //factors of 1.15 - prevent "very" sharp turns in exp functions.
	dr6[i] = 1;
    else dr6[i] = 0;
  } 

  trial[0] = N;                               //Resets flag on trh+prepares
  trial[1] = r;
  trial[2] = mm;

  for (int i = 0 ; i < 3 ; i++){  
    if (dr6[i] == 0) {
      dr4[i] = exp_get_c(i);                  //Solves equations for c (exp)
      dr3[i] = exp_get_b(i);                  //Solves equations for b (exp)
    }
    else {
      dr3[i] = poly_get_a(i);                 //Solves equations for a (poly)
      dr4[i] = poly_get_b(i);                 //Solves equations for b (poly)
    }
  }
  
   while (t<collapse_time) {
    l1 = l2 = l3 = l4 = l5 = true;
    tstep = params.frac*trh();
    t += tstep;
    precc_evolve(trial);
    N = trial[0]; r = trial[1]; mm = trial[2];
    output();
    trh_old = trh();                           //Stores previous relaxation time
    n_relax += (2.0*tstep)/(trh_old+trh());    //Count of elapsed trh
    }

  trial[0] = y[0] = Nstart = N;                //Sets for further evolution
  trial[1] = y[1] = r;
  trial[2] = y[2] = mm;
  trial[3] = y[3] = n_relax;
  interp = 0;                                  //End of interpolation
}   
   
double node::exp_get_c(int i){
  return dr2[i]*collapse_time/(dr1[i]*log(dr1[i]/y[i]));
}

double node::exp_get_b(int i){
  return dr2[i]/(dr4[i]*dr1[i]*pow((double)collapse_time,(double)(dr4[i]-1.0)));
}

double node::poly_get_a(int i){
  return (2.0/pow(collapse_time,3))*(y[i]+(dr2[i]*collapse_time)/2.0-dr1[i]);
}

double node::poly_get_b(int i){
  return (dr2[i]-3.0*dr3[i]*pow(collapse_time,2))/(2.0*collapse_time);
}

void node::precc_evolve(double out[]){
  for (int i = 0 ; i < 3 ; i++) {
    if ( dr6[i] == 1 ) out[i] = (dr3[i]*pow(t,3))+(dr4[i]*pow(t,2))+y[i];
    else out[i] = y[i]*exp(dr3[i]*pow((double)t,(double)dr4[i]));
    if (out[i] != out[i]) out[i] = y[i];       //NaN check
  }
}
