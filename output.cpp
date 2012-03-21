/*---------------------------------------------------------------------------*/
//Output functions of EMACSS
#include "emacss.h"

void node::output(){
  if (first) {
    /* Sets initial output properties - prints to cerr so output can be 
       separately piped into a file*/
    
    double rj = pow((init.N0*init.mm0)/(3.0*init.galaxy.M),  
		   (1.0/3.0))*init.galaxy.R;  //Jacobi radius at first output
    
    cerr << endl;                    
    cerr << "Initial cluster properties ";
    if (init.units == 0) cerr << "using N-body units:" << endl;
    if (init.units == 1) cerr << "using real units:" << endl;
    
    cerr << "N\tr\trh\trj\tM_gal\tr_gal\tmstar\tf_N\tf_r\tf_ms\tzeta" << endl; 
    //Labels
   
    cerr <<setprecision(6)<<init.N0<<'\t';
    if (init.units == 0) cerr <<setprecision(3)<<1.0<<'\t';
    if (init.units == 1) cerr <<setprecision(3)<<init.r0<<'\t';
    cerr <<setprecision(3)<<init.r0*0.77<<'\t';
    cerr <<setprecision(3)<<rj<<'\t';
    cerr <<setprecision(2)<<init.galaxy.M/init.mm0<<'\t';
    cerr <<setprecision(5)<<init.galaxy.R<<'\t';
    cerr <<setprecision(2)<<init.mm0<<'\t';
    cerr <<setprecision(3)<<init.f_N<<'\t';
    cerr <<setprecision(3)<<init.f_r<<'\t';
    cerr <<setprecision(3)<<init.f_mm<<'\t';
    cerr <<setprecision(3)<<init.zeta<<endl;  
    //Values input
    
    if (init.units == 1){            //Real Units
      cerr << "-\t[pc]\t[pc]\t[pc]\t[M_sun]\t[pc]\t[M_sun]\t-\t-\t-\t-" << endl;
    }
    
    else {                           //N-body units
      cerr << "-\t[NBODY]\t[NBODY]\t[NBODY]\t[NBODY]\t[NBODY]\t-\t-\t-\t-\t-" << endl;
    }
    
    cerr << endl; 
    cerr << "Output takes the form:" << endl;
    
    if (init.units == 1){             //Real Units  
      cerr << "1\t\t 2\t3\t4\t5\t6\t7\t\t8\t\t9\t\t10\t\t"; //Labels  
      cerr << "11\t\t12\t13\t14\t15\t\t16\t\t17\t18\t19\t20\t\t21\t\t22" << endl;
      cerr << "t\t\t N\tM_clus\tr\trj\tm_bar\tt_rh\t\tE\t\tdNdt\t\tdrdt\t";
      cerr << "\tt\t\tM_clus\tr\trj\tm_bar\t\tt_rh\t\tE\tzeta";
      cerr << "\tn_relax\tdNdt\t\tdrdt\t\tPre_Collapse?" << endl;
      cerr << "[Myr]\t\t -\t[M_sun]\t[pc]\t[pc]\t[M_sun]\t[Myr]\t\t[M_sun(pc/"; 
      cerr << "\t[/Myr]\t\t[pc/Myr]\t[NBODY]\t\t[NBODY]\t[NBODY]\t[NBODY]";
      cerr << "\t[NBODY]\t\t[NBODY]\t\t[NBODY]\t-\t-\t[NBODY]\t\t[NBODY]";
      cerr << "\t\t-" << endl;
      cerr <<"\t\t\t\t\t\t\t\t\tMyr)^2]" << endl;
    }
    else {                            //N-body units
      cerr << "1\t\t 2\t3\t4\t5\t6\t\t7\t\t8\t9\t10\t11\t\t12\t\t13" << endl;
      cerr << "t\t\t N\tM_clus\tr\trj\tm_bar\t\tt_rh\t\tE";     
      cerr << "\tzeta\tn_relax\tdNdt\t\tdrdt\t\tPre-collapse?" << endl;
      cerr << "[NBODY]\t\t -\t[NBODY]\t[NBODY]\t[NBODY]\t[NBODY]\t\t[NBODY]\t";
      cerr <<"\t[NBODY]\t-\t-\t[NBODY]\t\t[NBODY]\t\t-" << endl;
    }
    
    first = false;
  }
  
  /* Output - specifiers for particular formatting */
  if (init.units == 1){                //Real Units 
    printf("%6.3e\t%6.0f\t%5.0f\t%5.3f\t%5.2f\t%5.3f\t%6.3e\t%6.3e\t%6.3e\t%6.3e\t%6.3e\t%5.3f\t%5.3f\t%5.2f\t%6.3e\t%6.3e\t%5.4f\t%5.3f\t%5.2f\t%6.3e\t%6.3e\t%d\n",t*T_star,N,N*M_star/init.N0,r*R_star,r_jacobi()*R_star, mm*M_star,trh()*T_star,energy()*pow(M_star,2)*G_star/R_star,dNdt()/T_star,drdt()*R_star/T_star,t,N/init.N0,r,r_jacobi(),mm,trh(),energy(),zeta,n_relax,dNdt(),drdt(),interp);
  }
  
  else{                                 //N-body units
    printf("%6.3e\t%6.0f\t%4.4f\t%5.3f\t%5.2f\t%6.3e\t%3.2e\t%5.4f\t%5.5f\t%5.2f\t%6.3e\t%6.3e\t%d\n",t,N,N/init.N0,r,r_jacobi(),mm,trh(),energy(),zeta,n_relax,dNdt(),drdt(),interp);
  }
}
