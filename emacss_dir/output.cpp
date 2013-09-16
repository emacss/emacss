#include "../emacss.h"

void node::output(dynamics dyn){
    
    static bool first = true;
    
    if (first){
	first = false;
	
	cerr << "Cluster initialised with:" << endl;
	if (units > 0){
	    cerr << "[R] N = " <<setprecision(6)<<N<< "  ";
	    cerr << "mm = " <<setprecision(4)<<mm*M_star<< " [M_sun]  ";
	    cerr << "r = " <<setprecision(4)<<rh*R_star<< " [pc]  ";
	    cerr << "rJ = " <<setprecision(4)<<rj*R_star<< " [pc]  ";
	    cerr << "zeta = " <<setprecision(3)<<E.zeta<< endl; 
	    cerr << "Galaxy Conditions:" << endl;	
	    cerr << "RG = " <<setprecision(3)<<galaxy.R*R_star/1e3<< " [kpc]  ";
	    cerr << "vG = " <<setprecision(3)<<galaxy.v*(R_star/T_star)/pcMyr << " [kms]  ";
	    cerr << "MG = " <<setprecision(3)<<galaxy.M*M_star<< " [M_sun]  ";
	    cerr << "type = " <<galaxy.type<< endl;
	}
	if (units != 1){
	    cerr << "[N-body] N = " <<setprecision(6)<<N<< "  ";
	    cerr << "mm = " <<setprecision(4)<<mm<< "  ";
	    cerr << "r = " <<setprecision(4)<<rh<< "  ";
	    cerr << "rJ = " <<setprecision(4)<<rj<< "  ";
	    cerr << "zeta = " <<setprecision(3)<<E.zeta<< endl;  
	    cerr << "Galaxy Conditions:" << endl;	
	    cerr << "[N-body] RG = " <<setprecision(3)<<galaxy.R<<"  ";
	    cerr << "vG = " <<setprecision(3)<<galaxy.v<< " ";
	    cerr << "MG = " <<setprecision(3)<<galaxy.M<< " ";
	    cerr << "type = " <<galaxy.type<< endl; 
	} 
	
	if (units > 0){
	    cerr << "Scaling factors:" << endl;	
	    cerr << "G_scale = " <<setprecision(4)<<G_star<< "  ";
	    cerr << "M_scale = " <<setprecision(4)<<M_star<< "  ";
	    cerr << "R_scale = " <<setprecision(4)<<R_star<< "  ";
	    cerr << "T_scale = " <<setprecision(4)<<T_star<< endl;
	}

    	cerr << "Other parameters:" << endl;
	cerr << "units = " <<units<< "  " << endl;;
	cerr << endl;
	
	if (units != 1){
	  fprintf(stderr,"  %-12s %-9s %-9s %-9s %-9s %-9s %-9s %-5s %-1s %-9s %-9s %-10s %-10s %-10s %-10s %-10s \n",
		  "(1)","(2)","(3)","(4)","(5)","(6)","(7)","(8)","(9)","(10)","(11)","(12)","(13)","(14)","(15)","(16)"); 
	  fprintf(stderr,"  %-12s %-9s %-9s %-9s %-9s %-9s %-9s %-6s %-1s  %-9s %-9s  %-10s %-10s %-10s %-10s %-10s  \n",
		  "t","N","M","rhoc","rc","rh","rj","kappa","S",
		  "t_rc","t_rh",
		  "lambda","xi","mu","epsilon","delta"); 
	}
	
	if (units > 0){
	  fprintf(stderr,"  %-12s %-9s %-9s %-9s %-9s %-9s %-9s %-5s %-1s %-9s %-9s  %-10s %-10s %-10s %-10s  %-10s \n",
		  "(1)","(2)","(3)","(4)","(5)","(6)","(7)","(8)","(9)","(10)","(11)","(12)","(13)","(14)","(15)","(16)"); 
	  fprintf(stderr,"  %-12s %-9s %-9s %-9s %-9s %-9s %-9s %-6s %-1s  %-9s %-9s  %-10s %-10s %-10s %-10s %-10s  \n",
		  "t","N","M","rhoc","rc","rh","rj","kappa","S",
		  "t_rc","t_rh",
		  "lambda","xi","mu","epsilon","delta"); 
	}
    }
    if (units  != 1){
      printf("%12.6e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %6.4f %1d  %9.3e %9.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	time, N, N*mm,rhoc(),rc, rh,rj,kappa,E.source,t_rc,t_rh,
	dyn.lambda(),dyn.xi(),dyn.mu(),dyn.epsilon(),dyn.delta());  
    }
    if (units  > 0){
      printf("%12.6e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %6.4f %1d  %9.3e %9.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	time*T_star, N, N*mm*M_star,rhoc()*(M_star/pow(R_star,3)),rc*R_star,
	rh*R_star,rj*R_star,kappa,E.source,t_rc*T_star,t_rh*T_star,
	dyn.lambda(),dyn.xi(),dyn.mu(),dyn.epsilon(),dyn.delta());
	}

    }
