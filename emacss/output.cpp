#include "../emacss.h"

void node::output(stellar_evo se,dynamics dyn){
    
    static bool first = true;
    
    if (first){
	first = false;
	
	cerr << "Cluster initialised with:" << endl;
	if (units > 0){
	    cerr << "[R] N = " <<setprecision(6)<<N<< "  ";
	    cerr << "mm = " <<setprecision(4)<<mm.Msun<< " [M_sun]  ";
	    cerr << "r = " <<setprecision(4)<<r.pc<< " [pc]  ";
	    cerr << "rJ = " <<setprecision(4)<<rj.pc<< " [pc]  ";
	    cerr << "zeta = " <<setprecision(3)<<E.zeta<< endl;  
	}
	if (units != 1){
	    cerr << "[N-body] N = " <<setprecision(6)<<N<< "  ";
	    cerr << "mm = " <<setprecision(4)<<mm.nbody<< "  ";
	    cerr << "r = " <<setprecision(4)<<r.nbody<< "  ";
	    cerr << "rJ = " <<setprecision(4)<<rj.nbody<< "  ";
	    cerr << "zeta = " <<setprecision(3)<<E.zeta<< endl;  
	} 
	
	if (units > 0){
	    cerr << "Scaling factors:" << endl;	
	    cerr << "G_scale = " <<setprecision(4)<<G_star<< "  ";
	    cerr << "M_scale = " <<setprecision(4)<<M_star<< "  ";
	    cerr << "R_scale = " <<setprecision(4)<<R_star<< "  ";
	    cerr << "T_scale = " <<setprecision(4)<<T_star<< endl;
	}
	    
	cerr << "Galaxy Conditions:" << endl;	
	cerr << "RG = " <<setprecision(3)<<galaxy.R.pc/1e3<< " [kpc]  ";
	cerr << "vG = " <<setprecision(3)<<galaxy.v.kms<< " [kms]  ";
	cerr << "MG = " <<setprecision(3)<<galaxy.M.Msun<< " [M_sun]  ";
	cerr << "type = " <<galaxy.type<< endl;

    	cerr << "Other parameters:" << endl;
	cerr << "SE = " <<s<< "  ";
	cerr << "units = " <<units<< "  ";
	if (units > 0 && out_time.Myr < numeric_limits<double>::max())	
		cerr << "Output @ " <<out_time.Myr<<" [Myr]";
	if (units != 1 && out_time.Myr < numeric_limits<double>::max())	
		cerr << "Output @ " <<setprecision(4)<<out_time.nbody;
	cerr << endl;
	
        if (units > 0){
	    printf("#1r t[Myr]\tN\t\tM[M_sun]\tr[pc]\t\trj[pc]\n");
	    printf("#2r trh[Myr]\tn_relax\t\tE[Real]\t\tSource\t\tepsilon\n");  
	}
	if (units != 1){
	    printf("#1n t\t\tN\t\tM\t\tr\trj\n");
	    printf("#2n trh\t\tn_relax\t\tE\t\tSource\tepsilon\n");  
	}
	if (out_time.Myr == numeric_limits<double>::max())
	    printf("#3  xi\t\tmu\t\tgamma_se\tgamma_dyn\tlambda\n\n");
    }
    if (out_time.Myr == numeric_limits<double>::max()){
	if (units > 0){
	    printf("#1r %8.5e\t%8.5e\t%8.5e\t%8.5e\t%8.5e\n",
		  time.Myr,N,N*mm.Msun,r.pc,rj.pc);
	    printf("#2r %8.5e\t%8.5e\t%8.5e\t%1d\t\t%8.5e\n",
		  t_relax.Myr,trhelapsed,E.real,E.source,se.epsilon());
	}
	if (units != 1){
	    printf("#1n %8.5e\t%8.5e\t%8.5e\t%8.5e\t%8.5e\n",
		  time.nbody,N,N*mm.nbody,r.nbody,rj.nbody);
    	    printf("#2n %8.5e\t%8.5e\t%8.5e\t%1d\t\t%8.5e\n",
		  t_relax.nbody,trhelapsed,E.nbody,E.source,se.epsilon());
	}
	printf("#3  %8.5e\t%8.5e\t%8.5e\t%8.5e\t%8.5e\n\n",
		dyn.xi(),dyn.mu(),se.gamma_se(),dyn.gamma_dyn(),dyn.lambda());
    }
    if (time.nbody > out_time.nbody && units > 0){
	    printf("\n#1r  %8.5e\t%8d\t%8.5e\t%8.5e\t%8.5e\n",
		  time.Myr,N,N*mm.Msun,r.pc,rj.pc);
	    printf("#2r  %8.5e\t%8.5e\t%8.5e\t%1d\t\t%8.5e\n",
		  t_relax.Myr,trhelapsed,E.real,E.source,se.epsilon());
    }
    if (time.nbody > out_time.nbody && units != 1){
	    printf("\n#1n  %8.5e\t%8d\t%8.5e\t%8.5e\t%8.5e\n",
		    time.nbody,N,N*mm.nbody,r.nbody,rj.nbody);
	    printf("#2n  %8.5e\t%8.5e\t%8.5e\t%1d\t\t%8.5e\n",
		  t_relax.nbody,trhelapsed,E.nbody,E.source,se.epsilon());
    }
}
