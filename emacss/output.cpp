#include "../emacss.h"

void node::output(stellar_evo se,dynamics dyn){
    
    static bool first = true;
//    F = -se.sigma_R();
    if (first){
	first = false;
	
	cerr << "Cluster initialised with:" << endl;
	if (units > 0){
	    cerr << "[Real] N = " <<setprecision(6)<<N<< "  ";
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
	cerr << "vc = " <<setprecision(3)<<galaxy.vc.kms<< " [kms]  ";
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
	    printf("#1r t[Myr]\tN\t\tM[M_sun]\tr[pc]\t\trj[pc]\t\tRG[kpc]\t\tvG[km/s]\n");
	    printf("#2r trh[Myr]\tn_relax\t\tE[Real]\t\tMS-factor\tkappa\n");  
	}
	if (units != 1){
	    printf("#1n t\t\tN\t\tM\t\tr\trj\tRG\tvG\n");
	    printf("#2n trh\t\tn_relax\t\tE\t\tMS-factor\tkappa\n");  
	}
	if (out_time.Myr == numeric_limits<double>::max())
	    printf("#3  epsilon\txi\t\tmu\t\tgamma\tlambda\n\n");
    }
    cerr << pos[0]*R_star << ' ' << pos[1]*R_star << ' ' << vel[0]*V_star << ' ' << vel[1]*V_star << endl;
    if (out_time.Myr == numeric_limits<double>::max()){
	if (units > 0){
	    printf("#1r %8.4e\t%8.4e\t%8.4e\t%8.4e\t%8.4e\t%8.4e\t%8.4e\n",
		  time.Myr,N,N*mm.Msun,r.pc,rj.pc,galaxy.R.pc/1e3,galaxy.v.kms);
	    printf("#2r %8.4e\t%8.4e\t%8.4e\t%8.4e\t%8.4e\n",
		  t_relax.Myr,trhelapsed,E.real,mass_seg,kappa);
	}
	if (units != 1){
	    printf("#1n %8.4e\t%8.4e\t%8.4e\t%8.4e\t%8.4e\t%8.4e\t%8.4e\n",
		  time.nbody,N,N*mm.nbody,r.nbody,rj.nbody,galaxy.R.nbody,galaxy.v.nbody);
    	    printf("#2n %8.4e\t%8.4e\t%8.4e\t%8.4e\t\t%8.4e\n",
		  t_relax.nbody,trhelapsed,E.nbody,mass_seg,kappa);
	}
	printf("#3  %8.4e\t%8.4e\t%8.4e\t%8.4e\t%8.4e\n\n",
		se.epsilon(),dyn.xi(),dyn.mu(),dyn.gamma(),dyn.lambda());
    }
    if (time.nbody > out_time.nbody && units > 0){
	    printf("\n#1r  %8.4e\t%8.4e\t%8.4e\t%8.4e\t%8.4\n",
		  time.Myr,N,N*mm.Msun,r.pc,rj.pc);
	    printf("#2r  %8.4e\t%8.4e\t%8.4e\t%8.4e\t\t%8.4e\n",
		  t_relax.Myr,trhelapsed,E.real,mass_seg,kappa);
    }
    if (time.nbody > out_time.nbody && units != 1){
	    printf("\n#1n  %8.4e\t%8.4e\t%8.4e\t%8.4e\t%8.4\n",
		    time.nbody,N,N*mm.nbody,r.nbody,rj.nbody);
	    printf("#2n  %8.4e\t%8.4e\t%8.4e\t%8.4e\t\t%8.4e\n",
		  t_relax.nbody,trhelapsed,E.nbody,mass_seg,kappa);
    }
    if (time.nbody < out_time.nbody && N < 100 \
	&& units > 0 && out_time.nbody < numeric_limits<double>::max()){
	    printf("\n#1r  %8.4e\t0\t\t\t0\t\t0\t\n",out_time.Myr);
	    printf("#2r  0\t\t%8.4e\t0\t\t%8.4e\t0\n",trhelapsed,mass_seg);
    }
    if (time.nbody < out_time.nbody && N < 100 \
	&& units != 1 && out_time.nbody < numeric_limits<double>::max()){
	    printf("\n#1n  %8.4e\t0\t\t\t0\t\t0\t\n",out_time.nbody);
	    printf("#2n  0\t\t%8.4e\t0\t\t%8.4e\t0\n",trhelapsed,mass_seg);
    }
}
