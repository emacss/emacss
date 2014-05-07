#include "../emacss.h"

void node::output(stellar_evo se, dynamics dyn){
    if (units  == 1 && s == 1){
      printf("%12.6e %9.3e %9.3e %9.3e %9.3e %6.4f %4d %9.3e %9.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n",
	time.Myr, N, N*mm.Msun,r.pc,rj.pc,kappa,E.source,t_relax.Myr,
	dyn.lambda(),dyn.xi(),dyn.xi_ind(),dyn.xi()+dyn.xi_ind(),dyn.mu(),se.epsilon(),se.gamma_se(),dyn.gamma_dyn(),dyn.gamma());
    }
}
