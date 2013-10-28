import emacss_module
import numpy as np

f = emacss_module

t_in = 0
t_out = 12000
tcc0 = 0
trhp0 = 0

mm = 0.64
vG = 220

def main():

    RG = np.logspace(-1,2,10)
    N = np.logspace(3,7,10)
    rh = np.logspace(-1,1.5,10)

    input = open('input.dat','w')
    output = open('output.dat','w')
    print >> input, "t", "N", "rh", "R", "mm", "MG", "RG", "vG", "zeta"
    print >> output, "t", "N", "rh", "rj", "RG", "vG"

    for rG in RG:
      for n in N:
        for r in rh:
            f.initialise(t_in,t_out,n,r,mm,rG,vG,tcc0,trhp0)
            t_o, N_out, M_out, r_out, rj_out, dNdt, dmdt, dMdt, drdt, tcc, trhp = f.sim()
            print >> input, t_in, n, r, mm, rG, vG, tcc, trhp
            print >> output, t_o, N_out, M_out, r_out, tcc, trhp
    input.close()
    output.close()
 
if __name__ =='__main__':  
    main()
