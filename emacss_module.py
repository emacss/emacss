from ctypes import *
import numpy as np

lib = cdll.LoadLibrary('./emacss-code/libemacss.so')
InArray = c_double*9
OutArray = c_double*11
output = OutArray(0,0,0,0,0,0,0,0,0,0,0) #Initialise... for some reason

#Dummy setup, shouldn't be used apart from initialisation
test =  np.empty(9)

a = lib.emacss_new
a.restype = c_void_p

b = lib.input
b.argtypes = [c_void_p, POINTER(c_double)]
b.restype=None

c = lib.evaluate
c.argtypes = [c_void_p, POINTER(c_double)]
c.restype=None

class test_emacss(object):
    def __init__(self):
        self.emacss = a()

    def input(self,v):
        b(self.emacss,v)

    def evaluate(self,output):
        c(self.emacss,output)

f = test_emacss()

#Extracts data from sample file
def initialise(t,N,r,R,mm,MG,RG,vG,zeta):
    test[0] = t
    test[1] = N
    test[2] = r
    test[3] = R
    test[4] = mm
    test[5] = MG
    test[6] = RG
    test[7] = vG
    test[8] = zeta

def sim():
    v = InArray(test[0],test[1],test[2],test[3],test[4],test[5],test[6],\
                    test[7],test[8])
    f.input(v)   
    f.evaluate(output)
    found = [o for o in output]
    return found
