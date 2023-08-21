#Compute and output the static structure factor given an input trajectory
#The static structure factor is defined as:
#S(q) = (1/N) * \sum_{j,k} [exp(iq * (rj - rk))]
#where rj is the position of particle j
#and each sum is over all N particles.

#Input: Trajectory in h5md format (.h5)
#Output: S(q) vs allowed wavevectors q in ??? format

import numpy as np
import h5py
import sys
import numba

def main():

    myfile = sys.argv[1] #Expects .h5 input file

    pos, edges, times = load_traj(myfile) #Extract data

    #Compute allowed wavevectors

    #Compute S(q) for each wavevector

    #Output S(q) to file
    
    return 0

def load_traj(myfile):
    
    traj = h5py.File(myfile)

    return 0

main()