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

import ParticleIO
import MeasurementTools

def main():

    ### Load data ####

    myfile = sys.argv[1] #Expects .h5 input file
    traj = ParticleIO.load_traj(myfile) #Extract data
    eq_frac = 0.2 #cut off first 20% of data (equilibration)

    #### Compute allowed wavevectors ###

    #First get max magnitude and spacing.
    dq = 2*np.pi/np.max(traj['edges'])
    qmax = 2*np.pi #Assumes smallest relevant distance =1 for now
    #qvals = np.arange(-qmax, qmax+dq, dq)

    #Generate a grid of wavevectors (w/ lattice constant dq),
    #then select those within a sphere of radius qmax.
    qvals = np.arange(0, qmax+dq, dq)
    qx, qy, qz = np.meshgrid(qvals, qvals, qvals)
    qlist = []
    for kx in range(qx.shape[0]):
        for ky in range(qy.shape[1]):
            for kz in range(qz.shape[2]):
                qvec = np.array([qx[kx,ky,kz], qy[kx,ky,kz], qz[kx,ky,kz]])
                if np.linalg.norm(qvec)<=qmax:
                    qlist.append(qvec)

    #### Compute S(q) for each wavevector ####
    sqlist = []
    for q in qlist:
        print(q)
        sqlist.append(get_sq(traj['pos'], traj['edges'], q, eq_frac))

    #### Output S(q) to file ####
    print(sqlist)
    return 0

@numba.jit(nopython=True)
def get_sq(pos, edges, q, eq_frac):

    sq = 0. + 0.j
    N = pos.shape[1]
    traj_len = pos.shape[0]
    eq_len = int(eq_frac*traj_len)
    for t in range(eq_len, traj_len):
        print(t)
        for i in range(N):
            for j in range(N):
                sq += np.exp(-1j*np.dot(q, MeasurementTools.get_min_disp(pos[t,i,:], pos[t,j,:], edges)))
    sq *= (1.0/(N*(traj_len-eq_len)))

    return sq

if __name__ == '__main__':
    main()