#!/usr/bin/env python3

# script to generate batch files with all the various
# parameter combinations


import numpy as np

import collections
import itertools

the_lambda = [ 5000 ]

# background mortality
d =  [[ 1, 1]]
gamma =  [[ 0, 0 ]]
p = [ 0.5 ] #list(np.arange(0,1,0.1))
rho =  [ 0.9 ]
phi = [[ 1.0, 0.01, 0.01, 1.0], [ 1.0, 0.25, 0.25, 1.0 ], [ 1.0, 0.5, 0.5, 1.0 ]]


#mu_a =  list(np.arange(0.001,0.05,0.05/50))

mu_h = [ 0.01 ]
mu_a = [ 0.01, 0.05 ]

mu_b =  [ 0.01 ]
#sdmu = list(np.arange(0,0.1,0.1/50))
#sdmu = [ 0.01, 0.05, 0.1 ]
sdmu_a = [ 0.1 ]  
sdmu_b = [ 0.1 ]  
sdmu_h = [ 0.1 ]  


inits =  [[ 1000, 2000 ]]
initi =  [[ 1000, 0 ]]
kappa =  [0.3]
U =  [5]
#sigma_e =  list(np.arange(5/50,5,5/50))
sigma_e = [ 5.0 ] 
tau2 = [ 2.0 ]

init_v = 1.0

# maximum number of timesteps simulation should run
t_max =  10**8

nrep = range(0,2)

ctr = 0

# whether it should run in the background
background_job = True

exe = "./xhost_specialization_plasticity"

bg_string = ""

if background_job:
    bg_string += " &"

for lambda_i in the_lambda:
    for d_i in d:
        for gamma_i in gamma:
            for p_i in p:
                for rho_i in rho:
                    for phi_i in phi:
                        for mu_a_i in mu_a:
                            for mu_b_i in mu_b:
                                for mu_h_i in mu_h:
                                    for sdmu_a_i in sdmu_a:
                                        for sdmu_b_i in sdmu_b:
                                            for sdmu_h_i in sdmu_h:
                                                for inits_i in inits:
                                                    for initi_i in initi:
                                                        for kappa_i in kappa:
                                                            for U_i in U:
                                                                for sigma_e_i in sigma_e:
                                                                    for nrep_i in nrep:
                                                                        for tau2_i in tau2:

                                                                            ctr += 1
                                                                            print("echo " + str(ctr))

                                                                            print(exe + " "
                                                                                    + str(lambda_i) + " "
                                                                                    + str(d_i[0]) + " "
                                                                                    + str(d_i[1]) + " "
                                                                                    + str(gamma_i[0]) + " "
                                                                                    + str(gamma_i[1]) + " "
                                                                                    + "\t"
                                                                                    + str(p_i) + " "
                                                                                    + str(rho_i) + " "
                                                                                    + "\t"
                                                                                    + str(phi_i[0]) + " "
                                                                                    + str(phi_i[1]) + " "
                                                                                    + str(phi_i[2]) + " "
                                                                                    + str(phi_i[3]) + " "
                                                                                    + "\t"
                                                                                    + str(mu_a_i) + " "
                                                                                    + str(mu_b_i) + " "
                                                                                    + str(mu_h_i) + " "
                                                                                    + "\t"
                                                                                    + str(sdmu_a_i) + " "
                                                                                    + str(sdmu_b_i) + " "
                                                                                    + str(sdmu_h_i) + " "
                                                                                    + "\t"
                                                                                    + str(init_v) + " "
                                                                                    + "\t"
                                                                                    + str(inits_i[0]) + " "
                                                                                    + str(inits_i[1]) + " "
                                                                                    + str(initi_i[0]) + " "
                                                                                    + str(initi_i[1]) + " "
                                                                                    + "\t"
                                                                                    + str(t_max) + " "
                                                                                    + str(kappa_i) + " "
                                                                                    + str(U_i) + " "
                                                                                    + "\t"
                                                                                    + str(sigma_e_i) + " "
                                                                                    + str(tau2_i) + " "
                                                                                    + str(bg_string) + " "
                                                                                    )
