ising-2d-Gaussian-J
===================
disorder-gauss.cc generates a random no. from a given gaussian distribution. This we use as the random coupling.

The codes starting with "E"  basically compute the energies versus 1/temperature : 

E-vs-beta-sim-anl.cc computes the energy for the normal/unreplicated system for a given temperature. When you make a Monte Carlo move, the energy difference is only given by the nearest neighbours. That is incorporated in the code to reduce run-time. This has to be modified for nnn case.

We take 2 replicas: EmA-vs-beta.cc computes the enrgy for the constrained part of the replica, say region A. 

EmB-vs-beta.cc computesthe same but for the unconstrained region B of the replica.

Mutual-info-vs-beta.cc computes the mututal information from the date files generated from runing the energy codes for each temperature.

The make file just runs all the codes together.

It is better to modify the energy codes a bit to store the spin configuration for the latest beta, as in: https://github.com/ipsita84/bc-model/blob/master/E-print-spin.cc


Run using bash script mutualinfo.sh. For example\*

`./mutualinfo.sh <beta_min> <beta_max> <del_beta> `

\* *You might need to first make the script executable, i.e. chmod +x mutualinfo.sh*

Sister directories for the same problem are at https://github.com/ipsita84/ising-2d-random-sign and https://github.com/dartonias/2DRandomIsing

