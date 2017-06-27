ising-2d-Gaussian-J
===================
The codes starting with "E"  basically computes the energies versus 1/temperature : We take 2 replicas: EmA-vs-beta.cc computes the enrgy for the constrained part of the replica, say region A. When you make a Monte Carlo move, the energy difference is only given by the next nearest neighbours. Thta is incorporated in the code to reduce run-time. This has to be modified for nnn case.


EmB-vs-beta.cc computesthe same but for the unconstrained region B of the replica

disorder-gauss.cc generates a random no. from a given gaussian distribution. This we use as the random coupling.

Run using bash script mutualinfo.sh. For example\*

`./mutualinfo.sh <beta_min> <beta_max> <del_beta> `

\* *You might need to first make the script executable, i.e. chmod +x mutualinfo.sh*

