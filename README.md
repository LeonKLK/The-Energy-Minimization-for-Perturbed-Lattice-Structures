# Introduction

This is a computational physics projects at The Chinese University of Hong Kong(CUHK). In this project I focused on lattice structures with face-centered cubic structure, including Aluminium and Magnesium Oxide(MgO). C is the main language of this project.

I focused on the Coulomb–Buckingham potential(MgO) and Lennard-Jones potential(Aluminium), in the `unpertubed structure` file you can find the construction for both substances.

In the `pertubed structure` file, I generated random numbers within a range of -0.05 to 0.05 in a unit of the respective lattice constant for each atom. We use conjugate gradient method and steepest decent method for the minimization of total energy(average energy per atom). Notice that if the pertubation is too large, the atoms may not go back to the unperturbated positions(average energy not converge). Also, if the lattice is not large enough(I use 5x5x5 for all cases), the convergence test would fail(or having weird values in energy).

In the subfiles you can see `performance.csv`, it shows the average energy per atom at each iteration.



# Remark
If you have another substances with other structures, remember to reconstruct the struture first, also the parematers for different potentials.

Those parameters of potential are found in online resources, you can use any references you want.

