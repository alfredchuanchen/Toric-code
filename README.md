## Chern number for the lattice model
This is a project about studying the flux inside a Fermi liquid, which is based on dualization between a spin system and fermion-boson system.


The function of the files are the following:

- **parameterize.m**

establish the parameterization of the model, like (linear) system size, hopping amplitute, list of the $theta_i$'s conmbinations. After running it, it will
generate a *parateters.mat* file, which will be loaded in the consequent matlab codes.

- **H_BdG.m**

build up the BdG Hamiltonian, which depends on the position of the flux. In this model, we consider fermions has onsite energy $\Delta$, single-particle
hopping and pairing with amplitute $t$.

- **spectra.m**

calcualte the eigensystem of the $H_BdG$, in the form of u, v matrices, the eigenvalues are sorted in the form: ${E_1, E_2, ..., E_N, -E_1, -E_2, ..., -E_N}$,
with $E_1 > E_2 > ... > E_N >0$.

- **overlap.m**

calculate the overlap of the (normalized) ground states in two adjacent steps of the theta list, as for the final step with all $theta_i = \pi$, we calculate
$\langle \Omega_N| G_p | \Omega_1 \rangle$, with $G_p$ being the local gauge transformation at site (N/2+1, N/2).

- **berry_phase.m**

calculate the berry phase based on the output of *overlap.m* file.

- **pfaffian_LTL.m**

this is a code introduced by M. Wimmer [Efficient numerical computation of the Pfaffian for dense and banded skew-symmetric matrices](https://arxiv.org/abs/1102.3440), here I just slightly adapted it a bit: because we are sure that the matrices that we use is skew symmetric, so I commented
the command for checking whether the matrix is skew symmmetric or not.
*(One should note that this code is not perfect, as for some special cases, it can give rise to a wrong result, like 0, but it can be easily checked whether this happens or not)*

- **ticket_generator.nb**

this file generate the tickets for the slurm system, remember the convert the generated *.dat* file into *.sh* format before submit the job in slurm system.


