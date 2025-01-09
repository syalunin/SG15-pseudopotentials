
# SG15 ONCV pseudopotentials with core states

This repository contains the SG15 Optimised Norm-Conserving Vanderbilt (ONCV) pseudopotentials, generated using a modified version of the scalar-relativistic code by D. R. Hamann. I modified the code to ensure that the pseudopotentials now include the wavefunctions of the core states for all atoms. Otherwise, these pseudopotentials are identical to those obtained by M. Schlipf and F. Gigi, now known as the SG15 collection.

The core wavefunctions are provided in a separate GIPAW section of the pseudopotential file. In addition to the wavefunctions themselves, this section also contains information about the orbital angular momentum and the principal quantum number of the orbitals, so that Quantum ESPRESSO can read and use them in calculations. Theoretically, these wavefunctions can be used to calculate NMR shielding tensors, EFG tensors, and other magnetic resonance properties using, for example, the [QE-GIPAW package](https://github.com/dceresoli/qe-gipaw). However, I have not tested all of these features, so use them at your own risk, without any guarantee from my side.

The original code can be accessed on [Dr. Hamann's official webpage](http://www.quantum-simulation.org/potentials/sg15_oncv/) and the procedure for obtaining the pseudopotentials is described in the works:

[D. R. Hamann, Phys. Rev. B 88, 085117 (2013)](http://link.aps.org/doi/10.1103/PhysRevB.88.085117)

[M. Schlipf and F. Gygi, Computer Physics Communications 196, 36 (2015)](https://www.sciencedirect.com/science/article/pii/S0010465515001897)

[P. Scherpelz, M. Govoni, I. Hamada, G. Galli, J. Chem. Theory Comput. (2016)](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00114)

### Installation
comming soon...

### Building
To generate the executables, run the 'make' command in the root directory.

If you have any questions, please feel free to contact me via email at syalunin@users.noreply.github.com
