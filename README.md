
# SG15 ONCV pseudopotentials with core electrons

This repository contains the SG15 Optimised Norm-Conserving Vanderbilt (ONCV) pseudopotentials, generated using a modified version of the scalar-relativistic code by D. R. Hamann. The [modified code](https://github.com/syalunin/SG15-pseudopotentials/tree/main/src) retains all the functionality of the [original code](http://www.quantum-simulation.org/potentials/sg15_oncv/) and can also store the wave functions of the core electrons in the pseudopotential file. Thus allowing the generation of Schlipf-Gigi pseudopotentials (SG15 collection) with extended capabilities.

The core wavefunctions are stored in a separate section of the pseudopotential file in GIPAW format, enabling Quantum ESPRESSO to read and use them in calculations. This feature has been tested. The core wavefunctions can also potentially be used to calculate NMR shielding tensors, EFG tensors, and other magnetic resonance properties using the [QE-GIPAW package](https://github.com/dceresoli/qe-gipaw). Full compatibility with this package, however, has not yet been tested.

### Prerequisites
Before compiling this project, make sure that the following tools are installed:
- Intel Fortran and C++ compilers
- Intel MKL

### Building
To generate the executables, run the `make` command in the root directory.

### Runing
The script `run.sh` generates pseudopotentials with core electrons for selected chemical elements, and `run_all.sh` does the same for all chemical elements.

### Contact
If you are interested in continuing this project, do not hesitate to contact me via email:
syalunin@users.noreply.github.com
