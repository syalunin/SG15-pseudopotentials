
# SG15 ONCV pseudopotentials with core states

This repository contains the SG15 Optimised Norm-Conserving Vanderbilt (ONCV) pseudopotentials, generated using a modified version of the [scalar-relativistic code](http://www.quantum-simulation.org/potentials/sg15_oncv/) by D. R. Hamann. I modified the code to ensure that the pseudopotentials now include the wavefunctions of the core states for all atoms. Except for this, these pseudopotentials are identical to those obtained by M. Schlipf and F. Gigi, now known as the SG15 collection.

The core wavefunctions are written in a separate section of the pseudopotential file in GIPAW format, enabling Quantum ESPRESSO to read and use them in calculations. These wavefunctions can theoretically be used to calculate NMR shielding tensors, EFG tensors, and other magnetic resonance properties using tools such as the [QE-GIPAW package](https://github.com/dceresoli/qe-gipaw). However, full compatibility has not yet been tested.

### Installation
comming soon...

### Building
To generate the modified code, run the 'make' command in the root directory.

### Contact

If you have any questions, please feel free to contact me via email syalunin@users.noreply.github.com
