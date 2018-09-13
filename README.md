# MDACP (Molecular Dynamics code for Avogadro Challenge Project)

[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

## Summary

MDACP (Molecular Dynamics code for Avogadro Challenge Project) is an
efficient implementations of classical molecular dynamics (MD) method
for the Lennard-Jones particle systems.

## Files

- makefile.opt.fx10
  - make options for PRIMEHPC FX10
- mekfile.opt.mac
  - make options for Mac OS X

## Usage

	$ cp makefile.opt.mac makefile.opt
	$ make
	$ mpirun -np 1 ./mdacp

The latest information of MDACP is available at
http://mdacp.sourceforge.net/

## List of Developers

- Hiroshi Watanabe <hwatanabe_at_issp.u-tokyo.ac.jp>
  - Institute for Solid State Physics, University of Tokyo
  - (Corresponding Author)

- Masaru Suzuki
  - Department of Applied Quantum Physics and Nuclear Engineering, Faculty of Engineering, Kyushu University

- Nobuyasu Ito
  - Department of Applied Physics, School of Engineering, The University of Tokyo

## References

Details of algorithms used in MDACP are described in the following papars.

- Efficient Implementations of Molecular Dynamics Simulations for Lennard-Jones Systems
  - H. Watanabe, M. Suzuki, and N. Ito
  - Prog. Theor. Phys. 126 203-235 (2011)
  - [arXiv:1012.2677](https://arxiv.org/abs/1012.2677)
  - Describe Version 1.00 (flat MPI)
- Huge-scale Molecular Dynamics Simulation of Multibubble Nuclei
  - H. Watanabe, M. Suzuki, and N. Ito
  - Comput. Phys. Commun. 184 2775-2784 (2013)
  - [arXiv:1210.3450](https://arxiv.org/abs/1210.3450)
  - Describe Version 2.00 (MPI+OpenMP hybrid)

## Related Softwares

Here are the other MD programs which are OSS.

- [LAMMPS](https://lammps.sandia.gov/)
- [NAMD](https://www.ks.uiuc.edu/Research/namd/)

## History

- 2018-9-13 Ver 2.21: Fixd bug in Langevin dynamics.
- 2017-3-7 Ver 2.20: Added AVX2 implementations, contributed by [kohnakagawa](https://github.com/kohnakagawa).
- 2013-10-17 Ver 2.11: Fixed bugs in calculating observables.
- 2012-12-27 Ver 2.10:
  - Peformance for hybrid execution is improved.
  - Some examples of make options are included.
- 2012-10-15 Ver 2.00: MPI + OpenMP Hybrid parallelization.
- 2012-02-24 Ver 1.02: Fixed bug in calculating forces.
- 2011-01-15 Ver 1.01: Fixed bug in updating pairlists.
- 2010-12-12 Ver 1.00: First Release