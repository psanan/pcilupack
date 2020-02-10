# PETSc ILUPACK Plugin
Maintained by [Patrick Sanan](www.patricksanan.com)

This is a prototype of two plugins for [PETSc][1] which allows the use of preconditioners 
from [ILUPACK][2], by [Matthias Bollhoefer][3].

The first, `PCILUPACK`, wraps the multi-level ILU solver in ILUPACK.

The second, `PCILDL`, exposes a single-layer ILDL preconditioner,
with MC64-based symmetric weighted matching preprocessing.

## Acknowledgments

The makefile, install script, and registration code is adapted from 
[the DofColumns plugin][4] by [Tobin Isaac][5].

## Building and installation  

You need a working build of PETSc, configured with real, double precision scalars
and 32-bit integers (the defaults).
Visit [the PETSC website][1] for instructions.

You need a working ILUPACK library. See [the ILUPACK website][2].
Note that you need to provide your own compiled HSL object files
in `notdistributed/` (see example makefiles).

To build [and install] the plugin,
````
make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} PCILUPACK_OPENMP_FLAG=-fopenmp PCILUPACK_ILUPACK_DIR=/path/to/ilupack all [install]
````

To test that the shared library is working, try with a PETSc example, e.g. for `PCILUPACK`
````
$ cd $PETSC_DIR/src/ksp/ksp/examples/tutorials
$ make ex50
$ ./ex50 -dll_append /path/to/pcilupack/${PETSC_ARCH}/lib/libpcilupack.so -pc_type ilupack -ksp_monitor -ksp_view -pc_ilupack_droptol 1e-2 -pc_ilupack_droptolS 1e-2 -pc_ilupack_condest 10
factorization successful with 2 levels completed
final elbow space factor=    2.89
  0 KSP Residual norm 2.844097376456e-01
  1 KSP Residual norm 5.578652482730e-03
  2 KSP Residual norm 9.759582620006e-05
  3 KSP Residual norm 1.715185163440e-06
KSP Object: 1 MPI processes
  type: gmres
    GMRES: restart=30, using Classical (unmodified) Gram-Schmidt Orthogonalization with no iterative refinement
    GMRES: happy breakdown tolerance 1e-30
  maximum iterations=10000, initial guess is zero
  tolerances:  relative=1e-05, absolute=1e-50, divergence=10000.
  left preconditioning
  using PRECONDITIONED norm type for convergence test
PC Object: 1 MPI processes
  type: ilupack
    ILUPACK: droptol  : 0.01
    ILUPACK: droptolS : 0.01
    ILUPACK: condest  : 10.
  linear system matrix = precond matrix:
  Mat Object:   1 MPI processes
    type: seqaij
    rows=121, cols=121
    total: nonzeros=561, allocated nonzeros=561
    total number of mallocs used during MatSetValues calls =0
      has attached null space
      not using I-node routines
````

or for `PCILDL`

````
$ cd $PETSC_DIR/src/ksp/ksp/examples/tutorials
$ make ex50
$ ./ex50 -dll_append ~/pcilupack/arch-ubuntu-maint-double-debug/lib/libpcildl.so -pc_type ildl -pc_ildl_droptol 1e-3 -ksp_monitor -ksp_view
relative fill ILDL/A:  3.0e+00 (wrt 342 nz)
  0 KSP Residual norm 2.805274468165e-01
  1 KSP Residual norm 2.049842268559e-05
  2 KSP Residual norm 2.980826408179e-08
KSP Object: 1 MPI processes
  type: gmres
    GMRES: restart=30, using Classical (unmodified) Gram-Schmidt Orthogonalization with no iterative refinement
    GMRES: happy breakdown tolerance 1e-30
  maximum iterations=10000, initial guess is zero
  tolerances:  relative=1e-05, absolute=1e-50, divergence=10000.
  left preconditioning
  using PRECONDITIONED norm type for convergence test
PC Object: 1 MPI processes
  type: ildl
    ILDL: matching : 1
    ILDL: droptol : 0.001
    ILDL: ordering : metisn
  linear system matrix = precond matrix:
  Mat Object:   1 MPI processes
    type: seqaij
    rows=121, cols=121
    total: nonzeros=561, allocated nonzeros=561
    total number of mallocs used during MatSetValues calls =0
      has attached null space
      not using I-node routines

````
[1]: http://mcs.anl.gov/petsc
[2]: http://www.icm.tu-bs.de/~bolle/ilupack/
[3]: http://www.icm.tu-bs.de/~bolle/
[4]: https://github.com/tisaac/DofColumns
[5]: http://users.ices.utexas.edu/~tisaac/
