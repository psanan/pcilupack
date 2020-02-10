# Example ILUPACK installation instructions

2020.02.08

Explicit steps to obtain and test ILUPACK on a "standard" system,
here Ubuntu 18.04.3 LTS with gcc compilers, using the "GNU64" platform
for ILUPACK. See the ILUPACK documentation for general instructions.

1. Visit the [ILUPACK website](http://www.icm.tu-bs.de/~bolle/ilupack/)

2. Click on "Download" and then the first link ("zip file") to obtain ILUPACK v2.4 (05102016 version)

3. Copy the downloaded file to the desired destination

    mv $HOME/Downloads/ilupack05102016.zip $HOME/code

4. Unzip

    cd $HOME/code
    unzip ilupack05102016.zip

5. Obtain MC21 and MC64 from HSL, by visiting the [HSL website](http://www.hsl.rl.ac.uk/),
navigating to the sub-pages for MC21 and MC64, and requesting the libraries under an
appropriate license. When approved, download `mc21-1.0.0.tar.gz` and `mc64-1.6.0.tar.gz`
and unpack them

    mv $HOME/Downloads/mc21-1.0.0.tar.gz $HOME/code
    mv $HOME/Downloads/mc64-1.6.0.tar.gz $HOME/code
    cd $HOME/code
    tar xzvf mc21-1.0.0.tar.gz
    tar xzvf mc64-1.6.0.tar.gz

6. Compile the required HSL objects inside ILUPACK's directory

    cp mc64-1.6.0/src/mc64*.f ilupack/notdistributed/
    cp mc21-1.0.0/src/mc21*.f ilupack/notdistributed/
    cd ilupack/notdistributed
    gfortran -O3 -fPIC -c -o MC21D.o mc21d.f
    gfortran -O3 -fPIC -c -o MC64D.o mc64d.f
    gfortran -O3 -fPIC -c -o MC21S.o mc21s.f
    gfortran -O3 -fPIC -c -o MC64S.o mc64s.f
    # quicker: for file in *.f; do stem=${file%.f}; gfortran -O3 -fPIC -c -o ${stem^^}.o $file; done

 6. Test ILUPACK by building a 32-bit integer, double, symmetric example

    cd ilupack/simple_examples
    make -f MC64/makefile.GNU64 MAIN=dmainsym
     ./dmainsym.out

 Output:
    factorization successful with 2 levels completed
    final elbow space factor=    2.01
    approximate solution after one preconditioning step
     -4.2725543538362838e-02
     -1.5897898097154955e-02
      1.1053236900891336e-01
     -1.2041337668369712e-01
      2.6074727770036638e-02
     -1.2240764087025115e-01
      2.0476595964287550e-01
      2.1187527186133781e-01

    iteration successful completed after 5 steps
    approximate solution after completing the iterative process
     -4.2725543538362859e-02
     -1.5897898097154934e-02
      1.1053236900891339e-01
     -1.2041337668369716e-01
      2.6074727770036631e-02
     -1.2240764087025119e-01
      2.0476595964287558e-01
      2.1187527186133778e-01

    residual norm      0.0
