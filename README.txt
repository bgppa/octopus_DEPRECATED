OCTOPUS is a structured collection of mathematical utilities arising
during my PhD studies.

It is intended to be an exercise to improve my programming skills
and check empirically how many theoretical results behave in practice.

The **core** problem is here explained.
Given:
         - a function G: R^n -> R^m
         - a (noised) observation y \in R^m
Find:
        - element x \in R^n such that G(x) = y ("almost")

In other words: find a preimage of a noised y for a function G.
The technique here developed is bayesian.
I.e.: 

        1. we start on a "believe" on x
        (can be a random value, if we have no idea)
        believe = some probability measure, called "prior"
        
        2. we walk around such a belive, take a new_point,
        and observe how close G(new_point) is to y;
        Accept if enough, refuse otherwise.
        decision process = MCMC with Gibbs potential

        3. at the end we obtain a *probability*
        distribution for possible x values.
        called: posterior distribution.

        4. a kmeans algorithm (i.e. multidimensional histograms)
        tells us what is the most probably value for x,
        "solving" the problem of G(x) = y.
        called: MAP
        (Maximum A Posteriori estimator)

The algorithm above is sometimes called "bayesian inversion",
or "bayesian reconstruction", or similar.
Note that it involves:
        - linear algebra operations (of course...);
        - some probability (MCMC,...);
        - machine learning (kmeans);
        - the ability to implement G, the operator.

The last point is the hardest: G can be taken, for example, as a ordered
outcome of a PDE. Therefore, to evaluate and simulate G (step required)
can be challenging and costly.

Therefore my organization:
a subfolder called "mylib" that contains **all except G**,
i.e. tools like Linear Algebra, Simulations processes, File IO, etc..
mylib does not require any dependence and can be smoothly compiled
by using the script "compile.sh"
It's an independend step which must work by ist own.

Remark: yes, I could use LAPACK. But not now; the priority is on others.

A series of folders, using mylib as basic toolkit, where various
operators G are defined and the bayesian inverse algorithm is tested.


----    COMPLETE HIERARCHY DESCRIPTION  ----

mylib:  (folder)        basic toolkit

fNAME:  (many folders)  temporary files for developing a new
                        *f*unction named NAME that will be added
                        to the library. Examples:
                        fMultidimensionalGaussian
                        fPoissonProcess;

pNAME:  (many folders)  _complete_ *p*rojects of application
                        of the Bayesian Inverse Theory.

dNAME:  (many folders)  projects of application of the Bayesian Inverse Theory,
                        in *d*evelopment or further testing.

tests:  (folder)        tests for checking a specific library funcionality


REMARKS:
        Every f folder contains source code that can be
        independent of the main library, since it is supposed to
        be a new funcionality to add (usually, it uses Linear Algebra Lib).

        Every p folder refers to a project showing the usage of
        Bayesian principles, for example "pConvectionODE"
        deals with inverting the convection ODE.

        Every p folder must be completely independent of any other p folder,
        as well as of all others f folders.

        Every p folder **must** only depend on mylib, and contains:
        - README.txt    that explains the purpose of the project;
        - main.c        the main C source file ready to be compiled by a...
        - makefile      written once and valid for every project.

        [A project's makefile is very trivial: it searches for the
        [already-compiled mylib libraries, then compile main.c
        [linking it to these libraries (located in subfolders of mylib)]

        - other .c or .h files are admitted, as well as modifications
          to the makefile, the only requirement is to be well-documented.

Gold rule: **compiling a project MUST leave untouched the content of mylib**

Speaking about mylib, it contains three subfolders and a compile.sh script.
        - src:          where all the .c files are stored
        - include:      where all the .h headers are stored
        - obj:          where the compiled .o objects are stored
        - compile.sh    compile all src/*.c, move them to obj/

Every file.c in src has a file.h counterpart in include,
and they follow dependencies documented file-by-file.
Now they are ready to be used for f or p folders!
A project's makefile will expect obj/ to be ready,
without re-compiling or modifying anything.

d folders follow the same precise rules of p folders.

Summing up, here a pictorially organization:

|       README.txt
|       mylib           /compile.sh
                        /src/*.c
                        /include/*.h
                        /obj/*.o    

|       pNAME or dNAME  /main.c
                        /makefile
                        /*.c
                        /*.h
                        /README.txt
|       fNAME           / /* code that not necessarely must follow a rule */
|       tests           /various_tests.c
