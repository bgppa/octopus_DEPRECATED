This is the folder containing all my library for numerical purposes.
Here the hierarchy is listed and explained.

mylib:		(folder)	the actual complete (in progress) library
extra_py:	(folder)	some python3 scripts useful for plotting / data analysis;
old:		(folder)	old files and examples, potentially no more useful;
fNAME:		(many folders)	temporary files for developing a new *f*unction named NAME that will be added to the library (examples: fMultidimensionalGaussian, fPoissonProcess);
pNAME:		(many folders)	_complete_ *p*rojects of application of the Bayesian Inverse Theory (examples: pOneDimensionalHeatEquation, pBayesOptimization).
dNAME:		(many folders)	projects of application of the Bayesian Inverse Theory, in _development_ or further testing (examples: tKuantum).
tests		(folder)	contains test that can be done for checking a specific library funcionality

Every f folder contains source code that can be independent of the main library, since it is supposed to be a new funcionality to add (usually, it uses Linear Algebra Lib).

Every p folder refers to a project showing the usage of Bayesian principles, for example "pConvectionODE" deals with inverting the convection ODE, while "wKuantum" is related to Karen's quantum chemistry master's thesis.

Every p folder must be completely independent of any other p folder, as well as of f folders.

Every p folder must only depend on mylib, and contains:
	- README.txt that explains the purpose of the project, and ist current status (developed, but how much?).
	- main.c the main C source file ready to be compiled by a...
	- makefile, written once and valid for every project.  A project's malefile search for the already-compiled mylib libraries, so compile main.c linking it to these libraries (located in subfolders of mylib);
	- other .c or .h files are admitted, as well as modifications to the makefile, the only requirement is to be well-documented. Furthermore, compiling a project MUST leave untouched the content of mylib.

mylib contains three subfolders and a compile.sh script.

They are src (where all the .c files are stored), include (where all the .h headers are stored), obj (where the compiled .o objects are stored).
Every file.c in src has a file.h counterpart in include, and they follow a hierarcy documented file-by-file.
Then, compile.sh compile all the .c and store them in obj/.
Now they are ready to be used for f or p folders!
A project's makefile will check there, without re-compiling or modifying anything, following the motto "divide-et-impera".

d folders follow the same precise rules of p folders.

Summing up, here a pictorially organization:

|	README.txt
|	mylib		|- compile.sh
			|- src - *.c
			|- include - *.h
			|- obj - *.o	

|	extra_py 	| *.py
|	pNAME or dNAME	|-main.c
			|-makefile
			|-*.c
			|-*.h
			|-README.txt
|	fNAME		| /* code that not necessarely must follow a rule */
|	tests		| various_tests.c
|	old		
