/* Trying to deduce the heat equation starting condition from the
 * solution at time 0.1 */
#include<Python.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<assert.h>
#include<math.h>
#include<wchar.h>
#include"BASIC.h"
#include"RANVAR.h"
#include"PCNINV.h"
#include"FILEIO.h"

double glb_time_observed = 0.01;

/* Supplementary global python variables useful for the script */
PyObject* pFunc = NULL;
PyObject* py_a = NULL;
PyObject* py_basis_expansion = NULL; 
PyObject* py_y = NULL; 
PyObject* py_obs_number = NULL;
PyObject* pArgs = NULL;

/* ------ KEY REMARK ----
 * the difference between this script and dHeatReverse relies
 * on __where__ the functions alpha, phi, and solver are defined.
 * In this case, they are stored in a python script!
 * So my program read them, run a C conversion and then uses
 * the bayesian inverse solver as usual.
 * Why all of that?
 * First of all, in order to check that something like that is
 * in principle possible. The very sad drawback is given of course
 * in terms of performance. Now The inversion requires 2 minutes
 * instead of 20 seconds. But it means that if, given a certain problem,
 * the user DOES NOT HAVE A C solver, he can use e.g. a pre-made python
 * library. PRactical example: Fenics for solving PDEs!
 * Slow, but better than nothing. 
*/


/* Recall that the bayesian inversion of an operator G:R^n->R^m
 * requires, of course, the possibility of evaluating such
 * an operator. G is defined in C as the function "solver",
 * taking in total four parameters: x input, its dimension,
 * y output (vector where to write the evaluation) and its
 * dimension. In this version, solver is actually loaded from a
 * python script. py_solver takes precisely the same parameters,
 * covert them into python equivalent and evaluate the
 * function "solver" which is supposed to be contained in the attached
 * python script (when loaded, is stored into pFunc) */

void py_solver(const double* a, int basis_expansion, double*y, int obs_number){
	/* Evoke the solver from the python script */
	py_a = PyTuple_New(basis_expansion);
	py_basis_expansion = PyLong_FromLong(basis_expansion);
	py_y = PyList_New(obs_number);
	py_obs_number = PyLong_FromLong(obs_number);

	int i=0;
	for(i=0; i<basis_expansion; ++i){
		PyTuple_SetItem(py_a, i, PyFloat_FromDouble(a[i]));
	}
	for(i=0; i<obs_number; ++i){
		PyList_SetItem(py_y, i, PyFloat_FromDouble(y[i]));
	}

	/* Ok, now all the objects have been converted to a Python equivalent */
	/* Modify the global variable pArgs in order to have the new arguments */
	PyTuple_SetItem(pArgs, 0, py_a);
	PyTuple_SetItem(pArgs, 1, py_basis_expansion);
	PyTuple_SetItem(pArgs, 2, py_y);
	PyTuple_SetItem(pArgs, 3, py_obs_number);

	PyObject_CallObject(pFunc, pArgs);

	/* Ok, now the python data have to be re-converted in C */
	/* More specifically: py_y has to be copied into y */
	for(i=0; i<obs_number; ++i){
		y[i] = PyFloat_AsDouble(PyList_GetItem(py_y, i));
	}
	/* TO IMPLEMENT: the py calling counter? */
}


/* (the solver user to create toy data is given from python; py_solver above */
/* The user wants to produce toy-model data in order to test the inversion's
 * effectiveness. This function:
 * takes a noise intensity;
 * randomizes an array of basis coefficients a, that can later be re-read
 * in order to quantify the actual error;
 * solves to the euqation and does the observations by using "solver" and
 * storing the results into y;
 * finally, y is noised according to the previous given intensity */
void createToyData(double noise,
		double* a, int basis_expansion, double*y, int obs_number){
	int i=0;
	/* Randomize the parameters a */
	for(i=0; i<basis_expansion; ++i){
		a[i] = rndmUniformIn(-10, 10);
	}

	/* Solve the equation via python solver storing the results in y */
	py_solver((const double*) a, basis_expansion, y, obs_number);
	printf("\n** noise-free obs: \n");
	printVec(y, obs_number);
	printf("\n");
	
	/* Put a noise on each value of y */
	for(i=0; i<obs_number; ++i){
		y[i] += rndmGaussian(0, noise);
	}

	/* The "result" of this function is given by the coefficients
	 * now stored in a, and the noised solution in y which can be used
	 * as observations for reconstruction purposes */
}	
	

	
	

int main(int argc, char* argv[]){
	Py_Initialize();

	/* Convert the main arguments to Python3-compatible parameters */
        wchar_t** _argv = PyMem_Malloc(sizeof(wchar_t*)*argc);
        for (int i=0; i<argc; i++) {
                 wchar_t* arg = Py_DecodeLocale(argv[i], NULL);
                 _argv[i] = arg;
        }
	PySys_SetArgv(argc, _argv);

	/* Load the module pySolver.py appearing in the current local directory */
	/* 1) convert the name into a Python string */
	PyObject* pName = PyUnicode_FromString("pySolver");
	/* 2) then import the module itself */
	PyObject* pModule  = PyImport_Import(pName);
	if(pModule == NULL){
		printf("Error: unable to import pySolver\n");
		return -1;
	}
	/* Ok, now the module as been imported */
	PyObject* pDictOfFunctions = PyModule_GetDict(pModule);
	pFunc = PyDict_GetItemString(pDictOfFunctions, "solver");
	if(PyCallable_Check(pFunc)){
		printf("Function solver successfully loaded\n");
	}
	else{
		printf("Error: unable to load solver from python\n");
		return -1;
	}

	/* Prepare the argument list for pFunc. It will be initialized during
	 * each call to py_solver, by copying the C arrays and integers */
	pArgs = PyTuple_New(4);

	/* Ok, now pFunc contains a python reference to the python solver */

	srand(time(NULL));
//	srand(2);

	double data_noise = 1e-4; /* Temporarely disabled in createToyData */
	double mcmc_noise = 1e-3;

	int n2 = 10;
	int mcmc2 = 12;
	int expansion = 3;
	int num_observations = 11;

	if(argc == 5 ){
		/* Then expansion and num_observations are given
		 * in input */
		expansion = atoi(argv[3]);
		num_observations = atoi(argv[4]);
	}

	if(argc >= 3){
		n2 = atoi(argv[1]);
		mcmc2 = atoi(argv[2]);
	}		
		
	/* Arguments? later */
	double* true_params = malloc(sizeof(double) * expansion);
	double* observed = malloc(sizeof(double) * num_observations);

	createToyData(data_noise, true_params, expansion,
			observed, num_observations);

	printf("** true coeff: \n");
	printVec(true_params, expansion);
	printf("\n** noised obs: \n");
	printVec(observed, num_observations);	
	printf("\n");

	/* Now that the data are ready, set the bayes parameters */
	/* Output file where to write the posterior distribution */
	FILE* pfile = fopen("posterior_measure.txt", "w");
	int n = (int) pow(2, n2);
	int mcmc = (int) pow(2, mcmc2);
	/* Residual error produced by the bayesian inversion */
	double err = 0;
	int i, j;
	
	/* Estimated parameters */
	double* map = malloc(sizeof(double)*expansion);
	/* Covariance matrix for the gaussian */
	double* cov = malloc(sizeof(double)*expansion*expansion);
	/* Starting point where to start the chain */
	double* start = malloc(sizeof(double)*expansion);
	if(map == NULL || cov == NULL || start == NULL){
		fprintf(stderr, "malloc failed\n");
		return -1;
	}

	/* Reset map, set a random starting point, a small covariance matrix */
	for(i=0; i<expansion; ++i){
		map[i] = 0;
		start[i] = rndmUniformIn(-10., 10.);
		for(j=0; j<expansion; ++j){
			cov[i + j*expansion] = (i==j)? 0.9: 0.1;
		}
	}

	printf("** Starting point:\n");
	printVec(start, expansion);
	printf("\n%d samples, %d iterations per sample\n", n, mcmc);


	/* Proceed with the bayesian inversion:
	 * n = number of posterior samples to produces;
	 * mcmc = number of monte carlo iterations;
	 * map = the vector which will contain the most frequent sample = solution = MAP
	 * NULL = the value of the true solution, a,b, is not known
	 * solver = the linear interpolation defined above
	 * observed = vector containing the y_i
	 * expansion = domain's dimension
	 * observed = codomain's dimension
	 * noise during the mcmc chain = mcmc_noise
	 * 0.2 = beta, parameter for the pCN algorithm 0<beta<1
	 * cov = my covariance matrix, prior gaussian
	 * start = starting point for the chain
	 * pfile = file will contain the posterior distribution (values, probabilities)
	 * 0 = no verbose/debug mode */
	err=bayInv(n, mcmc, map, true_params, py_solver, observed,
		expansion, num_observations, mcmc_noise, 0.2, cov, start, pfile, 0);

	/* err contains the residual error */
	/* Print the results */
	printf("MAP: ");
	printVec(map, expansion);
	printf("RES ERR: %.3f%%\n", err);
	printf("Observed output:\n");
	printVec(observed, num_observations);
	printf("MAP output :\n");
	py_solver(map, expansion, observed, num_observations);
	printVec(observed, num_observations);

	/* Free all the allocated memory */
	free(true_params);
	free(observed);
	free(map);
	free(cov);
	free(start);
	fclose(pfile);
	PyMem_Free(_argv);
	Py_Finalize();
	return 0;
}
