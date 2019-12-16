/* In this specific example we have locally defined the operator G
 * acting like a Neural Network (not deep). Further information later */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include "myblas.h"
#include "ranvar.h"
#include "pcninv.h"
#define SECONDLAYER 30

/* This is the Neural Network example, which goal is the estimation
 * of its weights/biases by looking at its output.
 * The NN goes from R^z to {0,1}, the dimension z needs to be fixed
 * as a global variable in this temporarely script.
 * Given the NN we need to define the observation operator G
 * in a way to suit the Bayesian Inverse Problem.
 * Then remember: the observation operator must go from the set 
 * of parameters - weights/bias in our case, to the classification's value
 * of some set of points: we need so to specify them.
 * I did something easy: the user specify how many of them, say q.
 * Then q (uniform) random points in [-20,20]^z are generated one
 * time and stored in a global array. They will be classified!
 * That said, how does the G operator precisely work?
 * for each x_i beloging to the array above, it produces 0 or 1 according to:
 * weights -> NN(weights, x_i) -> {0,1}
 * Notice that so, when we have a bunch of 0, 1nes as output, inverting the
 * G operator allows to estimate the NN weights.
 * One you have its estimation, 1 - the residual error is basically
 * its precision rate. */

/* Comments about the network structure:
 * first layer: dimension z;
 * the second layer contains z * SECONDLAYER nodes
 * the third is the output
 * Therefore, once z is specified, can automatically deduce
 * all the others */

/* This is the z in the comments above, i.e. lenghts of input layer */
int glb_cube_dimension = 1;
/* how many weights_and_biases? */
int glb_w_and_b_dim = 0; /* Initialize later, in main */
double *glb_observed_pt = NULL; /* random points that will be classified,
                                i.e. the points on which G makes observation */

void init_glb_observed_pt(int how_many)
{
        /* There are how_many points to initialize, each of dim
         * glb_cube_dimension */
        glb_observed_pt = malloc(sizeof(double) * glb_cube_dimension
                                                        * how_many);
       for (int i = 0; i < how_many; ++i) {
                /* Initialize the i-th point as a multim uniform in -20,20 */
                for (int j = 0; j < glb_cube_dimension; ++j) {
                        glb_observed_pt[i * glb_cube_dimension + j] =
                                rndmUniformIn(-20., 20., NULL);
                }
        } 

        printf("Points to classify have been initialized:\n");
        for (int i = 0; i < how_many; ++i) {
                printVec(glb_observed_pt + i, glb_cube_dimension);
        }
}


/* Following functions auxhiliary for the Neural Network */
double heaviside(double t){
    return (t > 0 ? 1. : 0.);
}

/* Logistic is the sigmoid function used in the Watanabe text, too */
double logistic(double t){
    return 1./( 1. + exp(-t));
}

double ReLU(double t){ 
    return t > 0 ? t : 0.;
}

double identity(double t){
    return t;
}

/* Takes weights, their cardinality, a point, its dimension
 * and evaluate a simple NN in such a point */
double NN_model(const double *w_and_b, int len, double *z, int zdim)
{
        /* This function represents a non-deep Neural Network,
         * just two layers but, eh eh, good layers
         * For more about the notation you can use the slides from Jannik */
        int i, j;
        /* The dimension of the input layer coincides with the input dimension*/
        int dil = zdim;
        /* Just for simplicity, as a toy model, we set the hiddel layer
         * dimension as k times the input layer. CAN BE CUSTOMIZED, clearly */
        int k = SECONDLAYER; /* <- completely customizable constant */
        /* Dimension Hidden Layer */ 
        int dhl = k*zdim;
        /* Set the activation functions as logistics */
        double (*phi) (double);
        double (*ending_phi) (double);
//      phi = logistic;
        phi = ReLU;
        ending_phi = heaviside;

        /* Let's take into account some variables useful for counting the
         * dimensions;
         * the number of weight = dil * dhl + dhl
         * number of biases = dhl + 1
            So: w_and_b MUST have dimension equal to:
     dil*dhl + dhl + dhl + 1 , i.e.
    k*z*z + 2*k*z + 1
    */
        assert(len == (k * zdim * zdim + 2 * k * zdim + 1));
        int M = dil*dhl + dhl; /* Number of the weights */
        int bb = dhl + 1;   /* Number of the biases */
        double sum1 = 0;
        double sum2 = 0;

    /* Since all the weight w and biases b are included in the same variables w_and_b,
     * we spend some words to explain how they are represented
     * w_and_b[t] corresponds to...
     ... b(3) if t = bb-1 (so the last element in the array)
     ... TO COMPLETE, but I wrote a complete explanation on my notebook */

    for (j=0; j<dhl; ++j) {
        sum1 = 0;
        for (i = 0; i<dil; ++i) {
            sum1 += w_and_b[i * dhl + j] * z[i] + w_and_b[M + i];
        }
        sum2 += w_and_b[dil*dhl + j] * phi(sum1) + w_and_b[bb-1];
    }

    return ending_phi(sum2);
}



/* Let's define now the observation operator G */
void G(const double *wb, int wb_dim, double *y, int codomain)
{
        /* Simple: for each point in glb_observed,
                perform a NN evaluation by using the given weigths */
        for (int i = 0; i < codomain; ++i) {
                y[i] = NN_model(wb, wb_dim, glb_observed_pt + i, 
                                                glb_cube_dimension);
        }
}


/* Produce toy-model data. More precisely, it is assumed to
 * have G from R^domain_dim to R^codomain_dim.
 * The array x is initialized with random data, uniformly
 * between -10 and 10 (arbitrarely choice, no meaning).
 * Then y is created by applying G on x and is then perturbed by a noise.
 * So the aim of the main script will is to re-compute x having only y.
 * Since the true values of x are known, true error can be computed.
 * Parameters:
 - noise: covariance of gaussian's noise;
 - x: array of dimension domain_dim, elements set randomly;
 - y: array of dimension codomain_dim; elements set as observations. */
void createToyData(double noise, double *x, int domain_dim,
                        double *y, double *noise_free_y, int codomain_dim)
{
        int i = 0;
        /* Randomize the parameters x */
        for (i = 0; i < domain_dim; ++i) {
                x[i] = rndmUniformIn(-2, 2, NULL);
        }
        G((const double *) x, domain_dim, noise_free_y, codomain_dim);
        /* Put a noise on each value of y */
        for (i = 0; i < codomain_dim; ++i) {
                y[i] = noise_free_y[i] + rndmGaussian(0, noise, NULL);
        }
}

int main(int argc, char *argv[]) {
        /* Setup the number of weigths and biases, a global value,
         * according to the number of nodes and input dimension */
        glb_w_and_b_dim = SECONDLAYER * glb_cube_dimension * glb_cube_dimension
                                + 2 * SECONDLAYER * glb_cube_dimension + 1;
        srand(time(NULL));
        /* Noise used to produce the toy models data;
         * Noise introduced in the MCMC algorithm */
        double data_noise = 1e-2; 
        double mcmc_noise = 1e-2;

        /* The algorithm is very sensitive to the number of
         * produced samples, and how many monte carlo cycles
         * are used to produce each of it.
         * Default values: 2^10, 2^12 (powers set later) */
        int n = 10;
        int mcmc = 12;

        /* Default value for domain and codomain of G */
        int domain_dim = glb_w_and_b_dim; /* do not touch it */
        int num_observations = 50; /* <- THIS CAN BE FREELY SET, by coommand
                                        line, too */

        /* The values above can be modified via command arguments */
        if (argc >= 3){
                n = atoi(argv[1]);
                mcmc = atoi(argv[2]);
                if (argc == 5){
                /* Then also domain_dim and num_observations */
                        glb_cube_dimension = atoi(argv[3]);
                        glb_w_and_b_dim = SECONDLAYER * glb_cube_dimension
                                          * glb_cube_dimension + 
                                          2*SECONDLAYER*glb_cube_dimension + 1;
                        domain_dim = glb_w_and_b_dim;
                        num_observations = atoi(argv[4]);
                }
        }

        printf("Classifying %d random points belonging to the %d-dim\
                cube [-20, 20]. w and b required: %d\n",
                num_observations, glb_cube_dimension, glb_w_and_b_dim);
        getchar();

        /* Generate the random points we are going to classify */
        init_glb_observed_pt(num_observations);

        n = (int) pow(2, n);
        mcmc = (int) pow(2, mcmc);

        double *true_params = malloc(sizeof(double) * domain_dim);
        double *observed = malloc(sizeof(double) * num_observations);
        double *no_noise_observed = malloc(sizeof(double) * num_observations);
        assert(true_params != NULL && observed != NULL && no_noise_observed);

        createToyData(data_noise, true_params, domain_dim,
                        observed, no_noise_observed, num_observations);
        printf("** true coeff: \n");
        printVec(true_params, domain_dim);
        printf("\n** classification: \n");
        printVec(no_noise_observed, num_observations);   
        printf("\n");
        getchar();

        /* Now that the data are ready, set the bayes parameters */
        /* Output file where to write the posterior distribution */
        FILE *pfile = fopen("posterior.txt", "w");
        assert(pfile != NULL);
        FILE *ofile = fopen("Gposterior.txt", "w");
        assert(ofile != NULL);
        
        /* Covariance matrix for the gaussian */
        double *cov = malloc(sizeof(double) * domain_dim * domain_dim);
        /* Starting point where to start the chain */
        double *start = malloc(sizeof(double) * domain_dim);
        assert(cov != NULL && start != NULL);
        /* Set a random starting point, a small covariance matrix */
        for (int i = 0; i < domain_dim; ++i){
                start[i] = rndmUniformIn(-5., 5., NULL);
                for (int j = 0; j < domain_dim; ++j){
                        cov[i + j * domain_dim] = (i == j) ? 0.9 : 0.1;
                }
        }
        printf("** Starting point:\n");
        printVec(start, domain_dim);
        printf("\n%d samples, %d iterations per sample\n", n, mcmc);
        printf("--- press a key to continue ---\n");
        getchar();

        /* Create the seed for the parallelization */
        unsigned int *seed_r = malloc(sizeof(unsigned int) * n);
        seed_r[0] = time(NULL);
        for (int i = 1; i < n; ++i){
                seed_r[i] = seed_r[i-1] + 1;
        }

        /* Proceed with the bayesian inversion:
         * n = number of posterior samples to produces;
         * mcmc = number of monte carlo iterations;
         * true_params = true known parameters (toy model data)
         * G = the linear interpolation defined above
         * observed = vector containing the y_i
         * domain_dim = domain's dimension
         * observed = codomain's dimension
         * noise during the mcmc chain = mcmc_noise
         * 0.2 = beta, parameter for the pCN algorithm 0<beta<1
         * cov = my covariance matrix, prior gaussian
         * start = starting point for the chain
         * pfile = file to write posterior distribution (values, probabilities)
         * ofile
         * true G(weight) of the original weights
         * seed_r : seeds for paralelization MC
         * 0 = no verbose/debug mode */
        NNbayInv(n, mcmc, true_params, G, observed,
               domain_dim, num_observations,
               mcmc_noise, 0.2, cov, start, pfile, ofile,
               no_noise_observed, seed_r, 0);

        /* Free all the allocated memory */
        free(no_noise_observed);
        free(true_params);
        free(observed);
        free(cov);
        free(start);
        free(glb_observed_pt);
        free(seed_r);
        fclose(pfile);
        fclose(ofile);
        return 0;
}
