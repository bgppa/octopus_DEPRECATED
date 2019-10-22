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
        1. start on a "believe" on x
        (can be a random value, if have no idea)
        believe = some probability measure, called "prior"
        
        2. walk around such a belive, take a new_point,
        observe how close G(new_point) is to y;
        Accept if enough, refuse otherwise.
        decision process = MCMC with Gibbs potential

        3. at the end obtain a *probability*
        distribution for possible x values.
        called: posterior distribution.

        4. a kmeans algorithm (i.e. multidimensional histograms)
        tells what is the most probably value for x,
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

Therefore my motto: "Divide et impera":
A subfolder called "mylib" contains **all except G**,
i.e. tools like Linear Algebra, Simulations processes, File IO, etc..
It is a sort of exercise to recap some basic numerics, C, etc...
And I am happy to have re-invented the wheel, since it gives experience.

This subfolder does not require any dependence and can be smoothly compiled
by the elementary script "compile.sh"
Yes, I could have used LAPACK. But not now; the priority is somewhere else.

Then there are all the other folders, whose purpose is, in brief,
to implement new specific G operators and test the bayesian inversion on them.

****very important****
To implement "an operator G" is not a marginal task.
It might involve PDEs and Deep Neural Networks, for citing few examples.
The act of programming them is by itself an active learning step that
extend considerabily the usage of my library.

Here the name: Octopus, as sign of extension in multiple directions.

mylib is the head, every operator implementation a new tentacle.
