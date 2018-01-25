
Overview
--------

This folder contain the source code for ELMOEA/D (Extreme learning surrogate models in multi-objective optimization based on decomposition). The src folder contain the ELMOEA/D source files and the lib folder its dependencies.

Building Instructions
---------------------

The project is build with CMake version 2.8>. Other dependencies are a C++ compiler compatible with C++11 standard and gnuplot if plotting callback is used.

To build the executables execute the following in elmoead_benchmarks folder (assuming a Unix like system):

$ mkdir build && cd build
$ cmake ..
$ make

This will build the executables:

moead_surrogate - ELMOEA/D algorithm
moead_alg - Non-surrogate MOEA/D algorithm
test_all - unit tests

Running ELMOEA/D
----------------

To use ELMOEA/D run following command:

./moead_surrogate config_file.cnf final_pop.dat

where config_file.cnf is the parameters files and final_pop.dat the final population objectives data.

The configuration file have the following syntax:
<parameter name> = <value>

The supported parameters are:

seed <long> # random number generator seed

prob_name <string> # problem name (WFGn, ZDTn, UFn, DTLZn etc)
prob_no_vars <long> # number of variables
prob_no_objs <long> # number of objectives
hyp_ref <real real ...> # hypervolume reference point (prob_no_objs real variables)

de_operator <string> # differential evolution algorithm (DE/RAND/1/BIN, DE/RAND/2/BIN, DE/NONLINEAR)
de_cr <real> # differential evolution algorithm crossover rate
de_f <real> # differential evolution algorithm scaling parameter
weight_type <string> # moea/d weight generation (UNIFORM, RANDOM)
moead_type <string> # moea/d type (MOEA/D, MOEA/D-DE, MOEA/D-STM)
moead_aggr_func <string> # moead aggregation function (TCHEBYCHEFF, INVERTED_TCHEBYCHEFF PENALTY_BOUNDARY_INTERSECTION etc)
moead_no_partitions <long> # number of partitions in moead weight set
moead_no_neighbors <long> # number of individuals per neighborhood
moead_no_updates <long> # number of moead updates per generation
moead_real_b <real> # moea/d-de chance to choose from neighborhood
moead_sbx_xover <read> # moead sbx crossover rate
moead_sbx_eta <read> # moead sbx eta parameter
moead_pm_mut <read> # moead mutation rate
moead_pm_eta <read> # moead polynomial mutaion eta parameter

elm_act_func <string> # elm activation function (SIGMOID, RBF, GAUSSIAN, MULTIQUADRIC etc)
elm_C <real> # elm regularization parameter", 1.0);
elm_no_hidden <long> # elm number of hidden neurons
elm_norm_input <boolean> # normalize input variables flag
elm_norm_output <boolean> # normalize output variables flag

no_evals # maximum number of real evaluations
no_sel_partitions # elmoead number of partitions for selection weights
archive_eps # elmoead archive minimum difference

callback_type # callback function called every iteration (NONE, PLOT, LOG)

If the parameter is not provided in the configuration file, a default value will be set and printed on the execution.

Disclamer
---------

I do not own any of the libraries contained in the lib/ folder. Please respect their respective licences.

If you use this code please cite the paper:

Lucas M. Pavelski, Myriam R. Delgado, Carolina P. Almeida, Richard A. Gonçalves, Sandra M. Venske, Extreme Learning Surrogate Models in Multi-objective Optimization based on Decomposition, Neurocomputing, Available online 6 November 2015, ISSN 0925-2312, http://dx.doi.org/10.1016/j.neucom.2015.09.111.
(http://www.sciencedirect.com/science/article/pii/S0925231215016094)

This work was also based on the papers (please also cite them if relevant):

Saúl Zapotecas Martínez and Carlos A. Coello Coello. 2013. MOEA/D assisted by rbf networks for expensive multi-objective optimization problems. In Proceedings of the 15th annual conference on Genetic and evolutionary computation (GECCO '13), Christian Blum (Ed.). ACM, New York, NY, USA, 1405-1412. DOI=http://dx.doi.org/10.1145/2463372.2465805.

H. Li and Q. Zhang,  Multiobjective Optimization Problems with Complicated Pareto Sets,  MOEA/D and NSGA-II, IEEE Trans on Evolutionary Computation, vol. 12,  no 2,  pp 284-302, April/2009.

G.-B. Huang, H. Zhou, X. Ding, and R. Zhang, “Extreme Learning Machine for Regression and Multiclass Classification,” IEEE Transactions on Systems, Man, and Cybernetics - Part B: Cybernetics,  vol. 42, no. 2, pp. 513-529, 2012.

