#if !defined(SPSOWRAPPERHDR)
#define SPSOWRAPPERHDR

#include "sPSO.h"
#include "maxphase.h"
#include "ptapso.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

//Calls fitness function from sPSO algorithm
double evaluateFitness (size_t nDim, struct position x, ptrfunction function, void *ffParams, struct SS SS, double objective);

//sets up fitness function parameters such as xMin and xMax for sPSO
struct problem setupFitnessFunction(size_t nDim, void *ffParams);

//Fitness functions
double griewankFunction(struct position x);
double rastriginFunction(struct position x);

//sPSO algorithm
struct result sPSO (size_t nDim ,struct param param, struct problem problem, ptrfunction function, void *ffParams, struct psoParamStruct *psoParams);

void wrapper_sPSO(size_t, /* Dimensionality of fitness function */
            double (*)(gsl_vector *, void *), /* Pointer to fitness function */
		    void *, /* Fitness function parameter structure */
			struct psoParamStruct *, /* PSO parameters */
			struct returnData * /* Structure containing PSO output */
		   );

struct problem setupFitnessFunction(size_t nDim, void *ffParams);

#endif