#include "stdio.h"
#include "math.h"
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

void wrapper_sPSO(size_t, /* Dimensionality of fitness function */
            double (*)(gsl_vector *, void *), /* Pointer to fitness function */
		    void *, /* Fitness function parameter structure */
			struct psoParamStruct *, /* PSO parameters */
			struct returnData * /* Structure containing PSO output */
		   );
