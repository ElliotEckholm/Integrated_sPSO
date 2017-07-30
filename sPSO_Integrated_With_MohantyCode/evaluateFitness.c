// Evaluates the fitness value for the particle of rank s

double evaluateFitness(size_t nDim, struct position x, ptrfunction function, void *ffParams, struct SS SS, double objective)
{

		int d;
		double f = 0;

		//sPSO uses struct position to store particles position
		struct position xs;

		//will contain particles postion, has to be used because fitfunc requires parameter gsl_vector *
		gsl_vector *v = gsl_vector_alloc (nDim);

		if(SS.normalise>0)
		{
			// Back to the real search space
			xs.size=x.size;
			//Normalize the postion and load gsl_vector with particle's postion
			for(d=0;d<xs.size;d++) {
				xs.x[d]=SS.min[d]+(SS.max[d]-SS.min[d])*x.x[d]/SS.normalise;
				gsl_vector_set(v, d, xs.x[d]);
			}
		}
		//Particle is outside search space
		else {
			xs=x;
			//load gsl_vector with particle's postion
			for(d=0;d<xs.size;d++) {
				gsl_vector_set(v, d, xs.x[d]);
			}
		}

		//Call fitness function
		f = function(v, ffParams);

		gsl_vector_free (v);


		//calcualte fitness value with respect to sPSO objective, if there is one
		//Objective set to 0 by default
    f=fabs(f-objective);
    if(f<errMin) errMin=f; // For information
    if(f>errMax) {if(f<infinity) errMax=f; else errMax=infinity;} // For information
    return f;
	}
