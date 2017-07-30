#include "wrapper.h"
#include "maxphase.h"
#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>





void wrapper_sPSO(size_t nDim, /*!< Number of search dimensions */
            double (*fitfunc)(gsl_vector *, void *), /*!< Pointer to Fitness function */
			void *ffParams, /*!< Fitness function parameter structure */
            struct psoParamStruct *psoParams, /*!< PSO parameter structure */
			struct returnData *psoResults /*!< Output structure */){


    struct position bestBest; // Best position over all runs
  	int d;			// Current dimension
  	double D;
  	double error;			// Current error
  	double errorMean;		// Average error
  	double errorMin;		// Best result over all runs
  	double errorMeanBest[R_max];
  	double evalMean;		// Mean number of evaluations

  	int nFailure;		// Number of unsuccessful runs
  	double logProgressMean;
  	struct param param;
  	struct problem pb;
  	int randCase;
  	int run = 0;
    int runMax = 1;//, runMax;
  	struct result result;

  	int scanNb;
  	double Smean;
  	double success[funcMax];
  	double successRate;
  	int t;
  	double variance;
  	float z;
  	double zz;

  	E = exp ((long double) 1);
  	pi = acos ((long double) -1);
  	errMax=0;
  	nbRand=0;


  	//------------------------------------------------ PARAMETERS
  			// Bells and Whistles
  				// Not really part of the standard
  				// May improve the performance (not always)
  				// Sometimes not well mathematically founded (rules of thumbs)
  				// * => suggested value

  	param.BW[0]=0; 	//	0 => same swarm size for each run
  									//	1 => random swarm size around the given mean

  	param.BW[1]=0; 	//*	0 => when P=G, use the "standard" method
  									// 	1 => when P=G, a specific probabilistic method
  									//	2 => when P=G, a more conservative method
  									//  3 => when P=G, just look around
  									// 4 =>  different weights for X, P, G  (TEST)

  	param.BW[2]=-1;	// Randomness options
  	                // -2nn => Truncated KISS (simulated).
  									//			nn is the number of bits you use to define
  									//			each random number/ Example: -207 for 7 bits
  	                // -1 => "native" rand() of the C language
  									//* 0 => pseudo-random number generator KISS
  									//* 10 => pseudo-random number generator Mersenne 64 bits
  									// 1 => quasi-random numbers for initialisation, Sobol sequences
  									//      KISS after that
  									// 2 => quasi-random numbers for initialisation, Halton sequences
  									//      KISS after that
  									// 3nn => Read on a list of bits (f_rand_bin.txt).
  									//      Normally coming from a "true" (physical) random number generator
  									//      (quantic system, atmospheric noise ...)
  									//			nn is the number of bits you use to define
  									//			each random number/ Example: 307 for 7 bits
  									//      Warning: nn must be >=2
  									// 4 => Read on a list (f_rand_quasi.txt)


  	param.BW[3]=0;	// 1 => random numbering of the particles before each iteration
  									// 0 => always the same loop "particle 0 to S-1"

  	//--------

  	param.confin=0; 	// 0 => keep inside the search space (supposed to be a D-rectangle)
  										// 1 => no confinement
  										//   WARNING: may be very slow (and bad) for discrete problems

  	param.distrib=0; // -1 => uniform in the hypersphere
  										//* 0 => in the hypersphere, uniform along the radius
  										// 			(and therefore NOT uniform in the sphere)
  										// 1 => Gaussian (Box-Muller method). Warning: infinite loop possible
  										// 2 =>	Gaussian (CMS method)
  										// 3 => Other stable (CMS, experimental parameters)
  										// 4 => Slash distribution (Gaussian BM/Gaussian BM)
  				// Useful only if param.distrib>0;
  				param.mean=0.5; //Default: 0.5. For some functions 0 is better, though
  												//	Example: shifted Rosenbrock (code 102)
  				param.sigma=1./12; // Default: 1./12 (standard deviation of U(0,1))
  		// WARNING: the CMS method may not work with randomness option >=2
  	if(param.BW[2]>=2) param.distrib=	0;

  	Smean=40; //Swarm size or Mean swarm size (if BW[0]=1). Suggested: 40

  	param.K=3; 	// Parameter to compute the probability p for a particle to be an
  							// external informant. You may also directly define p (see below),
  							// but K is about the mean number of the informants of a particle.
  							// Default: 3

  	// Confidence coefficients. Default:
  	param.w = 1. / (2 * log ((double) 2)); // 0.721
  	param.c = 0.5 + log ((double) 2); // 1.193
  	param.topology = 0; // 0 => information links as in SPSO 2007 (quasi-random)
  											// 1 => variable random ring (EXPERIMENTAL)

    //-------------------------------------------------- False randomnesses
    switch (param.BW[2]) //
    {
      default: break;
      case 4: // Prepare a list of false random number, read on a file
      t=0;
      readRand:
      scanNb=fscanf(f_rand,"%f",&z);
      if(scanNb!=EOF)
      {
        randNumber[t]=z;
        t=t+1;
      goto readRand;
      }
      nCycleMax=t;
      printf("\n%i false random numbers read on a file",nCycleMax);

      break;
    }
  	// ----------------------------------------------- PROBLEM
  	param.trace=0; // If >0 more information is displayed/saved (f_trace.txt)
  	// Functions to optimise


  		// Define the problem
  		pb=setupFitnessFunction(nDim,ffParams);
  		//pb = problemDef(3,functionCode);
  		if(pb.SS.D>DMax) ERROR ("Can't solve it. You should increase DMax");

  		// ----------------------------------------------- RUNS
  		errorMean = 0;
  		evalMean = 0;
  		nFailure = 0;
  		D=pb.SS.D;

  		randCase=param.BW[2];
  		if(randCase>300) {nBit=randCase-300; randCase=3;} // Decode the number of bits
  		if(randCase<-200) {nBit=-randCase-200; randCase=-2;}

  		switch(randCase)
  		{
  		default:
  		 break;

  		case 0:
  		seed_rand_kiss(1294404794); // Initialise the RNG KISS for reproducible results
  		break;

  		case 10: // Mersenne 64 bits
  		init_genrand64(1294404794);
  		//init_genrand64(1234567890);
  		break;

  		case -2: // Truncated KISS (simulated)
  		rMax=pow(2,nBit)-1;
  		break;

  		case 3: // The file is a string of bits
      //nBit=3;
  		rMax=pow(2,nBit)-1; // For conversion of a nBit string into a number in [0,1]
  		break;

  		case 4: // The file directly contains the numbers
      nCycle=0;
      break;
  		}


  randRank=0; randChaos=0.02;

  	/*
  		seconds=time(NULL); // Initialise the RNG KISS more randomly
  		 printf("\n time %ld",seconds);
  		seed_rand_kiss(time(NULL));
  */
  	switch(randCase) // "Warm up" the RNG for pseudo-random numbers
  	{
  		default:
  		for (t=0;t<10000;t++) zz=alea(0,1,randCase);
  		break;

  		case 3:
  		case 4:
  		break;
  	}
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /*                      PSO  LOOP                              */
    //////////////////////////////////////////////////////////////////////////////////////////////////////
  	//	for (run = 0; run < runMax; run++)
  	//	{

  			if(param.BW[0]==0) param.S= psoParams->popsize; // Constant swarm size
  			else // Random swarm size "around" the mean value
  			param.S=(int)(0.5*(0.5+alea(Smean-D/2,Smean+D/2,0)+alea(Smean-D/2,Smean+D/2,0)));

  			param.p=1-pow(1-1./param.S,param.K); //for a "global best" PSO, directly set param.p=1

  			result = sPSO (nDim, param, pb, fitfunc, ffParams, psoParams);
  			error = result.error;

  			if (error > pb.epsilon) // Failure
  				nFailure = nFailure + 1;

  			if(pb.SS.normalise>0)
  			{
  				for(d=0;d<pb.SS.D;d++)
  					result.SW.P[result.SW.best].x[d]=
  					pb.SS.min[d]+(pb.SS.max[d]-pb.SS.min[d])*result.SW.P[result.SW.best].x[d]
  					/pb.SS.normalise;
  			}

  			// Memorize the best (useful if more than one run)
  			if(run==0) {
          bestBest=result.SW.P[result.SW.best];
          for (d=0;d<pb.SS.D;d++) {
              gsl_vector_set(psoResults->bestLocation, d, bestBest.x[d]);

          }
          psoResults->bestFitVal = bestBest.f;

        }
  			else {
  				if(error<bestBest.f) {
            bestBest=result.SW.P[result.SW.best];
            for (d=0;d<pb.SS.D;d++) {
                gsl_vector_set(psoResults->bestLocation, d, bestBest.x[d]);

            }
            psoResults->bestFitVal = bestBest.f;

          }
        }

  			// Result display
  			errorMean=errorMean+error;
  			printf ("\nRun %i. S %i,  Eval %f. Error %e ", run+1, param.S, result.nEval, error);
  			printf(" Mean %e",errorMean/(run+1));
  			zz=100*(1-(double)nFailure/(run+1));
  			printf("  Success  %.2f%%",zz);

  			// Best position display
  			//for (d=0;d<pb.SS.D;d++) printf(" %f",result.SW.P[result.SW.best].x[d]);

  			// Save result

  			// Save best position
  			//		for ( d = 0; d < pb.SS.D; d++ ) fprintf( f_run, " %f",  result.SW.P[result.SW.best].x[d] );

  			// Compute/save some statistical information
  			if (run == 0)
  				errorMin = error;
  			else if (error < errorMin)
  				errorMin = error;

  			evalMean = evalMean + result.nEval;
  			errorMeanBest[run] = error;
  			logProgressMean  = logProgressMean - log(error);
  	//	}		// End loop on "run"

      //Convert resutl outputs to psoResult struct
      psoResults->totalFuncEvals = result.nEval * runMax;
      psoResults->totalIterations = result.nIterations * runMax;



  		// ---------------------END

  		// Display some statistical information
  		evalMean = evalMean / (double) runMax;
  		errorMean = errorMean / (double) runMax;
  		logProgressMean = logProgressMean/(double) runMax;

/*
  		printf ("\n Eval. (mean)= %f", evalMean);
  		printf ("\n Error (mean) = %e", errorMean);

  		// Variance
  		variance = 0;

  		for (run = 0; run < runMax; run++)
  			variance = variance + pow (errorMeanBest[run] - errorMean, 2);

  		variance = sqrt (variance / runMax);
  		printf ("\n Std. dev. %e", variance);
  		printf("\n Log_progress (mean) = %f", logProgressMean);

  		// Success rate and minimum value
  		printf("\n Failure(s) %i  ",nFailure);

  		successRate = 100 * (1 - nFailure / (double) runMax);
  		printf ("Success rate = %.2f%%", successRate);
  		success[0]=successRate;

  		printf ("\n Best min value = %1.20e", errorMin);
  		printf ("\nPosition of the optimum: ");
  		for (d=0;d<pb.SS.D;d++) printf(" %.20f",bestBest.x[d]);




      printf("\n errMax : %f",errMax);

  		// Repeat informations
  	  printf("\n---------");
  	  printf("\n Function(s):");

  		printf("\n Confinement: ");
  		if(param.confin==0) printf("YES"); else printf("NO");
  		printf("\n Distribution: ");
  		switch(param.distrib)
  		{
  			case 0: printf(" uniform"); break;
  			case 1: printf(" Gaussian (%f,%f), Box-Muller",param.mean,param.sigma); break;
  			case 2: printf(" Gaussian (%f,%f), CMS",param.mean,param.sigma); break;
  			case 3: printf(" Stable (%f,%f)",param.mean,param.sigma); break;
  			case 4: printf(" Slash (%f,%f)",param.mean,param.sigma); break;
  		}

  		printf("\n BW = (%i, %i, %i, %i)",param.BW[0],param.BW[1],
  		       param.BW[2],param.BW[3]);
  		printf("\n Swarm size: ");
  		if(param.BW[0]==0) printf("%i",(int)psoParams->popsize); else printf(" mean %i",(int)Smean);
  		printf("\n K = %i",param.K);
  		printf("\n w = %f",param.w);
  		printf("\n c = %f",param.c);
  		printf("\n %e random numbers have been used",nbRand);
      */

}
