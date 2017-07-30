//Setup parameters for fitness function,
//this should be as general as possible i.e doesn't need to change for different fitness funcitons
struct problem setupFitnessFunction(size_t nDim, void *ffParams){
	int d;
	struct problem pb;

	struct fitFuncParams *inParams = (struct fitFuncParams *)ffParams; //use inParams to extract information from ffParams

	pb.epsilon = 0.00000;	// Acceptable error (default). May be modified below
	pb.objective = 0;       // Objective value (default). May be modified below
	pb.SS.quantisation=0;		// No quantisation needed (all variables are continuous)
	pb.SS.normalise=0; // Set to a value x.
										//  x>0 => Normalisation is applied (search space => [0,x]^D)

	pb.SS.D = nDim;

	// Set boundaries equal to standarized coords
	for (d = 0; d < pb.SS.D; d++)	{
		pb.SS.min[d] = 0;   // Real coord = gsl_vector_get(inParams->rmin, d);
		pb.SS.max[d] = 1;   //Real coord = -1*gsl_vector_get(inParams->rmin, d);
		pb.SS.q.q[d] = 0;
	}

	pb.evalMax =75000; //3200;
	pb.epsilon=50; //0.001;
	pb.objective=0;


	if(pb.SS.normalise>0) // Normalise the quanta
	{
			for (d = 0; d < pb.SS.D; d++)
				pb.SS.q.q[d]=pb.SS.q.q[d]*pb.SS.normalise/(pb.SS.max[d]-pb.SS.min[d]);
	}

	pb.SS.q.size = pb.SS.D;

	return pb;
}
