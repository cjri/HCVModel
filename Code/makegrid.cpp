#include <iostream>
#include <vector>
#include <string>

using namespace std;

#include "basicmodel.h"
#include "models.h"
#include "io.h"
#include "utilities.h"

int main(int argc, const char **argv) {
	
	run_params p;
	GetOptions(p,argc,argv);
	
	//Dimension
	int dim=p.dim;
	int model=p.model;
	//Use Runge-Kutta model
	int seed=(int) time(NULL);
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, seed);

    //Rates
    double d=0.01;	//Death rate; Set arbitrarily here
    double e0=p.e0;  //Rate of encountering CD81
	double s0=p.s0;  //Rate of encountering SRB1
    double p1=p.p1; //Probability of binding CD81 with no SRB1
    double p2=p.p2; //Probability of binding CD81 with SRB1
	double de=p.de; //Rate of downstream entry process
	
	cout << "Dimension " << dim << "\n";
	
	//Calculate initial value of h for Runge-Kutta method
	double lL=-1e10;
	double h=100;
	GetH(p,model,dim,h,d,1,1,e0,s0,p1,p2,de);
	cout << h << "\n";

    MakeGrid (h,0,model,dim,d,e0,s0,p1,p2,de);

    

}
