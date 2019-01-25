using namespace std;
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

struct run_params {
	string ref;
	int dim;	//Model dimension
	int model;	//Model index
	int rk;		//Flag for Runge-Kutta method
	int over;	//Flag for adding over-expression data
	int max_its;
	double e0;
	double s0;
	double p1;
	double p2;
	double de;
	double r;
	double N;
	double c1; //SRB1-independent binding of CD81
	double c2; //SRB1-dependent binding of CD81
    double c3; //Unbinding of CD81
	double s1; //Binding of SRB1
    double s2; //Unbinding of SRB1
    double Nc; //Number of CD81 proteins available
    double Ns; //Number of SR-B1 proteins available
    int Ne; //Number of E2 proteins involved in binding
	double lambda;
    double usegamma;
    int precision; //Makes the Runge-Kutta algorithm more accurate, and more consistent.  Use for final calculations
    int smallchange; //Start with small ds in optimisation
    double threshold; //Termination criterion for ds
	int grid;
	double gp1; //Used to specify particular grid points
	double gp2; //Used to specify particular grid points
    int nfix; //Flag to fix the parameter Ne in the optimisation
};

struct dat {
	double a; //Proportion of available receptor
	double sd; //Standard deviation in a
	double growth; //Growth rate of cells during experiment
	double count; //Number of cells at 48h
	int foci; //Number of cells infected
	int cd81; //CD81=1, SRB1=0
	int index; //Which proportion is this?
	double moi; //MOI from number of viruses
	vector<double> pvir; //Probability of k viruses attached to a cell
};


void GetH (int dim, double& h, double d, double m_best, double pe_best, double ps_best, double e0_best, double s0_best, double p1_wt, double p2_wt);

double BasicModelEuler (int dim, double d, double m, double pe, double ps, double e0, double s0, double p1, double p2);
void GetProps (vector<double>& pe_vec, vector<double>& ps_vec);
void GetOptions (run_params& p, int argc, const char **argv);
//void GetData (vector<double>& dat1_mean, vector<double>& dat1_sd, vector<double>& dat2_mean, vector<double>& dat2_sd);
//void FindLikelihood (int model, double& lL, int verb, int dim, double d, double m, double e0, double s0, double p1, double p2, double de, vector<double> dat1_mean, vector<double> dat1_sd, vector<double> dat2_mean, vector<double> dat2_sd, vector<double> pe_vec, vector<double> ps_vec);
//void FindLikelihoodRK (int model, double& h, double& lL, int verb, int verb_model, int dim, double d, double m, double e0, double s0, double p1, double p2, double de, vector<double> dat1_mean, vector<double> dat1_sd, vector<double> dat2_mean, vector<double> dat2_sd, vector<double> pe_vec, vector<double> ps_vec);


double BasicModelRK1 (int dim, int verb, double h, double d, double m, double pe, double ps, double s0, double p1, double p2);

