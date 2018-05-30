#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

#include <gsl/gsl_randist.h>

void GetH (run_params p, int model, int dim, double& h, double d, double pe_best, double ps_best, double s0_best, double p1_best, double p2_best, double de_best);
void FindLikelihoodRK (run_params p, double& h, double theta, double& lL, int maxv, int verb, int verb_model, int model, int dim, double d, double s0, double p1, double p2, double de, vector<double> cd81_val, vector<double> srb_val, vector< vector<dat> >& all_dat, vector<double> fact_store);
void MakeGrid (double& h, int verb_model, int model, int dim, double d, double e0, double s0, double p1, double p2, double de);
void ConstructPVals (int maxv, vector<double> evals_cd81, vector<double> evals_srb1, vector< vector<double> >& pvals_cd81, vector< vector<double> >& pvals_srb1);
double CalcLike(int verb, double theta, vector< vector<double> >& pvals_cd81, vector< vector<double> >& pvals_srb1, vector< vector<dat> >& all_dat, vector<double>& fact_store);
double BasicModelRK0 (int dim, int verb, double h, double d, double pe, double e0, double p1, double de); //Not actually implemented
double BasicModelRK1 (int dim, int verb, double h, double d, double pe, double ps, double s0, double p1, double p2);
double BasicModelRK2 (int dim, int verb, double h, double d, double pe, double ps, double s0, double p1, double p2, double de);


