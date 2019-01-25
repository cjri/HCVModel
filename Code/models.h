using namespace std;
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

#include <gsl/gsl_randist.h>

void GetH (run_params p, int dim, double& h, double pe_best, double ps_best, vector<double>& params);
void CheckH (run_params p, int dim, double& h, double& h_out, double pe_best, double ps_best, vector<double>& params);
double BasicModelRK0 (int dim, int verb, int quickprop, double h, double pe, double ps, double d, double c1, int Ne, double de);
double BasicModelRK6 (int dim, int verb, int quickprop, double h, double pe, double ps, double d, int Ne, double c1, double c2, double s1, double de);
double BasicModelRK21 (int dim, int verb, int quickprop, double h, double pe, double ps, double d, int Ne, double c1, double c2, double c3, double s1, double de);
double BasicModelRK22 (int dim, int verb, int quickprop, double h, double pe, double ps, double d, int Ne, double c1, double c2, double s1, double s2, double de);


void MatrixSetup (int dim, vector<double>& M, vector<double>& dM, vector<double>& zero);
void MatrixSetup6 (int dim, int Ne, vector< vector<double> >& M, vector< vector<double> >& dM, vector< vector<double> >& zero, vector< vector<int> >& mask);
void HalfTemp(vector<double>& temp);
void HalfTempMat(vector< vector<double> >& temp);
void MakeMask (double cut, vector< vector<double> >& M, vector< vector<int> >& mask);
void CleanupMatrix (int& change, double cut, double N, double pec1, double psc2, double r, vector<int>& cs, vector< vector<double> >& M);
void FinalCalc (int& rs, double& d, double& de, double& D, double& E, vector<double>& M);
void FinalCalcMat (int& rs, vector<int>& cs, double& d, double& de, double& D, double& E, vector< vector<double> >& M);
void UpdateColumns (int& rs, double cut, vector<double>& M);
void CleanupMatrix (int& change, double cut, double Ne, double pec1, double pec2, double pss1, vector<int>& cs, vector< vector<double> >& M);
void MakeMask (double cut, vector< vector<double> >& M, vector< vector<int> >& mask);
void UpdateColumnsMat (int& change, int& rs, vector<int>& cs, vector< vector<double> >& M);




void FindLikelihoodRK (run_params p, double cons, double t, double k, double& h, double theta, double& lL, int maxv, int verb, int verb_model, int dim, vector<double> params, vector<double> cd81_val, vector<double> srb_val, vector< vector<dat> >& all_dat, vector<double> fact_store);


void ConstructPVals (int maxv, vector<double> evals_cd81, vector<double> evals_srb1, vector< vector<double> >& pvals_cd81, vector< vector<double> >& pvals_srb1);
double CalcLike(int verb, double theta, vector< vector<double> >& pvals_cd81, vector< vector<double> >& pvals_srb1, vector< vector<dat> >& all_dat, vector<double>& fact_store);
double BasicModelRK0 (int dim, int verb, double h, double d, double pe, double e0, double p1, double de); //Not actually implemented
double BasicModelRK1 (int dim, int verb, double h, double d, double pe, double ps, double s0, double p1, double p2);
double BasicModelRK2 (int dim, int verb, double h, double d, double pe, double ps, double s0, double p1, double p2, double de);


