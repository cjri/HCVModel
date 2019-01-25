#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>


void ProcessData (vector< vector<dat> >& all_dat);
void CorrectOEData (vector< vector<dat> >& all_dat);
void GetMaxV (int& maxv, vector< vector<dat> >& all_dat);
void GetReductions (vector< vector<dat> > all_dat, vector<double>& cd81_val, vector<double>& srb_val);
void IndexData (vector<double> cd81_val, vector<double> srb_val, vector< vector<dat> >& all_dat);
void GetFactStore(vector<double>& fact_store,int N);
void GetParams (run_params p, vector<double>& params);
void CollectParams (run_params p, vector<double>& params, double& d, double& c1, double& c2, double& c3, double& s1, double& s2, double& Nc, double& Ns, int& Ne, double& de);

void ChangeParameters (run_params p, vector<double>& params, vector<double>& params_fix, vector<double> ds, gsl_rng *rgen);

double CalcDoublePoisson(double mu, double theta, int y, vector<double> fact_store);
double CalcGamma (double& cons, double& t, double& k, double x);


