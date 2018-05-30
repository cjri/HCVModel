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

void ChangeParameters (int model, double& s0, double& p1, double& p2, double& de, double ds, gsl_rng *rgen);
void ChangeParametersFix (int model, double lambda, double& s0, double& p1, double& p2, double& de, double ds, gsl_rng *rgen);

double CalcDoublePoisson(double mu, double theta, int y, vector<double> fact_store);


