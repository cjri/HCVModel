#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>


void GetData (string name, int cd, vector<dat>& dat1);
void GetAllDat (vector< vector<dat> >& all_dat);
void ConstructOutputVectors(vector< vector<dat> >& all_dat, vector<double>& cd81_indices, vector<double>& srb_indices, vector<double>& cd81_x, vector<double>& srb_x);
void ModifyOutputVectors (int i, int j, vector< vector<dat> >& all_dat, double p, vector<double>& cd81_indices, vector<double>& srb_indices, vector<double>& cd81_x, vector<double>& srb_x);
void OutputData (vector< vector<dat> >& all_dat, vector<double>& cd81_x, vector<double>& srb_x);
void OutputFreqs (vector<double>& cd81_indices, vector<double>& srb_indices, vector<double>& cd81_x, vector<double>& srb_x);
void PrintParameters (double lL, vector<double> params, ofstream& rec_file);
