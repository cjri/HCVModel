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
	
	//Import data
	vector< vector<dat> > all_dat;
	GetAllDat(all_dat);
	//Reduce total number of cells to total at 6h using growth factor
	//Calculate MOI given total number of viruses
	ProcessData(all_dat);
	
	CorrectOEData(all_dat);  //Linear correction of over-expression data
	int maxv=0;
	GetMaxV(maxv,all_dat);
	//cout << maxv << "\n";
	double theta=0.359234; //Parameter characterising likelihood variance
	vector<double> fact_store;
	GetFactStore(fact_store,1000);
	
	//Question here - what simulations need to be run
	vector<double> cd81_val;
	vector<double> srb_val;
	GetReductions(all_dat,cd81_val,srb_val);
	IndexData(cd81_val,srb_val,all_dat);

	
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
	double s0=p.s0;  //Rate of encountering SRB1
    double p1=p.p1; //Probability of binding CD81 with no SRB1 x Probability of encountering CD81
	double p2=p.p2; //Probability of binding CD81 with SRB1 x Probability of encountering CD81
	double de=p.de; //Rate of downstream entry process
	
	if (p.lambda!=0) {
		p1=s0/p.lambda;
	}
	
	cout << "Dimension " << dim << "\n";
	
	//Calculate initial value of h for Runge-Kutta method
	double lL=-1e10;
	double h=100;
	if (p.grid<=2) {
		GetH(p,model,dim,h,d,1,1,s0,p1,p2,de);
	} else {
		GetH(p,model,dim,h,d,p.gp1,p.gp2,s0,p1,p2,de);
	}
	cout << h << "\n";

	if (p.grid==2) {
		for (int i=0;i<20;i++) {
			double e_inf=BasicModelRK2 (dim,0,h,d,i,i,s0,p1,p2,de);
			cout << i << " " << i << " " << e_inf << "\n";
		}
		return 0;
	}
	
	if (p.grid==3) {
		double e_inf=BasicModelRK2 (dim,0,h,d,p.gp1,p.gp2,s0,p1,p2,de);
		cout << p.gp1 << " " << p.gp2 << " " << e_inf << "\n";
		return 0;
	}

	
	if (p.grid==1) {
		double i=0;
		for (int ic=0;ic<=20;ic++) {
			double j=0;
			for (int jc=0;jc<=20;jc++) {
				double e_inf=BasicModelRK2 (dim,0,h,d,i,j,s0,p1,p2,de);
				cout << i << " " << j << " " << e_inf << "\n";
				j=j+0.1;
			}
			i=i+0.1;
		}
		return 0;
	} else {
		
		int first=1;
		double lL_best=-1e10;
		double s0_best=s0;
		double p1_best=p1;
		double p2_best=p2;
		double de_best=de;
	
		double n_acc=0;
		double n_try=0;
		double rate=0;
		double ds=0.01;
		if (p.smallchange==1) {
			ds=0.001;
		}
		int show_all=p.show_all;
		ofstream rec_file;
		rec_file.open("Record.out");
		//Find WT likelihood
		for (int it=0;it<=p.max_its;it++) {
			//cout << "Iteration " << it << "\n";
			if (first==0) {
				if (lL>lL_best) {
					PrintParameters (model,lL,s0,p1,p2,de,rec_file);
					lL_best=lL;
					s0_best=s0;
					p1_best=p1;
					p2_best=p2;
					de_best=de;
					if (p.precision>0) {
						FindLikelihoodRK(p,h,theta,lL,maxv,1,0,model,dim,d,s0_best,p1_best,p2_best,de_best,cd81_val,srb_val,all_dat,fact_store);
					}
					n_acc++;
				} else {
					s0=s0_best;
					p1=p1_best;
					p2=p2_best;
					de=de_best;
				}
		
				if (it>2&&it%1000==1) {
					rate=(n_acc+0.)/(n_try+0.);
					cout << "Old ds " << ds << " " << n_acc << " " << n_try << "\n";
					ds=ds*(rate+0.95);
					cout << "Rate " << rate << "\n";
					cout << "New ds " << ds << "\n";
					n_acc=0;
					n_try=0;
				}
			}
			//Change parameters
			n_try++;
			if (p.lambda==0) {
				ChangeParameters (model,s0,p1,p2,de,ds,rgen);
			} else {
				ChangeParametersFix (model,p.lambda,s0,p1,p2,de,ds,rgen);
			}
			first=0;
			lL=0;
			if (it==0) {
				FindLikelihoodRK(p,h,theta,lL,maxv,2,show_all,model,dim,d,s0_best,p1_best,p2_best,de_best,cd81_val,srb_val,all_dat,fact_store);
			} else {
				FindLikelihoodRK(p,h,theta,lL,maxv,0,0,model,dim,d,s0,p1,p2,de,cd81_val,srb_val,all_dat,fact_store);
				//cout << "lL " << lL << " " << lL_best << "\n";
			}
		
		}
	}
	
}
