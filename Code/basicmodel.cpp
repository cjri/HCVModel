#include <iostream>
#include <vector>
#include <string>

//Model 1: Standard model
//Model 2: With SR-B1 unbinding
//Model 3: With limited CD81
//Model 4: With limited SR-B1
//Model 5: With CD81 unbinding

using namespace std;

#include "basicmodel.h"
#include "models.h"
#include "io.h"
#include "utilities.h"

int main(int argc, const char **argv) {

	run_params p;
	GetOptions(p,argc,argv);

	//Gamma distributional fit
	double cons=0.000414216785953741;
	double k=2.903557980326837;
	double t=11.868631620231913;
	
	cout << "P.dim is " << p.dim << "\n";
	
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
	
    //Set up data
    vector<double> cd81_val;
	vector<double> srb_val;
	GetReductions(all_dat,cd81_val,srb_val);
	IndexData(cd81_val,srb_val,all_dat);
	
	//Dimension
	int dim=p.dim;
    cout << "Dimension " << dim << "\n";
	//Use Runge-Kutta model - initialise random seed
	int seed=(int) time(NULL);
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, seed);

    //Rates - have new parameters here
    vector<double> params;
    vector<double> params_fix;
    GetParams(p,params);
    
	//Calculate initial value of h for Runge-Kutta method
	double lL=-1e10;
	double h=100;
	if (p.grid<=2) {
        GetH (p,dim,h,1,1,params);
	} else {
		GetH (p,dim,h,p.gp1,p.gp2,params);
	}
	cout << h << "\n";
    double h2=h;
	
	
	if (p.grid==2) {
        double d=0.1;
        double c1=0.1;
        double c2=0.1;
        double c3=0.1;
        double s1=0.1;
        double s2=0.1;
        double Nc=100;
        double Ns=100;
        int Ne=100;
        double de=0.1;
		for (int i=0;i<20;i++) {
            double e_inf=-10;
            CollectParams(p,params,d,c1,c2,c3,s1,s2,Nc,Ns,Ne,de);
			if (p.model==0) {
				e_inf=BasicModelRK0(dim,0,1,h,i,i,d,c1,Ne,de);
            } else if (p.model==6) {
                e_inf=BasicModelRK6(dim,0,1,h,i,i,d,Ne,c1,c2,s1,de);
			} else if (p.model==21) {
				e_inf=BasicModelRK21(dim,0,1,h,i,i,d,Ne,c1,c2,c3,s1,de);
			} else if (p.model==22) {
				e_inf=BasicModelRK22(dim,0,1,h,i,i,d,Ne,c1,c2,s1,s2,de);
			}
			cout << i << " " << i << " " << e_inf << "\n";
		}
		return 0;
	}
	
	if (p.grid==3) {
        double d=0.1;
        double c1=0.1;
        double c2=0.1;
        double c3=0.1;
        double s1=0.1;
        double s2=0.1;
        double Nc=100;
        double Ns=100;
        int Ne=100;
        double de=0.1;
        double e_inf=-10;
        CollectParams(p,params,d,c1,c2,c3,s1,s2,Nc,Ns,Ne,de);
		if (p.model==0) {
			e_inf=BasicModelRK0(dim,0,1,h,p.gp1,p.gp2,d,c1,Ne,de);
        } else if (p.model==6) {
            e_inf=BasicModelRK6(dim,0,1,h,p.gp1,p.gp2,d,Ne,c1,c2,s1,de);
		} else if (p.model==21) {
			e_inf=BasicModelRK21(dim,0,1,h,p.gp1,p.gp2,d,Ne,c1,c2,c3,s1,de);
		} else if (p.model==22) {
			e_inf=BasicModelRK22(dim,0,1,h,p.gp1,p.gp2,d,Ne,c1,c2,s1,s2,de);
		}
        cout << p.gp1 << " " << p.gp2 << " " << e_inf << "\n";
		return 0;
	}
	
	if (p.grid==1) {
        double d=0.1;
        double c1=0.1;
        double c2=0.1;
        double c3=0.1;
        double s1=0.1;
        double s2=0.1;
        double Nc=100;
        double Ns=100;
        int Ne=100;
        double de=0.1;
        CollectParams(p,params,d,c1,c2,c3,s1,s2,Nc,Ns,Ne,de);
		double i=0;
		h=h/8;
		for (int ic=0;ic<=20;ic++) {
			double j=0;
			for (int jc=0;jc<=20;jc++) {
                double e_inf=-10;
				if (p.model==0) {
					e_inf=BasicModelRK0(dim,0,1,h,i,j,d,c1,Ne,de);
                } else if (p.model==6) {
                    e_inf=BasicModelRK6(dim,0,1,h,i,j,d,Ne,c1,c2,s1,de);
				} else if (p.model==21) {
					e_inf=BasicModelRK21(dim,0,1,h,i,j,d,Ne,c1,c2,c3,s1,de);
				} else if (p.model==22) {
					e_inf=BasicModelRK22(dim,0,1,h,i,j,d,Ne,c1,c2,s1,s2,de);
				}
				cout << i << " " << j << " " << e_inf << "\n";
				j=j+0.1;
			}
			i=i+0.1;
		}
		return 0;
	} else {
		
		int first=1;
        vector<double> params_best=params;
		double lL_best=-1e10;
	
		double n_acc=0;
		double n_try=0;
		double rate=0;
        vector<double> ds;
        vector<double> ds_fix;
        for (int i=0;i<params.size();i++) {
            ds.push_back(0.01);
            params_fix.push_back(0);
            ds_fix.push_back(0);
        }
        if (p.smallchange==1) {
            for (int i=0;i<params.size();i++) {
                ds[i]=0.001;
            }
        }
        if (p.model==0||p.model==6) {
            cout << "Here\n";
            ds[4]=5;
            ds_fix[4]=1;
            if (p.model==6&&p.nfix==1) {
                params_fix[4]=1;
            }
        }
		if (p.model==21||p.model==22) {
			ds[5]=5;
			ds_fix[5]=1;
			if (p.nfix==1) {
				params_fix[5]=1;
			}
		}
		
        //Change ds as required for other models
		ofstream rec_file;
		rec_file.open("Record.out");
		//Find WT likelihood
		for (int it=0;it<=p.max_its;it++) {
			//cout << "Iteration " << it << "\n";
			if (first==0) {
				if (lL>lL_best) {
                    double h_store=h;
                    CheckH (p,dim,h,h2,1,1,params);
                    if (h2!=h_store) {
                        cout << "Change h to " << h2 << "\n";
                        h=h2;
                        FindLikelihoodRK(p,cons,t,k,h,theta,lL,maxv,0,0,dim,params,cd81_val,srb_val,all_dat,fact_store);
                    }

					PrintParameters (lL,params,rec_file);
					lL_best=lL;
                    params_best=params;
					if (p.precision>0) {
						FindLikelihoodRK(p,cons,t,k,h,theta,lL,maxv,1,0,dim,params_best,cd81_val,srb_val,all_dat,fact_store);
					}
					n_acc++;
				} else {
                    params=params_best;
				}
		
				if (it>2&&it%1000==1) {
					rate=(n_acc+0.)/(n_try+0.);
					cout << "Old ds " << ds[0] << " " << n_acc << " " << n_try << "\n";
                    for (int i=0;i<ds.size();i++) {
                        if (ds_fix[i]!=1) {
                            ds[i]=ds[i]*(rate+0.975);
                            if (ds[i]>0.02) {
                                ds[i]=0.02;
                            }
                        }
                    }
					cout << "Rate " << rate << "\n";
					cout << "New ds " << ds[0] << "\n";
					n_acc=0;
					n_try=0;
                    if (ds[0]<p.threshold) {
                        break;
                    }
				}
                n_try++;
                ChangeParameters (p,params,params_fix,ds,rgen);
			}
			//Change parameters
            first=0;
			lL=0;
			if (it==0) {
				FindLikelihoodRK(p,cons,t,k,h,theta,lL,maxv,1,0,dim,params_best,cd81_val,srb_val,all_dat,fact_store);
			} else {
				FindLikelihoodRK(p,cons,t,k,h,theta,lL,maxv,0,0,dim,params,cd81_val,srb_val,all_dat,fact_store);
			}
		}
	}
}
