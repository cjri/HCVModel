#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "basicmodel.h"

void ProcessData (vector< vector<dat> >& all_dat) {
	//Calculate numbers of cells at 6h
	for (int i=0;i<all_dat.size();i++) {
		for (int j=0;j<all_dat[i].size();j++) {
			all_dat[i][j].count=all_dat[i][j].count/all_dat[i][j].growth;
		}
	}
	
	//Total number of viruses
	double n_vir=11125;
	n_vir=n_vir*0.454; //Correct for proportion of cells observed
	
	//Calculate Poisson model of number of viruses per cell
	for (int i=0;i<all_dat.size();i++) {
		for (int j=0;j<all_dat[i].size();j++) {
			all_dat[i][j].moi=0.767247352;
			for (int k=0;k<100;k++) {
				double p=gsl_ran_poisson_pdf(k,all_dat[i][j].moi);
				if (p>1e-10) {
					all_dat[i][j].pvir.push_back(p);
				}
			}
		}
	}
}

void CorrectOEData (vector< vector<dat> >& all_dat) {
	double s11=0.0135929/0.0122426;
	double s12=0.0147563/0.0110657;
	double s21=0.0135929/0.0173837;
	double s22=0.0147563/0.0153883;
	
	for (int i=4;i<=5;i++) {
		for (int j=0;j<all_dat[i].size();j++) {
			if (all_dat[i][j].cd81==1) {
				all_dat[i][j].count=all_dat[i][j].count/s11;
			} else {
				all_dat[i][j].count=all_dat[i][j].count/s12;
			}
		}
	}
	for (int i=6;i<=7;i++) {
		for (int j=0;j<all_dat[i].size();j++) {
			if (all_dat[i][j].cd81==1) {
				all_dat[i][j].count=all_dat[i][j].count/s21;
			} else {
				all_dat[i][j].count=all_dat[i][j].count/s22;
			}
		}
	}
}


void GetMaxV (int& maxv, vector< vector<dat> >& all_dat) {
	//Get maximum length of vector from dataset
	for (int i=0;i<all_dat.size();i++) {
		for (int j=0;j<all_dat[i].size();j++) {
			if (all_dat[i][j].pvir.size()>maxv) {
				maxv=all_dat[i][j].pvir.size();
			}
		}
	}
}

void GetReductions (vector< vector<dat> > all_dat, vector<double>& cd81_val, vector<double>& srb_val) {
	for (int i=0;i<all_dat.size();i++) {
		for (int j=0;j<all_dat[i].size();j++) {
			if (all_dat[i][j].cd81==1) {
				int seen=0;
				for (int k=0;k<cd81_val.size();k++) {
					if (cd81_val[k]==all_dat[i][j].a) {
						seen=1;
						break;
					}
				}
				if (seen==0) {
					cd81_val.push_back(all_dat[i][j].a);
				}
			} else {
				int seen=0;
				for (int k=0;k<srb_val.size();k++) {
					if (srb_val[k]==all_dat[i][j].a) {
						seen=1;
						break;
					}
				}
				if (seen==0) {
					srb_val.push_back(all_dat[i][j].a);
				}
			}
		}
	}
	sort(srb_val.begin(),srb_val.end());
	sort(cd81_val.begin(),cd81_val.end());
}

void IndexData (vector<double> cd81_val, vector<double> srb_val, vector< vector<dat> >& all_dat) {
	for (int i=0;i<all_dat.size();i++) {
		for (int j=0;j<all_dat[i].size();j++) {
			if (all_dat[i][j].cd81==1) {
				for (int k=0;k<cd81_val.size();k++) {
					if (all_dat[i][j].a==cd81_val[k]) {
						all_dat[i][j].index=k;
					}
				}
			} else {
				for (int k=0;k<srb_val.size();k++) {
					if (all_dat[i][j].a==srb_val[k]) {
						all_dat[i][j].index=k;
					}
				}
			}
		}
	}
}

void GetFactStore(vector<double>& fact_store,int N){
	double logN=0;
	fact_store.push_back(0);
	for (int i=1;i<=N;i++) {
		logN=logN+log(i);
		fact_store.push_back(logN);
	}
}

void GetParams (run_params p, vector<double>& params) {
    double d=0.01;	//Death rate; Set arbitrarily here
    params.push_back(d);
	if (p.model==0) {
		params.push_back(p.c1);
		params.push_back(p.de);
		if (p.usegamma==1) {
			params.push_back(p.s1);
		}
		params.push_back(p.Ne);
	}
    if (p.model==1) {
        params.push_back(p.c1);
        params.push_back(p.c2);
        params.push_back(p.s1);
        params.push_back(p.de);
    }
    if (p.model==2) {
        params.push_back(p.c1);
        params.push_back(p.c2);
        params.push_back(p.s1);
        params.push_back(p.s2);
        params.push_back(p.de);
    }
    if (p.model==3) {
        params.push_back(p.c1);
        params.push_back(p.c2);
        params.push_back(p.s1);
        params.push_back(p.Nc);
        params.push_back(p.de);
    }
    if (p.model==4) {
        params.push_back(p.c1);
        params.push_back(p.c2);
        params.push_back(p.s1);
        params.push_back(p.Ns);
        params.push_back(p.de);
    }
    if (p.model==5) {
        params.push_back(p.c1);
        params.push_back(p.c2);
        params.push_back(p.c3);
        params.push_back(p.s1);
        params.push_back(p.de);
    }
    if (p.model==6) {
        params.push_back(p.c1);
        params.push_back(p.c2);
        params.push_back(p.s1);
        params.push_back(p.Ne);
        params.push_back(p.de);
    }
	if (p.model==60) {
		params.push_back(p.c1);
		params.push_back(p.c2);
		params.push_back(p.s1);
		params.push_back(p.de);
	}
	if (p.model==7) {
		params.push_back(p.c1);
		params.push_back(p.c2);
		params.push_back(p.c3);
		params.push_back(p.s1);
		params.push_back(p.s2);
		params.push_back(p.de);
	}
	if (p.model==8) {
		params.push_back(p.c1);
		params.push_back(p.c2);
		params.push_back(p.c3);
		params.push_back(p.s1);
		params.push_back(p.Nc);
		params.push_back(p.de);
	}
	if (p.model==9) {
		params.push_back(p.c1);
		params.push_back(p.c2);
		params.push_back(p.c3);
		params.push_back(p.s1);
		params.push_back(p.Ns);
		params.push_back(p.de);
	}
	if (p.model==10) {
		params.push_back(p.c1);
		params.push_back(p.c3);
		params.push_back(p.de);
		if (p.usegamma==1) {
			params.push_back(p.s1);
		}
	}
	if (p.model==11) {
		params.push_back(p.c1);
		params.push_back(p.c2);
		params.push_back(p.s1);
		params.push_back(p.s2);
		params.push_back(p.Ne);
		params.push_back(p.de);
	}
	if (p.model==12) {
		params.push_back(p.c1);
		params.push_back(p.c2);
		params.push_back(p.c3);
		params.push_back(p.s1);
		params.push_back(p.Ne);
		params.push_back(p.de);
	}
	if (p.model==13) {
		params.push_back(p.c1);
		params.push_back(p.c2);
		params.push_back(p.c3);
		params.push_back(p.s1);
		params.push_back(p.s2);
		params.push_back(p.Ne);
		params.push_back(p.de);
	}
	if (p.model==21) {
		params.push_back(p.c1);
		params.push_back(p.c2);
		params.push_back(p.c3);
		params.push_back(p.s1);
		params.push_back(p.Ne);
		params.push_back(p.de);
	}
	if (p.model==22) {
		params.push_back(p.c1);
		params.push_back(p.c2);
		params.push_back(p.s1);
		params.push_back(p.s2);
		params.push_back(p.Ne);
		params.push_back(p.de);
	}
	if (p.model==23) {
		params.push_back(p.c1);
		params.push_back(p.c2);
		params.push_back(p.c3);
		params.push_back(p.s1);
		params.push_back(p.s2);
		params.push_back(p.Ne);
		params.push_back(p.de);
	}

}

void CollectParams (run_params p, vector<double>& params, double& d, double& c1, double& c2, double& c3, double& s1, double& s2, double& Nc, double& Ns, int& Ne, double& de) {
	if (p.model==0) {
		d=params[0];
		c1=params[1];
		de=params[2];
		if (p.usegamma==1) {
			s1=params[3];
		}
		double Ntemp=params[4]+0.0001;
		Ne=(int)Ntemp;
	}
    if (p.model==1) {
        d=params[0];
        c1=params[1];
        c2=params[2];
        s1=params[3];
        de=params[4];
    }
    if (p.model==2) {
        d=params[0];
        c1=params[1];
        c2=params[2];
        s1=params[3];
        s2=params[4];
        de=params[5];
    }
    if (p.model==3) {
        d=params[0];
        c1=params[1];
        c2=params[2];
        s1=params[3];
        Nc=params[4];
        de=params[5];
    }
    if (p.model==4) {
        d=params[0];
        c1=params[1];
        c2=params[2];
        s1=params[3];
        Ns=params[4];
        de=params[5];
    }
    if (p.model==5) {
        d=params[0];
        c1=params[1];
        c2=params[2];
        c3=params[3];
        s1=params[4];
        de=params[5];
    }
    if (p.model==6) {
        d=params[0];
        c1=params[1];
        c2=params[2];
        s1=params[3];
        double Ntemp=params[4]+0.0001;
        Ne=(int)Ntemp;
        de=params[5];
    }
	if (p.model==60) {
		d=params[0];
		c1=params[1];
		c2=params[2];
		s1=params[3];
		de=params[4];
	}
	if (p.model==7) {
		d=params[0];
		c1=params[1];
		c2=params[2];
		c3=params[3];
		s1=params[4];
		s2=params[5];
		de=params[6];
	}
	if (p.model==8) {
		d=params[0];
		c1=params[1];
		c2=params[2];
		c3=params[3];
		s1=params[4];
		Nc=params[5];
		de=params[6];
	}
	if (p.model==9) {
		d=params[0];
		c1=params[1];
		c2=params[2];
		c3=params[3];
		s1=params[4];
		Ns=params[5];
		de=params[6];
	}
	if (p.model==10) {
		d=params[0];
		c1=params[1];
		c3=params[2];
		de=params[3];
		if (p.usegamma==1) {
			s1=params[4];
		}
	}
	if (p.model==11) {
		d=params[0];
		c1=params[1];
		c2=params[2];
		s1=params[3];
		s2=params[4];
		double Ntemp=params[5]+0.0001;
		Ne=(int)Ntemp;
		de=params[6];
	}
	if (p.model==12) {
		d=params[0];
		c1=params[1];
		c2=params[2];
		c3=params[3];
		s1=params[4];
		double Ntemp=params[5]+0.0001;
		Ne=(int)Ntemp;
		de=params[6];
	}
	if (p.model==13) {
		d=params[0];
		c1=params[1];
		c2=params[2];
		c3=params[3];
		s1=params[4];
		s2=params[5];
		double Ntemp=params[6]+0.0001;
		Ne=(int)Ntemp;
		de=params[7];
	}
	if (p.model==21) {
		d=params[0];
		c1=params[1];
		c2=params[2];
		c3=params[3];
		s1=params[4];
		double Ntemp=params[5]+0.0001;
		Ne=(int)Ntemp;
		de=params[6];
	}
	if (p.model==22) {
		d=params[0];
		c1=params[1];
		c2=params[2];
		s1=params[3];
		s2=params[4];
		double Ntemp=params[5]+0.0001;
		Ne=(int)Ntemp;
		de=params[6];
	}
	if (p.model==23) {
		d=params[0];
		c1=params[1];
		c2=params[2];
		c3=params[3];
		s1=params[4];
		s2=params[5];
		double Ntemp=params[6]+0.0001;
		Ne=(int)Ntemp;
		de=params[7];
	}


}

void ChangeParameters (run_params p, vector<double>& params, vector<double>& params_fix, vector<double> ds, gsl_rng *rgen) {
	int rnd=0;
	rnd=floor(gsl_rng_uniform(rgen)*(params.size()-1))+1;
    if (params_fix[rnd]==0) {
        params[rnd]=params[rnd]+(gsl_rng_uniform(rgen)*2*ds[rnd])-ds[rnd];
    }
    if (p.model==6) {
        params[4]=floor(params[4]);
        if (params[4]<p.dim) {
            params[4]=p.dim;
        }
    }
	if (p.model==11||p.model==12) {
		params[5]=floor(params[5]);
		if (params[5]<p.dim) {
			params[5]=p.dim;
		}
	}
    if (params[rnd]<0) {
        params[rnd]=-params[rnd];
    }
}

double CalcDoublePoisson(double mu, double theta, int y, vector<double> fact_store) {
	double c1=(1-theta)/(12*theta*mu);
	double c2=1+(1/(mu*theta));
	double c=1/(1+(c1*c2));
	double L=log(pow(theta,0.5));
	L=L-(theta*mu);
	if (y>0) {
		L=L-fact_store[y];
		L=L-y;
		L=L+(y*log(y));
	}
	double L2=0;
	if (y>0) {
		L2=1+log(mu);
		L2=L2-log(y);
		L2=L2*(theta*y);
	}
	L=L+L2+log(c);
	return L;
}

double CalcGamma (double& cons, double& t, double& k, double x) {
	double L=pow(x,k-1);
	L=L*exp(-x/t);
	L=L*cons;
	L=log(L);
	return L;
}
