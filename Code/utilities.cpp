#include <vector>
#include <cmath>
#include <iostream>

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


void ChangeParameters (int model, double& s0, double& p1, double& p2, double& de, double ds, gsl_rng *rgen) {
	int r=0;
	if (model==1) {
		r=floor(gsl_rng_uniform(rgen)*3)+1;
	}
	if (model==2) {
		r=floor(gsl_rng_uniform(rgen)*4)+1;
	}
	if (r==1) {
		s0=s0+(gsl_rng_uniform(rgen)*2*ds)-ds;
		if (s0<0) {
			s0=-s0;
		}
	} else if (r==2) {
		p1=p1+(gsl_rng_uniform(rgen)*2*ds)-ds;
		if (p1<0) {
			p1=-p1;
		}
	} else if (r==3) {
		p2=p2+(gsl_rng_uniform(rgen)*2*ds)-ds;
		if (p2<0) {
			p2=-p2;
		}
	} else if (r==4) {
		de=de+(gsl_rng_uniform(rgen)*2*ds)-ds;
		if (de<0) {
			de=-de;
		}
	}
}

void ChangeParametersFix (int model, double lambda, double& s0, double& p1, double& p2, double& de, double ds, gsl_rng *rgen) {
	int r=0;
	if (model==1) {
		r=floor(gsl_rng_uniform(rgen)*2)+1;
	}
	if (model==2) {
		r=floor(gsl_rng_uniform(rgen)*3)+1;
	}
	if (r==1) {
		s0=s0+(gsl_rng_uniform(rgen)*2*ds)-ds;
		if (s0<0) {
			s0=-s0;
		}
		p1=s0/lambda;
	} else if (r==2) {
		p2=p2+(gsl_rng_uniform(rgen)*2*ds)-ds;
		if (p2<0) {
			p2=-p2;
		}
	} else if (r==3) {
		de=de+(gsl_rng_uniform(rgen)*2*ds)-ds;
		if (de<0) {
			de=-de;
		}
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

