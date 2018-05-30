using namespace std;

#include "basicmodel.h"
#include "models.h"
#include "utilities.h"
#include "io.h"

void GetH (run_params p, int model, int dim, double& h, double d, double pe_best, double ps_best, double s0_best, double p1_best, double p2_best, double de_best) {
	cout << "Doing GetH\n";
	//Identify step size that gives an accurate-enough solution of the R-K equations
	h=0.5;
	vector<double> e_vals;
	vector<double> h_vals;
	for (int i=0;i<12;i++) {
		double e_inf=0;
		if (model==1) {
			e_inf=BasicModelRK1 (dim,0,h,d,pe_best,ps_best,s0_best,p1_best,p2_best);
		}
		if (model==2) {
			e_inf=BasicModelRK2 (dim,0,h,d,pe_best,ps_best,s0_best,p1_best,p2_best,de_best);
		}
		
		//cout << h << " " << e_inf << "\n";
		if (e_inf>0) {
			e_vals.push_back(e_inf);
			h_vals.push_back(h);
		}
		h=h/2;
	}
	
	for (int i=0;i<e_vals.size()-1;i++) {
		double r=e_vals[i]/e_vals[i+1];
		r=abs(r-1);
        if (p.precision>0) {
            cout << "Ratio " << r << "\n";
            if (r<pow(10,-p.precision)) {
                cout << "Match " << i << "\n";
                h=h_vals[i];
                break;
            }
        } else {
            if (r<1e-4) {
                cout << "Match " << i << "\n";
                h=h_vals[i];
                break;
            }
        }
	}
	cout << "h= " << h << "\n";
}

void FindLikelihoodRK (run_params p, double& h, double theta, double& lL, int maxv, int verb, int verb_model, int model, int dim, double d, double s0, double p1, double p2, double de, vector<double> cd81_val, vector<double> srb_val, vector< vector<dat> >& all_dat, vector<double> fact_store) {
	
	double e_inf=-10;
	h=h*2;
	while (e_inf==-10) {
		e_inf=BasicModelRK1 (dim,0,h,d,1,1,s0,p1,p2);
		h=h/2;
	}
	
	if (verb_model==1) { //Print progress only for the default case
		if (model==1) {
			double e=BasicModelRK1(dim,verb_model,h,d,1,1,s0,p1,p2);
		}
		if (model==2) {
			double e=BasicModelRK2(dim,verb_model,h,d,1,1,s0,p1,p2,de);
		}
	}
	
	//Calculate values across range of CD81, keeping SRB1
	vector<double> evals_cd81;
	for (int i=0;i<cd81_val.size();i++) {
		double e=-10;
		h=h*2;
		while (e==-10) {
			e=-10;
			if (model==1) {
				e=BasicModelRK1(dim,0,h,d,cd81_val[i],1,s0,p1,p2);
			}
			if (model==2) {
				e=BasicModelRK2(dim,0,h,d,cd81_val[i],1,s0,p1,p2,de);
			}
			h=h/2;
		}
		evals_cd81.push_back(e); //Gives proportion of viruses that get into a cell
	}
	vector<double> evals_srb1;
	
	for (int i=0;i<srb_val.size();i++) {
		double e=-10;
		h=h*2;
		while (e==-10) {
			e=-10;
			if (model==1) {
				e=BasicModelRK1(dim,0,h,d,1,srb_val[i],s0,p1,p2);
			}
			if (model==2) {
				e=BasicModelRK2(dim,0,h,d,1,srb_val[i],s0,p1,p2,de);
			}
			h=h/2;
		}
		evals_srb1.push_back(e);
	}
	//Find probabilities of entry given k input viruses
	vector< vector<double> > pvals_cd81;
	vector< vector<double> > pvals_srb1;
	ConstructPVals(maxv,evals_cd81,evals_srb1,pvals_cd81,pvals_srb1);
	lL=CalcLike(verb,theta,pvals_cd81,pvals_srb1,all_dat,fact_store);
}

void MakeGrid (double& h, int verb_model, int model, int dim, double d, double s0, double p1, double p2, double de) {
    double e_inf=-10;
    h=h*2;
    while (e_inf==-10) {
        e_inf=BasicModelRK1 (dim,0,h,d,1,1,s0,p1,p2);
        h=h/2;
    }
    for (double i=0;i<=2;i=i+0.1) {
        for (double j=0;j<=2;j=j+0.1) {
            double e=0;
            if (model==1) {
                e=BasicModelRK1(dim,verb_model,h,d,i,j,s0,p1,p2);
            }
            if (model==2) {
                e=BasicModelRK2(dim,verb_model,h,d,i,j,s0,p1,p2,de);
            }
            cout << i << " " << j << " " << e << "\n";
        }
    }
}

void ConstructPVals (int maxv, vector<double> evals_cd81, vector<double> evals_srb1, vector< vector<double> >& pvals_cd81, vector< vector<double> >& pvals_srb1) {
	/*cout << "Construct\n";
	for (int i=0;i<evals_cd81.size();i++) {
		cout << evals_cd81[i] << " ";
	}
	cout << "\n";
	for (int i=0;i<evals_srb1.size();i++) {
		cout << evals_srb1[i] << " ";
	}
	cout << "\n";*/
	for (int i=0;i<evals_cd81.size();i++) {
		vector<double> v;
		for (int j=0;j<maxv;j++) {
			double p=1-pow(1-evals_cd81[i],j);
			v.push_back(p);
		}
		pvals_cd81.push_back(v);
	}
	for (int i=0;i<evals_srb1.size();i++) {
		vector<double> v;
		for (int j=0;j<maxv;j++) {
			double p=1-pow(1-evals_srb1[i],j);
			v.push_back(p);
		}
		pvals_srb1.push_back(v);
	}
}

double CalcLike(int verb, double theta, vector< vector<double> >& pvals_cd81, vector< vector<double> >& pvals_srb1, vector< vector<dat> >& all_dat, vector<double>& fact_store) {
	vector<double> cd81_indices;
	vector<double> srb_indices;
	vector<double> cd81_x;
	vector<double> srb_x;
	if (verb>=1) {
		ConstructOutputVectors(all_dat,cd81_indices,srb_indices,cd81_x,srb_x);
	}
	double lL=0;
	for (int i=0;i<all_dat.size();i++) {
		for (int j=0;j<all_dat[i].size();j++) {
			double p=0;
			if (all_dat[i][j].cd81==1) {
				for (int k=0;k<all_dat[i][j].pvir.size();k++) {
					p=p+all_dat[i][j].pvir[k]*pvals_cd81[all_dat[i][j].index][k];
					//cout << i << " " << j << " " << k << " " << all_dat[i][j].pvir[k] << " " << pvals_cd81[all_dat[i][j].index][k] << "\n";
				}
			} else {
				for (int k=0;k<all_dat[i][j].pvir.size();k++) {
					p=p+all_dat[i][j].pvir[k]*pvals_srb1[all_dat[i][j].index][k];
				}
			}
			if (verb>=1) {
				ModifyOutputVectors(i,j,all_dat,p,cd81_indices,srb_indices,cd81_x,srb_x);
			}
			//Multiply by number of cells
			if (p>0) {  //No information if nothing gets through...
				p=p*all_dat[i][j].count;
				double val=CalcDoublePoisson(p,theta,all_dat[i][j].foci,fact_store);
				/*if (val<1e-300) {
					val=1e-300;
				}*/
				//cout << "p " << p << " " << all_dat[i][j].foci << " " << val << " " << lL << "\n";
				//Calculate Poisson likelihood from data
				lL=lL+val;
			}
		}
	}
	if (verb==2) {
		OutputData(all_dat,cd81_x,srb_x);
	}
	if (verb>=1) {
		OutputFreqs(cd81_indices,srb_indices,cd81_x,srb_x);
	}
	//cout << "L " << lL << "\n";
	return lL;
}

double BasicModelRK0 (int dim, int verb, double h, double d, double pe, double e0, double p1, double de) {
	//Model has no SRB1 effect - linear model with downstream entry
	//pe = Proportion of CD81 not bound by antibody - fixed by antibody input
	//ps = Proportion of SRB1 not bound by antibody - fixed by antibody input
	//Variables
	double D=0;  //Dead viruses - get here from every other state at Poisson rate
	
	vector<double> M;
	vector<double> dM;
	vector<double> zero;
	for (int i=0;i<=dim;i++) {
		M.push_back(0);
	}
	dM=M;
	zero=M;
	M[0]=1; //Population starts at M[0][0]; attached to cell, but not bound to anything
	double E=0;
	double pep=pe*e0;
	vector<double> dM_k1=dM;
	vector<double> dM_k2=dM;
	vector<double> dM_k3=dM;
	vector<double> dM_k4=dM;
	double D_k1=0;
	double E_k1=0;
	double D_k2=0;
	double E_k2=0;
	double D_k3=0;
	double E_k3=0;
	double D_k4=0;
	double E_k4=0;
	int i=0;
	int j=0;
	
	double chk=0;
	while (chk<0.999) {
		
		dM=zero;
		dM_k1=dM;
		dM_k2=dM;
		dM_k3=dM;
		dM_k4=dM;
		
		//Step one
		D_k1=0;
		E_k1=0;
		for (i=0;i<=dim;i++) {
			//Non-SRB1-mediated binding of CD81
			if (i<dim) { //Out
				dM_k1[i]=dM_k1[i]-(h*pep*p1*M[i]);
			}
			if (i>0) { //In
				dM_k1[i]=dM_k1[i]-(h*pep*p1*M[i-1]);
			}
			
			//Death from all points on membrane
			D_k1=D_k1+(h*d*M[i]);
			dM_k1[i]=dM_k1[i]-(h*d*M[i]);
			
			//Entry process given sufficient CD81
			if (i==dim) {
				E_k1=E_k1+(h*de*M[i]); //In
				dM_k1[i]=dM_k1[i]-(h*de*M[i]); //Out
			}

		}
		
		//Step two
		D_k2=0;
		E_k2=0;
		for (i=0;i<=dim;i++) {
			//Non-SRB1-mediated binding of CD81
			if (i<dim) {
				dM_k2[i]=dM_k2[i]-(h*pep*p1*(M[i]+dM_k1[i]/2)); //Out
				dM_k2[i+1]=dM_k2[i+1]+(h*pep*p1*(M[i]+dM_k1[i]/2)); //In
			}
			//Death
			D_k2=D_k2+(h*d*(M[i]+(dM_k1[i]/2)));
			dM_k2[i]=dM_k2[i] - (h*d*(M[i]+dM_k1[i]/2));
			
			//Viral entry
			if (i==dim) {
				E_k2=E_k2+(h*de*(M[i]+dM_k1[i]/2));
				dM_k2[i]=dM_k2[i]-(h*de*(M[i]+dM_k1[i]/2));
			}
		}
		
		//Step three
		D_k3=0;
		E_k3=0;
		for (i=0;i<=dim;i++) {
			//Non-SRB1-mediated binding of CD81
			if (i<dim) {
				dM_k3[i]=dM_k3[i]-(h*pep*p1*(M[i]+dM_k2[i]/2)); //Out
				dM_k3[i+1]=dM_k3[i+1]+(h*pep*p1*(M[i]+dM_k2[i]/2)); //In
			}
			//Death
			D_k3=D_k3+(h*d*(M[i]+(dM_k2[i]/2)));
			dM_k3[i]=dM_k3[i] - (h*d*(M[i]+dM_k2[i]/2));
			
			//Viral entry
			if (i==dim) {
				E_k3=E_k3+(h*de*(M[i]+dM_k2[i]/2));
				dM_k3[i]=dM_k3[i]-(h*de*(M[i]+dM_k2[i]/2));
			}
		}
		
		//Step four
		D_k4=0;
		E_k4=0;
		for (i=0;i<=dim;i++) {
			//Non-SRB1-mediated binding of CD81
			if (i<dim) {
				dM_k4[i]=dM_k4[i]-(h*pep*p1*(M[i]+dM_k3[i])); //Out
				dM_k4[i+1]=dM_k4[i+1]+(h*pep*p1*(M[i]+dM_k3[i])); //In
			}
			//Death
			D_k4=D_k4+(h*d*(M[i]+(dM_k3[i])));
			dM_k4[i]=dM_k4[i]-(h*d*(M[i]+dM_k3[i]));
			
			//Viral entry
			if (i==dim) {
				E_k4=E_k4+(h*de*(M[i]+dM_k3[i]));
				dM_k4[i]=dM_k4[i]-(h*de*(M[i]+dM_k3[i]));
			}
		}
		
		//Final step
		D=D+(1./6.)*(D_k1+(2*D_k2)+(2*D_k3)+D_k4);
		for (int i=0;i<=dim;i++) {
			M[i]=M[i]+(1./6.)*(dM_k1[i]+(2*dM_k2[i])+(2*dM_k3[i])+dM_k4[i]);
		}
		//Calculate entry - all viruses that have the required number of CD81 receptors.
		E=E+(1./6.)*(E_k1+(2*E_k2)+(2*E_k3)+E_k4);

		//Calculate proportion of viruses that have gained entry
		double B = 0;
		for (int i=0;i<dim;i++) {
			B=B+M[i];
		}
		int err=0;
		double T=D+B;
		D=D/T;
		if (D<-1e-10) {
			err=1;
			//cout << "Error at D " << D << "\n";
		}
		for (int i=0;i<=dim;i++) {
			M[i]=M[i]/T;
			if (M[i]<-1e-10) {
				//cout << "Error at M[" << i << "][" << j << "]: " << M[i][j] << "\n";
				err=1;
			}
		}
		if (err==1) {
			E=-10;
			return E;
		}
		// cout << "Gen: " << t
		if (verb==1) {
			cout << "D " << D << " M ";
			for (int i=0;i<=dim;i++) {
				cout << M[i] << " ";
			}
			cout << " E " << E << " Total " << D+B << "\n";
		}
		//cout << D << " " << M[0][0] << " " << M[0][1] << " " << E << "\n";
		chk=D+E;
	}
	//cout << "E = " << E << "\n";
	return E;
	
}


double BasicModelRK1 (int dim, int verb, double h, double d, double pe, double ps, double s0, double p1, double p2) {
	//Non-explicit variables in code:
	//pe = Proportion of CD81 not bound by antibody - fixed by antibody input
	//ps = Proportion of SRB1 not bound by antibody - fixed by antibody input
	//Variables
	double D=0;  //Dead viruses - get here from every other state at Poisson rate

	vector< vector<double> > M;
	vector< vector<double> > dM;
	vector< vector<double> > zero;
	for (int i=0;i<=dim;i++) {
		vector<double> mm;
		for (int j=0;j<=dim;j++) {
			mm.push_back(0);
		}
		M.push_back(mm);
	}
	dM=M;
	zero=M;
	M[0][0]=1; //Population starts at M[0][0]; attached to cell, but not bound to anything
	double E=0;
	double pep=pe;
	double psp=ps*s0;
	vector< vector<double> > dM_k1=dM;
	vector< vector<double> > dM_k2=dM;
	vector< vector<double> > dM_k3=dM;
	vector< vector<double> > dM_k4=dM;
	double D_k1=0;
	double E_k1=0;
	double D_k2=0;
	double E_k2=0;
	double D_k3=0;
	double E_k3=0;
	double D_k4=0;
	double E_k4=0;
	int i=0;
	int j=0;

	double chk=0;
	while (chk<0.999) {
		
		dM=zero;
		dM_k1=dM;
		dM_k2=dM;
		dM_k3=dM;
		dM_k4=dM;
		
		//Step one
		D_k1=0;
		E_k1=0;
		for (i=0;i<=dim;i++) {
			for (j=0;j<=dim;j++) {
				//Non-SRB1-mediated binding of CD81
				if (i<dim) { //Out
					dM_k1[i][j]=dM_k1[i][j]-(h*pep*p1*M[i][j]);
				}
				if (i>0) { //In
					dM_k1[i][j]=dM_k1[i][j]-(h*pep*p1*M[i-1][j]);
				}
				
				//SRB1-mediated binding of CD81
				if (i<dim&&j>0) { //Out
					dM_k1[i][j]=dM_k1[i][j]-(h*pep*p2*M[i][j]);
				}
				if (i>0&&j<dim) { //In
					dM_k1[i][j]=dM_k1[i][j]+(h*pep*p2*M[i-1][j+1]);
				}
				
				//Binding to SRB1
				if (j<dim) { //Out
					dM_k1[i][j]=dM_k1[i][j]-(h*psp*M[i][j]);
				}
				if (j>0) { //In
					dM_k1[i][j]=dM_k1[i][j]+(h*psp*M[i][j-1]);
				}
				
				//Death from all points on membrane: by final point have entry
				if (i<dim) {
					D_k1=D_k1+(h*d*M[i][j]);
					dM_k1[i][j]=dM_k1[i][j] - (h*d*M[i][j]);
				}
			}
		}
	
		//Step two
		D_k2=0;
		E_k2=0;
		for (i=0;i<=dim;i++) {
			for (j=0;j<=dim;j++) {
				//Non-SRB1-mediated binding of CD81
				if (i<dim) {
					dM_k2[i][j]=dM_k2[i][j]-(h*pep*p1*(M[i][j]+dM_k1[i][j]/2)); //Out
					dM_k2[i+1][j]=dM_k2[i+1][j]+(h*pep*p1*(M[i][j]+dM_k1[i][j]/2)); //In
				}
				//SRB1-mediated binding of CD81
				if (i<dim&&j>0) {
					dM_k2[i][j]=dM_k2[i][j]-(h*pep*p2*(M[i][j]+dM_k1[i][j]/2)); //Out
					dM_k2[i+1][j-1]=dM_k2[i+1][j-1]+(h*pep*p2*(M[i][j]+dM_k1[i][j]/2)); //In
				}
				//Binding to SRB1
				if (j<dim) {
					dM_k2[i][j]=dM_k2[i][j]-(h*psp*(M[i][j]+dM_k1[i][j]/2)); //Out
					dM_k2[i][j+1]=dM_k2[i][j+1]+(h*psp*(M[i][j]+dM_k1[i][j]/2)); //In
				}
				//Death
				if (i<dim) {
					D_k2=D_k2+(h*d*(M[i][j]+(dM_k1[i][j]/2)));
					dM_k2[i][j]=dM_k2[i][j] - (h*d*(M[i][j]+dM_k1[i][j]/2));
				}

			}
		}
		
		//Step three
		D_k3=0;
		E_k3=0;
		for (i=0;i<=dim;i++) {
			for (j=0;j<=dim;j++) {
				//Non-SRB1-mediated binding of CD81
				if (i<dim) {
					dM_k3[i][j]=dM_k3[i][j]-(h*pep*p1*(M[i][j]+dM_k2[i][j]/2)); //Out
					dM_k3[i+1][j]=dM_k3[i+1][j]+(h*pep*p1*(M[i][j]+dM_k2[i][j]/2)); //In
				}
				//SRB1-mediated binding of CD81
				if (i<dim&&j>0) {
					dM_k3[i][j]=dM_k3[i][j]-(h*pep*p2*(M[i][j]+dM_k2[i][j]/2)); //Out
					dM_k3[i+1][j-1]=dM_k3[i+1][j-1]+(h*pep*p2*(M[i][j]+dM_k2[i][j]/2)); //In
				}
				//Binding to SRB1
				if (j<dim) {
					dM_k3[i][j]=dM_k3[i][j]-(h*psp*(M[i][j]+dM_k2[i][j]/2)); //Out
					dM_k3[i][j+1]=dM_k3[i][j+1]+(h*psp*(M[i][j]+dM_k2[i][j]/2)); //In
				}
				//Death
				if (i<dim) {
					D_k3=D_k3+(h*d*(M[i][j]+(dM_k2[i][j]/2)));
					dM_k3[i][j]=dM_k3[i][j] - (h*d*(M[i][j]+dM_k2[i][j]/2));
				}
				
			}
		}

		//Step four
		D_k4=0;
		E_k4=0;
		for (i=0;i<=dim;i++) {
			for (j=0;j<=dim;j++) {
				//Non-SRB1-mediated binding of CD81
				if (i<dim) {
					dM_k4[i][j]=dM_k4[i][j]-(h*pep*p1*(M[i][j]+dM_k3[i][j])); //Out
					dM_k4[i+1][j]=dM_k4[i+1][j]+(h*pep*p1*(M[i][j]+dM_k3[i][j])); //In
				}
				//SRB1-mediated binding of CD81
				if (i<dim&&j>0) {
					dM_k4[i][j]=dM_k4[i][j]-(h*pep*p2*(M[i][j]+dM_k3[i][j])); //Out
					dM_k4[i+1][j-1]=dM_k4[i+1][j-1]+(h*pep*p2*(M[i][j]+dM_k3[i][j])); //In
				}
				//Binding to SRB1
				if (j<dim) {
					dM_k4[i][j]=dM_k4[i][j]-(h*psp*(M[i][j]+dM_k3[i][j])); //Out
					dM_k4[i][j+1]=dM_k4[i][j+1]+(h*psp*(M[i][j]+dM_k3[i][j])); //In
				}
				//Death
				if (i<dim) {
					D_k4=D_k4+(h*d*(M[i][j]+(dM_k3[i][j])));
					dM_k4[i][j]=dM_k4[i][j]-(h*d*(M[i][j]+dM_k3[i][j]));
				}
			}
		}
		
		//Final step
		D=D+(1./6.)*(D_k1+(2*D_k2)+(2*D_k3)+D_k4);
		for (int i=0;i<=dim;i++) {
			for (int j=0;j<=dim;j++) {
				M[i][j]=M[i][j]+(1./6.)*(dM_k1[i][j]+(2*dM_k2[i][j])+(2*dM_k3[i][j])+dM_k4[i][j]);
			}
		}
		//Calculate entry - all viruses that have the required number of CD81 receptors.
		E=0;
		for (int j=0;j<=dim;j++) {
			E=E+M[dim][j];
		}
		//Calculate proportion of viruses that have gained entry
		double B = 0;
		for (int i=0;i<dim;i++) {
			for (int j=0;j<=dim;j++) {
				B=B+M[i][j];
			}
		}
		int err=0;
		double T=D+B;
		D=D/T;
		if (D<-1e-10) {
			err=1;
			//cout << "Error at D " << D << "\n";
		}
		for (int i=0;i<=dim;i++) {
			for (int j=0;j<=dim;j++) {
				M[i][j]=M[i][j]/T;
				if (M[i][j]<-1e-10) {
					//cout << "Error at M[" << i << "][" << j << "]: " << M[i][j] << "\n";
					err=1;
				}
			}
		}
		if (err==1) {
			E=-10;
			return E;
		}
		// cout << "Gen: " << t
		if (verb==1) {
			cout << "D " << D << " M ";
			for (int i=0;i<=dim;i++) {
				for (int j=0;j<=dim;j++) {
					cout << M[i][j] << " ";
				}
			}
			cout << " E " << E << " Total " << D+B << "\n";
		}
		//cout << D << " " << M[0][0] << " " << M[0][1] << " " << E << "\n";
		chk=D+E;
	}
	//cout << "E = " << E << "\n";
	return E;
	
}



double BasicModelRK2 (int dim, int verb, double h, double d, double pe, double ps, double s0, double p1, double p2, double de) {
	//Non-explicit variables in code:
	//pe = Proportion of CD81 not bound by antibody - fixed by antibody input
	//ps = Proportion of SRB1 not bound by antibody - fixed by antibody input
	//Variables
	double D=0;  //Dead viruses - get here from every other state at Poisson rate
	
	vector< vector<double> > M;
	vector< vector<double> > dM;
	vector< vector<double> > zero;
	for (int i=0;i<=dim;i++) {
		vector<double> mm;
		for (int j=0;j<=dim;j++) {
			mm.push_back(0);
		}
		M.push_back(mm);
	}
	dM=M;
	zero=M;
	M[0][0]=1; //Population starts at M[0][0]; attached to cell, but not bound to anything
	double E=0;
	double pep=pe;
	double psp=ps*s0;
	vector< vector<double> > dM_k1=dM;
	vector< vector<double> > dM_k2=dM;
	vector< vector<double> > dM_k3=dM;
	vector< vector<double> > dM_k4=dM;
	double D_k1=0;
	double E_k1=0;
	double D_k2=0;
	double E_k2=0;
	double D_k3=0;
	double E_k3=0;
	double D_k4=0;
	double E_k4=0;
	int i=0;
	int j=0;
	
	double chk=0;
	while (chk<0.999) {
		
		dM=zero;
		dM_k1=dM;
		dM_k2=dM;
		dM_k3=dM;
		dM_k4=dM;
		
		//Step one
		D_k1=0;
		E_k1=0;
		for (i=0;i<=dim;i++) {
			for (j=0;j<=dim;j++) {
				//Non-SRB1-mediated binding of CD81
				if (i<dim) { //Out
					dM_k1[i][j]=dM_k1[i][j]-(h*pep*p1*M[i][j]);
				}
				if (i>0) { //In
					dM_k1[i][j]=dM_k1[i][j]-(h*pep*p1*M[i-1][j]);
				}
				
				//SRB1-mediated binding of CD81
				if (i<dim&&j>0) { //Out
					dM_k1[i][j]=dM_k1[i][j]-(h*pep*p2*M[i][j]);
				}
				if (i>0&&j<dim) { //In
					dM_k1[i][j]=dM_k1[i][j]+(h*pep*p2*M[i-1][j+1]);
				}
				
				//Binding to SRB1
				if (j<dim) { //Out
					dM_k1[i][j]=dM_k1[i][j]-(h*psp*M[i][j]);
				}
				if (j>0) { //In
					dM_k1[i][j]=dM_k1[i][j]+(h*psp*M[i][j-1]);
				}
				
				//Death from all points on membrane
				D_k1=D_k1+(h*d*M[i][j]);
				dM_k1[i][j]=dM_k1[i][j] - (h*d*M[i][j]);
				
				//Entry process given sufficient CD81
				if (i==dim) {
					E_k1=E_k1+(h*de*M[i][j]); //In
					dM_k1[i][j]=dM_k1[i][j]-(h*de*M[i][j]); //Out
				}
				
			}
		}
		
		//Step two
		D_k2=0;
		E_k2=0;
		for (i=0;i<=dim;i++) {
			for (j=0;j<=dim;j++) {
				//Non-SRB1-mediated binding of CD81
				if (i<dim) {
					dM_k2[i][j]=dM_k2[i][j]-(h*pep*p1*(M[i][j]+dM_k1[i][j]/2)); //Out
					dM_k2[i+1][j]=dM_k2[i+1][j]+(h*pep*p1*(M[i][j]+dM_k1[i][j]/2)); //In
				}
				//SRB1-mediated binding of CD81
				if (i<dim&&j>0) {
					dM_k2[i][j]=dM_k2[i][j]-(h*pep*p2*(M[i][j]+dM_k1[i][j]/2)); //Out
					dM_k2[i+1][j-1]=dM_k2[i+1][j-1]+(h*pep*p2*(M[i][j]+dM_k1[i][j]/2)); //In
				}
				//Binding to SRB1
				if (j<dim) {
					dM_k2[i][j]=dM_k2[i][j]-(h*psp*(M[i][j]+dM_k1[i][j]/2)); //Out
					dM_k2[i][j+1]=dM_k2[i][j+1]+(h*psp*(M[i][j]+dM_k1[i][j]/2)); //In
				}
				//Death
				D_k2=D_k2+(h*d*(M[i][j]+(dM_k1[i][j]/2)));
				dM_k2[i][j]=dM_k2[i][j] - (h*d*(M[i][j]+dM_k1[i][j]/2));
				
				//Viral entry
				if (i==dim) {
					E_k2=E_k2+(h*de*(M[i][j]+dM_k1[i][j]/2));
					dM_k2[i][j]=dM_k2[i][j]-(h*de*(M[i][j]+dM_k1[i][j]/2));
				}

			}
		}
		
		//Step three
		D_k3=0;
		E_k3=0;
		for (i=0;i<=dim;i++) {
			for (j=0;j<=dim;j++) {
				//Non-SRB1-mediated binding of CD81
				if (i<dim) {
					dM_k3[i][j]=dM_k3[i][j]-(h*pep*p1*(M[i][j]+dM_k2[i][j]/2)); //Out
					dM_k3[i+1][j]=dM_k3[i+1][j]+(h*pep*p1*(M[i][j]+dM_k2[i][j]/2)); //In
				}
				//SRB1-mediated binding of CD81
				if (i<dim&&j>0) {
					dM_k3[i][j]=dM_k3[i][j]-(h*pep*p2*(M[i][j]+dM_k2[i][j]/2)); //Out
					dM_k3[i+1][j-1]=dM_k3[i+1][j-1]+(h*pep*p2*(M[i][j]+dM_k2[i][j]/2)); //In
				}
				//Binding to SRB1
				if (j<dim) {
					dM_k3[i][j]=dM_k3[i][j]-(h*psp*(M[i][j]+dM_k2[i][j]/2)); //Out
					dM_k3[i][j+1]=dM_k3[i][j+1]+(h*psp*(M[i][j]+dM_k2[i][j]/2)); //In
				}
				//Death
				D_k3=D_k3+(h*d*(M[i][j]+(dM_k2[i][j]/2)));
				dM_k3[i][j]=dM_k3[i][j] - (h*d*(M[i][j]+dM_k2[i][j]/2));

				//Viral entry
				if (i==dim) {
					E_k3=E_k3+(h*de*(M[i][j]+dM_k2[i][j]/2));
					dM_k3[i][j]=dM_k3[i][j]-(h*de*(M[i][j]+dM_k2[i][j]/2));
				}

				
			}
		}
		
		//Step four
		D_k4=0;
		E_k4=0;
		for (i=0;i<=dim;i++) {
			for (j=0;j<=dim;j++) {
				//Non-SRB1-mediated binding of CD81
				if (i<dim) {
					dM_k4[i][j]=dM_k4[i][j]-(h*pep*p1*(M[i][j]+dM_k3[i][j])); //Out
					dM_k4[i+1][j]=dM_k4[i+1][j]+(h*pep*p1*(M[i][j]+dM_k3[i][j])); //In
				}
				//SRB1-mediated binding of CD81
				if (i<dim&&j>0) {
					dM_k4[i][j]=dM_k4[i][j]-(h*pep*p2*(M[i][j]+dM_k3[i][j])); //Out
					dM_k4[i+1][j-1]=dM_k4[i+1][j-1]+(h*pep*p2*(M[i][j]+dM_k3[i][j])); //In
				}
				//Binding to SRB1
				if (j<dim) {
					dM_k4[i][j]=dM_k4[i][j]-(h*psp*(M[i][j]+dM_k3[i][j])); //Out
					dM_k4[i][j+1]=dM_k4[i][j+1]+(h*psp*(M[i][j]+dM_k3[i][j])); //In
				}
				//Death
				D_k4=D_k4+(h*d*(M[i][j]+(dM_k3[i][j])));
				dM_k4[i][j]=dM_k4[i][j]-(h*d*(M[i][j]+dM_k3[i][j]));

				//Viral entry
				if (i==dim) {
					E_k4=E_k4+(h*de*(M[i][j]+dM_k3[i][j]));
					dM_k4[i][j]=dM_k4[i][j]-(h*de*(M[i][j]+dM_k3[i][j]));
				}

			}
		}
		
		//Final step
		D=D+(1./6.)*(D_k1+(2*D_k2)+(2*D_k3)+D_k4);
		for (int i=0;i<=dim;i++) {
			for (int j=0;j<=dim;j++) {
				M[i][j]=M[i][j]+(1./6.)*(dM_k1[i][j]+(2*dM_k2[i][j])+(2*dM_k3[i][j])+dM_k4[i][j]);
			}
		}
		//Calculate entry - all viruses that have the required number of CD81 receptors.
		E=E+(1./6.)*(E_k1+(2*E_k2)+(2*E_k3)+E_k4);

		//Calculate proportion of viruses that have gained entry
		
		//Proportion of viruses attached to the membrane
		double B = 0;
		for (int i=0;i<=dim;i++) {
			for (int j=0;j<=dim;j++) {
				B=B+M[i][j];
			}
		}
		
		int err=0;
		//Correct for machine precision
		double T=D+B+E;
		D=D/T;
		E=E/T;
		B=B/T;
		for (int i=0;i<=dim;i++) {
			for (int j=0;j<=dim;j++) {
				M[i][j]=M[i][j]/T;
				if (M[i][j]<-1e-10) {
					err=1;
				}
			}
		}
		if (D<-1e-10) {
			err=1;
		}
		if (err==1) {
			E=-10;
			return E;
		}
		
	//	if (verb==1) {
	//		cout << D << " " << E << " " << B << " " << D+E+B << "\n";
	//	}

		if (verb==1) {
			cout << D << " " << E << " " << B << " " << D+E << " ";
			for (int i=0;i<M.size();i++) {
				for (int j=0;j<M[i].size();j++) {
					cout << M[i][j] << " ";
				}
			}
			cout << "\n";
		 }

		chk=D+E;
	}
	//cout << "E = " << E << "\n";
	return E;
	
}




