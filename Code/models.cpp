#include "basicmodel.h"
#include "models.h"
#include "utilities.h"
#include "io.h"


void GetH (run_params p, int dim, double& h, double pe_best, double ps_best, vector<double>& params) {
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
    cout << "Collected Ne is " << Ne << "\n";
    
	cout << "Doing GetH\n";
	//Identify step size that gives an accurate-enough solution of the R-K equations
	h=32;
	vector<double> e_vals;
	vector<double> h_vals;
	for (int i=0;i<18;i++) {
        double e_inf=-10;
		if (p.model==0) {
			e_inf=BasicModelRK0(dim,0,1,h,pe_best,ps_best,d,c1,Ne,de);
        } else if (p.model==6) {
            e_inf=BasicModelRK6(dim,0,1,h,pe_best,ps_best,d,Ne,c1,c2,s1,de);
		} else if (p.model==21) {
			e_inf=BasicModelRK21(dim,0,1,h,pe_best,ps_best,d,Ne,c1,c2,c3,s1,de);
		} else if (p.model==22) {
			e_inf=BasicModelRK22(dim,0,1,h,pe_best,ps_best,d,Ne,c1,c2,s1,s2,de);
		}
		cout << h << " " << e_inf << "\n";
		if (e_inf>0) {
			e_vals.push_back(e_inf);
			h_vals.push_back(h);
		}
		h=h/2;
	}
	
	for (int i=0;i<e_vals.size()-1;i++) {
		double rat=e_vals[i]/e_vals[i+1];
		rat=abs(rat-1);
        if (p.precision>0) {
            cout << "Ratio " << rat << "\n";
            if (rat<pow(10,-p.precision)) {
                cout << "Match " << i << "\n";
                h=h_vals[i];
                break;
            }
        } else {
            if (rat<1e-4) {
                cout << "Match " << i << "\n";
                h=h_vals[i];
                break;
            }
        }
	}
	cout << "h= " << h << "\n";
}

void CheckH (run_params p, int dim, double& h, double& h_out, double pe_best, double ps_best, vector<double>& params) {
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
    
    //Identify step size that gives an accurate-enough solution of the R-K equations
    h=32;
    vector<double> e_vals;
    vector<double> h_vals;
    int index=0;
    for (int i=0;i<18;i++) {
        double e_inf=-10;
		if (p.model==0) {
			e_inf=BasicModelRK0(dim,0,1,h,1,1,d,c1,Ne,de);
        } else if (p.model==6) {
            e_inf=BasicModelRK6(dim,0,1,h,1,1,d,Ne,c1,c2,s1,de);
		} else if (p.model==21) {
			e_inf=BasicModelRK21(dim,0,1,h,1,1,d,Ne,c1,c2,c3,s1,de);
		} else if (p.model==22) {
			e_inf=BasicModelRK22(dim,0,1,h,1,1,d,Ne,c1,c2,s1,s2,de);
		}
       // cout << i << " " << h << " " << e_inf << "\n";
        if (e_inf>0) {
            e_vals.push_back(e_inf);
            h_vals.push_back(h);
            index++;
        }
        if (index>1) {
            double rat=e_vals[index-2]/e_vals[index-1];
            rat=abs(rat-1);
            if (rat<1e-4) {
                h=h_vals[index-2];
                break;
            }
        }
        
        h=h/2;
    }
    h_out=h;
    //cout << "H " << h_in << " " << h_out << "\n";
}

double BasicModelRK0 (int dim, int verb, int quickprop, double h, double pe, double ps, double d, double c1, int Ne, double de) {

	//Set up initial variables
	double D=0;  //Dead viruses - get here from every other state at Poisson rate
	double E=0;	 //Viruses which have gained entry
	vector<double> M;  //In this model we just have one row, counting number of CD81
	vector<double> dM;
	vector<double> zero;
	MatrixSetup(dim,M,dM,zero);
	int r_check=0; //Check for having reached a state where only the final column is occupied
	vector<double> rates;
	for (int i=0;i<=Ne;i++) {
		double r=(i+0.)/(Ne+0.);
		rates.push_back(r);
	}

	
	int rs=0;  //Tracking for row start
	
	M[0]=1; //Population starts at M[0]; attached to cell, but not bound to anything

	//Set up Runge-Kutta variables
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
	double pec1=pe*c1;
	double chk=0;
	double cut=1e-20;
	
	while (chk<0.999) {
		
		//Proportion of viruses attached to the membrane
		double B = 0;
		for (int i=rs;i<M.size();i++) {
			B=B+M[i];
		}
		
		if (rs==M.size()-1) {
			//cout << "Only one column left\n";
			//At this point, don't need to model the complete matrix any more.  Only outputs are to E and D
			r_check=1;
			if (quickprop==1) {
				FinalCalc (rs,d,de,D,E,M);
				return E;
			}
		}
		
		dM=zero;
		dM_k1=dM;
		dM_k2=dM;
		dM_k3=dM;
		dM_k4=dM;
		
		//Step one
		D_k1=0;
		E_k1=0;
		for (int i=rs;i<M.size();i++) {
			double dm=0;
			double k1=0;
			if (Ne-i>=0) {
				k1=rates[Ne-i];
			}
			if (i<M.size()-1) {  //Gain of CD81
				dm=h*(pec1*k1*M[i]);
				dM_k1[i]=dM_k1[i]-dm;
				dM_k1[i+1]=dM_k1[i+1]+dm;
			}
			//Death from membrane
			dm=h*d*M[i];
			D_k1=D_k1+dm;
			dM_k1[i]=dM_k1[i]-dm;
			//Viral entry
			if (i==dim) {
				dm=h*de*M[i];
				E_k1=E_k1+dm; //In
				dM_k1[i]=dM_k1[i]-dm; //Out
			}
		}
		
		//Step two
		D_k2=0;
		E_k2=0;
		vector<double> temp;
		temp=dM_k1;
		HalfTemp(temp);
		for (int i=rs;i<M.size();i++) {
			double dm=0;
			double k1=0;
			if (Ne-i>=0) {
				k1=rates[Ne-i];
			}
			if (i<M.size()-1) {  //Gain of CD81
				dm=h*(pec1*k1*(M[i]+temp[i]));
				dM_k2[i]=dM_k2[i]-dm;
				dM_k2[i+1]=dM_k2[i+1]+dm;
			}
			//Death from membrane
			dm=h*d*(M[i]+temp[i]);
			D_k2=D_k2+dm;
			dM_k2[i]=dM_k2[i]-dm;
			//Viral entry
			if (i==dim) {
				dm=h*de*(M[i]+temp[i]);
				E_k2=E_k2+dm; //In
				dM_k2[i]=dM_k2[i]-dm; //Out
			}
		}
		
		//Step three
		D_k3=0;
		E_k3=0;
		temp=dM_k2;
		HalfTemp(temp);
		for (int i=rs;i<M.size();i++) {
			double dm=0;
			double k1=0;
			if (Ne-i>=0) {
				k1=rates[Ne-i];
			}
			if (i<M.size()-1) {  //Gain of CD81
				dm=h*(pec1*k1*(M[i]+temp[i]));
				dM_k3[i]=dM_k3[i]-dm;
				dM_k3[i+1]=dM_k3[i+1]+dm;
			}
			//Death from membrane
			dm=h*d*(M[i]+temp[i]);
			D_k3=D_k3+dm;
			dM_k3[i]=dM_k3[i]-dm;
			//Viral entry
			if (i==dim) {
				dm=h*de*(M[i]+temp[i]);
				E_k3=E_k3+dm; //In
				dM_k3[i]=dM_k3[i]-dm; //Out
			}
		}
		
		
		//Step four
		D_k4=0;
		E_k4=0;
		temp=dM_k3;
		for (int i=rs;i<M.size();i++) {
			double dm=0;
			double k1=0;
			if (Ne-i>=0) {
				k1=rates[Ne-i];
			}
			if (i<M.size()-1) {  //Gain of CD81
				dm=h*(pec1*k1*(M[i]+temp[i]));
				dM_k4[i]=dM_k4[i]-dm;
				dM_k4[i+1]=dM_k4[i+1]+dm;
			}
			//Death from membrane
			dm=h*d*(M[i]+temp[i]);
			D_k4=D_k4+dm;
			dM_k4[i]=dM_k4[i]-dm;
			//Viral entry
			if (i==dim) {
				dm=h*de*(M[i]+temp[i]);
				E_k4=E_k4+dm; //In
				dM_k4[i]=dM_k4[i]-dm; //Out
			}
		}
		
		//Final step
		D=D+(1./6.)*(D_k1+(2*D_k2)+(2*D_k3)+D_k4);
		//Calculate entry - all viruses that have the required number of CD81 receptors.
		E=E+(1./6.)*(E_k1+(2*E_k2)+(2*E_k3)+E_k4);
		for (int i=rs;i<M.size();i++) {
			M[i]=M[i]+(1./6.)*(dM_k1[i]+(2*dM_k2[i])+(2*dM_k3[i])+dM_k4[i]);
		}
		
		//Calculate proportion of viruses that have gained entry
		
		//Proportion of viruses attached to the membrane
		B = 0;
		for (int i=rs;i<M.size();i++) {
			B=B+M[i];
		}
		//cout << "B " << B << "\n";
		
		int err=0;
		//Correct for machine precision errors
		double T=D+B+E;
		//cout << "T " << T << "\n";
		D=D/T;
		E=E/T;
		B=B/T;
		for (int i=rs;i<M.size();i++) {
			M[i]=M[i]/T;
			if (M[i]<-1e-10) {
				err=1;
			}
		}
		if (D<-1e-10) {
			err=1;
		}
		if (err==1) {
			E=-10;
			return E;
		}
		
		if (verb==1) {
			cout << D << " " << E << " " << B << " " << D+E << "\n";
			for (int i=0;i<M.size();i++) {
				cout << M[i] << " ";
			}
			cout << "\n";
		}
		
		chk=D+E;
		
		//Update the value of rs - exclude columns once they are cleaned up
		if (M[rs]<1e-20) {
			M[rs+1]=M[rs+1]+M[rs];
			rs++;
		}
	}
	UpdateColumns(rs,cut,M);
	return E;
	
}

double BasicModelRK6 (int dim, int verb, int quickprop, double h, double pe, double ps, double d, int Ne, double c1, double c2, double s1, double de) {
	
    //Set up initial variables
    double D=0;  //Dead viruses - get here from every other state at Poisson rate
    double E=0;	 //Viruses which have gained entry
    vector< vector<double> > M;  //From M[0][0] up.  Question - how far to go with SRB1?
    vector< vector<double> > dM;
    vector< vector<double> > zero;
    vector< vector<int> > mask;
    MatrixSetup6(dim,Ne,M,dM,zero,mask);
	vector<double> rates;
	for (int i=0;i<=Ne;i++) {
		double r=(i+0.)/(Ne+0.);
		rates.push_back(r);
	}
	
    //Row start
    int rs=0;
	
    //Column start vector
    vector<int> cs;
    for (int i=0;i<M.size();i++) {
        cs.push_back(0);
    }
	
    M[0][0]=1; //Population starts at M[0][0]; attached to cell, but not bound to anything
	
    double cut=1e-20;
    //MakeMask(cut,M,mask);
	
    //Set up Runge-Kutta variables
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
    double pec1=pe*c1;
    double pec2=pe*c2;
    double pss1=ps*s1;
    double chk=0;
	
    int domask=1;
	
	
    while (chk<0.999) {
		
		
        //Proportion of viruses attached to the membrane
        double B = 0;
        for (int i=rs;i<M.size();i++) {
            for (int j=cs[i];j<M[i].size();j++) {
                B=B+M[i][j];
            }
        }
        
        if (rs==M.size()-1) {
            //cout << "Only one column left\n";
            //At this point, don't need to model the complete matrix any more.  Only outputs are to E and D
            if (quickprop==1) {
                FinalCalcMat (rs,cs,d,de,D,E,M);
                return E;
            }
        }
        
        //Mask elements - don't need to include all of them in the calcuation
        if (domask==1) {
            MakeMask(cut,M,mask);
        }
        //PrintMatrix(mask);
        
        dM=zero;
        dM_k1=dM;
        dM_k2=dM;
        dM_k3=dM;
        dM_k4=dM;
        
        //Note: Seems to be some problem in the differential equations - T is getting too high.
        
        //Step one
        D_k1=0;
        E_k1=0;
        for (int i=rs;i<M.size();i++) {
            for (int j=cs[i];j<M[i].size();j++) {
                if (domask==0||mask[i][j]>0) {
                    double dm=0;
					double k1=0;
					if (Ne-i-j>=0) {
						k1=rates[Ne-i-j];
					}
					double k2=rates[j];;
                    if (i<M.size()-1) { //SR-B1 independent binding
                        if (j<Ne) {
                            dm=h*pec1*k1*M[i][j];
                            dM_k1[i][j]=dM_k1[i][j]-dm;
                            dM_k1[i+1][j]=dM_k1[i+1][j]+dm;
                        }
                        if (j>0) { //SR-B1 dependent binding
                            dm=h*pec2*k2*M[i][j];
                            dM_k1[i][j]=dM_k1[i][j]-dm;
                            dM_k1[i+1][j-1]=dM_k1[i+1][j-1]+dm;
                        }
                    }
                    if (j<M[i].size()-1) {  //Binding of SR-B1
                        dm=h*pss1*k1*M[i][j];
                        dM_k1[i][j]=dM_k1[i][j]-dm;
                        dM_k1[i][j+1]=dM_k1[i][j+1]+dm;
                    }
                    
                    //Death from all points on membrane
                    dm=h*d*M[i][j];
                    D_k1=D_k1+dm;
                    dM_k1[i][j]=dM_k1[i][j]-dm;
                    
                    //Entry process given sufficient CD81
                    if (i==dim) {
                        dm=h*de*M[i][j];
                        E_k1=E_k1+dm; //In
                        dM_k1[i][j]=dM_k1[i][j]-dm; //Out
                    }
                }
            }
        }
        
        //Model is coded to here...
        
        //Step two
        D_k2=0;
        E_k2=0;
        vector< vector<double> > temp;
        temp=dM_k1;
        HalfTempMat(temp);
        for (int i=rs;i<M.size();i++) {
            for (int j=cs[i];j<M[i].size();j++) {
                if (domask==0||mask[i][j]>0) {
                    double dm=0;
					double k1=0;
					if (Ne-i-j>=0) {
						k1=rates[Ne-i-j];
					}
					double k2=rates[j];;
                    if (i<M.size()-1) { //SR-B1 independent binding
                        if (j<Ne) {
                            dm=h*pec1*k1*(M[i][j]+temp[i][j]);
                            dM_k2[i][j]=dM_k2[i][j]-dm;
                            dM_k2[i+1][j]=dM_k2[i+1][j]+dm;
                        }
                        if (j>0) { //SR-B1 dependent binding
                            dm=h*pec2*k2*(M[i][j]+temp[i][j]);
                            dM_k2[i][j]=dM_k2[i][j]-dm;
                            dM_k2[i+1][j-1]=dM_k2[i+1][j-1]+dm;
                        }
                    }
                    if (j<M[i].size()-1) {  //Binding of SR-B1
                        dm=h*pss1*k1*(M[i][j]+temp[i][j]);
                        dM_k2[i][j]=dM_k2[i][j]-dm;
                        dM_k2[i][j+1]=dM_k2[i][j+1]+dm;
                    }
                    
                    //Death from all points on membrane
                    dm=h*d*(M[i][j]+temp[i][j]);
                    D_k2=D_k2+dm;
                    dM_k2[i][j]=dM_k2[i][j]-dm;
                    
                    //Entry process given sufficient CD81
                    if (i==dim) {
                        dm=h*de*(M[i][j]+temp[i][j]);
                        E_k2=E_k2+dm; //In
                        dM_k2[i][j]=dM_k2[i][j]-dm; //Out
                    }
                }
            }
        }

        //Step three
        D_k3=0;
        E_k3=0;
        temp=dM_k2;
        HalfTempMat(temp);
        for (int i=rs;i<M.size();i++) {
            for (int j=cs[i];j<M[i].size();j++) {
                if (domask==0||mask[i][j]>0) {
                    double dm=0;
					double k1=0;
					if (Ne-i-j>=0) {
						k1=rates[Ne-i-j];
					}
					double k2=rates[j];;
                    if (i<M.size()-1) { //SR-B1 independent binding
                        if (j<Ne) {
                            dm=h*pec1*k1*(M[i][j]+temp[i][j]);
                            dM_k3[i][j]=dM_k3[i][j]-dm;
                            dM_k3[i+1][j]=dM_k3[i+1][j]+dm;
                        }
                        if (j>0) { //SR-B1 dependent binding
                            dm=h*pec2*k2*(M[i][j]+temp[i][j]);
                            dM_k3[i][j]=dM_k3[i][j]-dm;
                            dM_k3[i+1][j-1]=dM_k3[i+1][j-1]+dm;
                        }
                    }
                    if (j<M[i].size()-1) {  //Binding of SR-B1
                        dm=h*pss1*k1*(M[i][j]+temp[i][j]);
                        dM_k3[i][j]=dM_k3[i][j]-dm;
                        dM_k3[i][j+1]=dM_k3[i][j+1]+dm;
                    }
                    
                    //Death from all points on membrane
                    dm=h*d*(M[i][j]+temp[i][j]);
                    D_k3=D_k3+dm;
                    dM_k3[i][j]=dM_k3[i][j]-dm;
                    
                    //Entry process given sufficient CD81
                    if (i==dim) {
                        dm=h*de*(M[i][j]+temp[i][j]);
                        E_k3=E_k3+dm; //In
                        dM_k3[i][j]=dM_k3[i][j]-dm; //Out
                    }
                }
            }
        }

        //Step four
        D_k4=0;
        E_k4=0;
        temp=dM_k3;
        for (int i=rs;i<M.size();i++) {
            for (int j=cs[i];j<M[i].size();j++) {
                if (domask==0||mask[i][j]>0) {
                    double dm=0;
					double k1=0;
					if (Ne-i-j>=0) {
						k1=rates[Ne-i-j];
					}
					double k2=rates[j];;
                    if (i<M.size()-1) { //SR-B1 independent binding
                        if (j<Ne) {
                            dm=h*pec1*k1*(M[i][j]+temp[i][j]);
                            dM_k4[i][j]=dM_k4[i][j]-dm;
                            dM_k4[i+1][j]=dM_k4[i+1][j]+dm;
                        }
                        if (j>0) { //SR-B1 dependent binding
                            dm=h*pec2*k2*(M[i][j]+temp[i][j]);
                            dM_k4[i][j]=dM_k4[i][j]-dm;
                            dM_k4[i+1][j-1]=dM_k4[i+1][j-1]+dm;
                        }
                    }
                    if (j<M[i].size()-1) {  //Binding of SR-B1
                        dm=h*pss1*k1*(M[i][j]+temp[i][j]);
                        dM_k4[i][j]=dM_k4[i][j]-dm;
                        dM_k4[i][j+1]=dM_k4[i][j+1]+dm;
                    }
                    
                    //Death from all points on membrane
                    dm=h*d*(M[i][j]+temp[i][j]);
                    D_k4=D_k4+dm;
                    dM_k4[i][j]=dM_k4[i][j]-dm;
                    
                    //Entry process given sufficient CD81
                    if (i==dim) {
                        dm=h*de*(M[i][j]+temp[i][j]);
                        E_k4=E_k4+dm; //In
                        dM_k4[i][j]=dM_k4[i][j]-dm; //Out
                    }
                }
            }
        }

        
        //Final step
        D=D+(1./6.)*(D_k1+(2*D_k2)+(2*D_k3)+D_k4);
        for (int i=rs;i<M.size();i++) {
            for (int j=cs[i];j<M[i].size();j++) {
                M[i][j]=M[i][j]+(1./6.)*(dM_k1[i][j]+(2*dM_k2[i][j])+(2*dM_k3[i][j])+dM_k4[i][j]);
            }
        }
        //Calculate entry - all viruses that have the required number of CD81 receptors.
        E=E+(1./6.)*(E_k1+(2*E_k2)+(2*E_k3)+E_k4);
        
        //Calculate proportion of viruses that have gained entry
        
        //Proportion of viruses attached to the membrane
        B = 0;
        for (int i=rs;i<M.size();i++) {
            for (int j=cs[i];j<M[i].size();j++) {
                B=B+M[i][j];
            }
        }
        //cout << "B " << B << "\n";
        
        int err=0;
        //Correct for machine precision errors
        double T=D+B+E;
        //cout << "T " << T << "\n";
        D=D/T;
        E=E/T;
        B=B/T;
        for (int i=rs;i<M.size();i++) {
            for (int j=cs[i];j<M[i].size();j++) {
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
        
        if (verb==1) {
            cout << D << " " << E << " " << B << " " << D+E << " ";
            for (int i=0;i<M.size();i++) {
                for (int j=0;j<M[i].size();j++) {
                    cout << M[i][j] << " ";
                }
            }
            cout << "\n";
        }
        //cout << "D " << D << "\n";
        //cout << "E " << E << "\n";
        
        chk=D+E;
        
        int change=0;
        //Clear up the matrix M.  Remove things less than 1e-20.
        CleanupMatrix (change,cut,Ne,pec1,pec2,pss1,cs,M);
        
        //Update the value of rs - exclude columns once they are cleaned up
        UpdateColumnsMat (change,rs,cs,M);
    }
    //cout << "E = " << E << "\n";
    return E;
    
}



double BasicModelRK21 (int dim, int verb, int quickprop, double h, double pe, double ps, double d, int Ne, double c1, double c2, double c3, double s1, double de) {
	
	//Set up initial variables
	double D=0;  //Dead viruses - get here from every other state at Poisson rate
	double E=0;	 //Viruses which have gained entry
	vector< vector<double> > M;  //From M[0][0] up.  Question - how far to go with SRB1?
	vector< vector<double> > dM;
	vector< vector<double> > zero;
	vector< vector<int> > mask;
	MatrixSetup6(dim,Ne,M,dM,zero,mask);
	vector<double> rates;
	for (int i=0;i<=Ne;i++) {
		double r=(i+0.)/(Ne+0.);
		rates.push_back(r);
	}
	
	//Row start
	int rs=0;
	
	//Column start vector
	vector<int> cs;
	for (int i=0;i<M.size();i++) {
		cs.push_back(0);
	}
	
	M[0][0]=1; //Population starts at M[0][0]; attached to cell, but not bound to anything
	
	double cut=1e-20;
	//MakeMask(cut,M,mask);
	
	//Set up Runge-Kutta variables
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
	double pec1=pe*c1;
	double pec2=pe*c2;
	double pss1=ps*s1;
	double chk=0;
	
	int domask=1;
	
	
	while (chk<0.999) {
		
		
		//Proportion of viruses attached to the membrane
		double B = 0;
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
				B=B+M[i][j];
			}
		}
		
		if (rs==M.size()-1) {
			//cout << "Only one column left\n";
			//At this point, don't need to model the complete matrix any more.  Only outputs are to E and D
			if (quickprop==1) {
				FinalCalcMat (rs,cs,d,de,D,E,M);
				return E;
			}
		}
		
		//Mask elements - don't need to include all of them in the calcuation
		if (domask==1) {
			MakeMask(cut,M,mask);
		}
		//PrintMatrix(mask);
		
		dM=zero;
		dM_k1=dM;
		dM_k2=dM;
		dM_k3=dM;
		dM_k4=dM;
		
		//Note: Seems to be some problem in the differential equations - T is getting too high.
		
		//Step one
		D_k1=0;
		E_k1=0;
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
				if (domask==0||mask[i][j]>0) {
					double dm=0;
					double k1=0;
					if (Ne-i-j>=0) {
						k1=rates[Ne-i-j];
					}
					double k2=rates[j];;
					if (i<M.size()-1) { //SR-B1 independent binding
						if (j<Ne) {
							dm=h*((pec1*k1*M[i][j])-(c3*(i+1)*M[i+1][j]));
							dM_k1[i][j]=dM_k1[i][j]-dm;
							dM_k1[i+1][j]=dM_k1[i+1][j]+dm;
						}
						if (j>0) { //SR-B1 dependent binding
							dm=h*pec2*k2*M[i][j];
							dM_k1[i][j]=dM_k1[i][j]-dm;
							dM_k1[i+1][j-1]=dM_k1[i+1][j-1]+dm;
						}
					}
					if (j<M[i].size()-1) {  //Binding of SR-B1
						dm=h*pss1*k1*M[i][j];
						dM_k1[i][j]=dM_k1[i][j]-dm;
						dM_k1[i][j+1]=dM_k1[i][j+1]+dm;
					}
					
					//Death from all points on membrane
					dm=h*d*M[i][j];
					D_k1=D_k1+dm;
					dM_k1[i][j]=dM_k1[i][j]-dm;
					
					//Entry process given sufficient CD81
					if (i==dim) {
						dm=h*de*M[i][j];
						E_k1=E_k1+dm; //In
						dM_k1[i][j]=dM_k1[i][j]-dm; //Out
					}
				}
			}
		}
		
		//Model is coded to here...
		
		//Step two
		D_k2=0;
		E_k2=0;
		vector< vector<double> > temp;
		temp=dM_k1;
		HalfTempMat(temp);
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
				if (domask==0||mask[i][j]>0) {
					double dm=0;
					double k1=0;
					if (Ne-i-j>=0) {
						k1=rates[Ne-i-j];
					}
					double k2=rates[j];;
					if (i<M.size()-1) { //SR-B1 independent binding
						if (j<Ne) {
							dm=h*((pec1*k1*(M[i][j]+temp[i][j]))-(c3*(i+1)*(M[i+1][j]+temp[i+1][j])));
							dM_k2[i][j]=dM_k2[i][j]-dm;
							dM_k2[i+1][j]=dM_k2[i+1][j]+dm;
						}
						if (j>0) { //SR-B1 dependent binding
							dm=h*pec2*k2*(M[i][j]+temp[i][j]);
							dM_k2[i][j]=dM_k2[i][j]-dm;
							dM_k2[i+1][j-1]=dM_k2[i+1][j-1]+dm;
						}
					}
					if (j<M[i].size()-1) {  //Binding of SR-B1
						dm=h*pss1*k1*(M[i][j]+temp[i][j]);
						dM_k2[i][j]=dM_k2[i][j]-dm;
						dM_k2[i][j+1]=dM_k2[i][j+1]+dm;
					}
					
					//Death from all points on membrane
					dm=h*d*(M[i][j]+temp[i][j]);
					D_k2=D_k2+dm;
					dM_k2[i][j]=dM_k2[i][j]-dm;
					
					//Entry process given sufficient CD81
					if (i==dim) {
						dm=h*de*(M[i][j]+temp[i][j]);
						E_k2=E_k2+dm; //In
						dM_k2[i][j]=dM_k2[i][j]-dm; //Out
					}
				}
			}
		}
		
		//Step three
		D_k3=0;
		E_k3=0;
		temp=dM_k2;
		HalfTempMat(temp);
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
				if (domask==0||mask[i][j]>0) {
					double dm=0;
					double k1=0;
					if (Ne-i-j>=0) {
						k1=rates[Ne-i-j];
					}
					double k2=rates[j];;
					if (i<M.size()-1) { //SR-B1 independent binding
						if (j<Ne) {
							dm=h*((pec1*k1*(M[i][j]+temp[i][j]))-(c3*(i+1)*(M[i+1][j]+temp[i+1][j])));
							dM_k3[i][j]=dM_k3[i][j]-dm;
							dM_k3[i+1][j]=dM_k3[i+1][j]+dm;
						}
						if (j>0) { //SR-B1 dependent binding
							dm=h*pec2*k2*(M[i][j]+temp[i][j]);
							dM_k3[i][j]=dM_k3[i][j]-dm;
							dM_k3[i+1][j-1]=dM_k3[i+1][j-1]+dm;
						}
					}
					if (j<M[i].size()-1) {  //Binding of SR-B1
						dm=h*pss1*k1*(M[i][j]+temp[i][j]);
						dM_k3[i][j]=dM_k3[i][j]-dm;
						dM_k3[i][j+1]=dM_k3[i][j+1]+dm;
					}
					
					//Death from all points on membrane
					dm=h*d*(M[i][j]+temp[i][j]);
					D_k3=D_k3+dm;
					dM_k3[i][j]=dM_k3[i][j]-dm;
					
					//Entry process given sufficient CD81
					if (i==dim) {
						dm=h*de*(M[i][j]+temp[i][j]);
						E_k3=E_k3+dm; //In
						dM_k3[i][j]=dM_k3[i][j]-dm; //Out
					}
				}
			}
		}
		
		//Step four
		D_k4=0;
		E_k4=0;
		temp=dM_k3;
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
				if (domask==0||mask[i][j]>0) {
					double dm=0;
					double k1=0;
					if (Ne-i-j>=0) {
						k1=rates[Ne-i-j];
					}
					double k2=rates[j];;
					if (i<M.size()-1) { //SR-B1 independent binding
						if (j<Ne) {
							dm=h*((pec1*k1*(M[i][j]+temp[i][j]))-(c3*(i+1)*(M[i+1][j]+temp[i+1][j])));
							dM_k4[i][j]=dM_k4[i][j]-dm;
							dM_k4[i+1][j]=dM_k4[i+1][j]+dm;
						}
						if (j>0) { //SR-B1 dependent binding
							dm=h*pec2*k2*(M[i][j]+temp[i][j]);
							dM_k4[i][j]=dM_k4[i][j]-dm;
							dM_k4[i+1][j-1]=dM_k4[i+1][j-1]+dm;
						}
					}
					if (j<M[i].size()-1) {  //Binding of SR-B1
						dm=h*pss1*k1*(M[i][j]+temp[i][j]);
						dM_k4[i][j]=dM_k4[i][j]-dm;
						dM_k4[i][j+1]=dM_k4[i][j+1]+dm;
					}
					
					//Death from all points on membrane
					dm=h*d*(M[i][j]+temp[i][j]);
					D_k4=D_k4+dm;
					dM_k4[i][j]=dM_k4[i][j]-dm;
					
					//Entry process given sufficient CD81
					if (i==dim) {
						dm=h*de*(M[i][j]+temp[i][j]);
						E_k4=E_k4+dm; //In
						dM_k4[i][j]=dM_k4[i][j]-dm; //Out
					}
				}
			}
		}
		
		
		//Final step
		D=D+(1./6.)*(D_k1+(2*D_k2)+(2*D_k3)+D_k4);
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
				M[i][j]=M[i][j]+(1./6.)*(dM_k1[i][j]+(2*dM_k2[i][j])+(2*dM_k3[i][j])+dM_k4[i][j]);
			}
		}
		//Calculate entry - all viruses that have the required number of CD81 receptors.
		E=E+(1./6.)*(E_k1+(2*E_k2)+(2*E_k3)+E_k4);
		
		//Calculate proportion of viruses that have gained entry
		
		//Proportion of viruses attached to the membrane
		B = 0;
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
				B=B+M[i][j];
			}
		}
		//cout << "B " << B << "\n";
		
		int err=0;
		//Correct for machine precision errors
		double T=D+B+E;
		//cout << "T " << T << "\n";
		D=D/T;
		E=E/T;
		B=B/T;
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
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
		
		if (verb==1) {
			cout << D << " " << E << " " << B << " " << D+E << " ";
			for (int i=0;i<M.size();i++) {
				for (int j=0;j<M[i].size();j++) {
					cout << M[i][j] << " ";
				}
			}
			cout << "\n";
		}
		//cout << "D " << D << "\n";
		//cout << "E " << E << "\n";
		
		chk=D+E;
		
		int change=0;
		//Clear up the matrix M.  Remove things less than 1e-20.
		CleanupMatrix (change,cut,Ne,pec1,pec2,pss1,cs,M);
		
		//Update the value of rs - exclude columns once they are cleaned up
		UpdateColumnsMat (change,rs,cs,M);
	}
	//cout << "E = " << E << "\n";
	return E;
	
}


double BasicModelRK22 (int dim, int verb, int quickprop, double h, double pe, double ps, double d, int Ne, double c1, double c2, double s1, double s2, double de) {
	
	//Set up initial variables
	double D=0;  //Dead viruses - get here from every other state at Poisson rate
	double E=0;	 //Viruses which have gained entry
	vector< vector<double> > M;  //From M[0][0] up.  Question - how far to go with SRB1?
	vector< vector<double> > dM;
	vector< vector<double> > zero;
	vector< vector<int> > mask;
	MatrixSetup6(dim,Ne,M,dM,zero,mask);
	vector<double> rates;
	for (int i=0;i<=Ne;i++) {
		double r=(i+0.)/(Ne+0.);
		rates.push_back(r);
	}
	
	//Row start
	int rs=0;
	
	//Column start vector
	vector<int> cs;
	for (int i=0;i<M.size();i++) {
		cs.push_back(0);
	}
	
	M[0][0]=1; //Population starts at M[0][0]; attached to cell, but not bound to anything
	
	double cut=1e-20;
	//MakeMask(cut,M,mask);
	
	//Set up Runge-Kutta variables
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
	double pec1=pe*c1;
	double pec2=pe*c2;
	double pss1=ps*s1;
	double chk=0;
	
	int domask=1;
	
	
	while (chk<0.999) {
		
		
		//Proportion of viruses attached to the membrane
		double B = 0;
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
				B=B+M[i][j];
			}
		}
		
		if (rs==M.size()-1) {
			//cout << "Only one column left\n";
			//At this point, don't need to model the complete matrix any more.  Only outputs are to E and D
			if (quickprop==1) {
				FinalCalcMat (rs,cs,d,de,D,E,M);
				return E;
			}
		}
		
		//Mask elements - don't need to include all of them in the calcuation
		if (domask==1) {
			MakeMask(cut,M,mask);
		}
		//PrintMatrix(mask);
		
		dM=zero;
		dM_k1=dM;
		dM_k2=dM;
		dM_k3=dM;
		dM_k4=dM;
		
		//Note: Seems to be some problem in the differential equations - T is getting too high.
		
		//Step one
		D_k1=0;
		E_k1=0;
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
				if (domask==0||mask[i][j]>0) {
					double dm=0;
					double k1=0;
					if (Ne-i-j>=0) {
						k1=rates[Ne-i-j];
					}
					double k2=rates[j];;
					if (i<M.size()-1) { //SR-B1 independent binding
						if (j<Ne) {
							dm=h*pec1*k1*M[i][j];
							dM_k1[i][j]=dM_k1[i][j]-dm;
							dM_k1[i+1][j]=dM_k1[i+1][j]+dm;
						}
						if (j>0) { //SR-B1 dependent binding
							dm=h*pec2*k2*M[i][j];
							dM_k1[i][j]=dM_k1[i][j]-dm;
							dM_k1[i+1][j-1]=dM_k1[i+1][j-1]+dm;
						}
					}
					if (j<M[i].size()-1) {  //Binding of SR-B1
						dm=h*((pss1*k1*M[i][j])-(s2*(j+1)*M[i][j+1]));
						dM_k1[i][j]=dM_k1[i][j]-dm;
						dM_k1[i][j+1]=dM_k1[i][j+1]+dm;
					}
					
					//Death from all points on membrane
					dm=h*d*M[i][j];
					D_k1=D_k1+dm;
					dM_k1[i][j]=dM_k1[i][j]-dm;
					
					//Entry process given sufficient CD81
					if (i==dim) {
						dm=h*de*M[i][j];
						E_k1=E_k1+dm; //In
						dM_k1[i][j]=dM_k1[i][j]-dm; //Out
					}
				}
			}
		}
		
		//Model is coded to here...
		
		//Step two
		D_k2=0;
		E_k2=0;
		vector< vector<double> > temp;
		temp=dM_k1;
		HalfTempMat(temp);
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
				if (domask==0||mask[i][j]>0) {
					double dm=0;
					double k1=0;
					if (Ne-i-j>=0) {
						k1=rates[Ne-i-j];
					}
					double k2=rates[j];;
					if (i<M.size()-1) { //SR-B1 independent binding
						if (j<Ne) {
							dm=h*pec1*k1*(M[i][j]+temp[i][j]);
							dM_k2[i][j]=dM_k2[i][j]-dm;
							dM_k2[i+1][j]=dM_k2[i+1][j]+dm;
						}
						if (j>0) { //SR-B1 dependent binding
							dm=h*pec2*k2*(M[i][j]+temp[i][j]);
							dM_k2[i][j]=dM_k2[i][j]-dm;
							dM_k2[i+1][j-1]=dM_k2[i+1][j-1]+dm;
						}
					}
					if (j<M[i].size()-1) {  //Binding of SR-B1
						dm=h*((pss1*k1*(M[i][j]+temp[i][j]))-(s2*(j+1)*(M[i][j+1]+temp[i][j+1])));
						dM_k2[i][j]=dM_k2[i][j]-dm;
						dM_k2[i][j+1]=dM_k2[i][j+1]+dm;
					}
					
					//Death from all points on membrane
					dm=h*d*(M[i][j]+temp[i][j]);
					D_k2=D_k2+dm;
					dM_k2[i][j]=dM_k2[i][j]-dm;
					
					//Entry process given sufficient CD81
					if (i==dim) {
						dm=h*de*(M[i][j]+temp[i][j]);
						E_k2=E_k2+dm; //In
						dM_k2[i][j]=dM_k2[i][j]-dm; //Out
					}
				}
			}
		}
		
		//Step three
		D_k3=0;
		E_k3=0;
		temp=dM_k2;
		HalfTempMat(temp);
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
				if (domask==0||mask[i][j]>0) {
					double dm=0;
					double k1=0;
					if (Ne-i-j>=0) {
						k1=rates[Ne-i-j];
					}
					double k2=rates[j];;
					if (i<M.size()-1) { //SR-B1 independent binding
						if (j<Ne) {
							dm=h*pec1*k1*(M[i][j]+temp[i][j]);
							dM_k3[i][j]=dM_k3[i][j]-dm;
							dM_k3[i+1][j]=dM_k3[i+1][j]+dm;
						}
						if (j>0) { //SR-B1 dependent binding
							dm=h*pec2*k2*(M[i][j]+temp[i][j]);
							dM_k3[i][j]=dM_k3[i][j]-dm;
							dM_k3[i+1][j-1]=dM_k3[i+1][j-1]+dm;
						}
					}
					if (j<M[i].size()-1) {  //Binding of SR-B1
						dm=h*((pss1*k1*(M[i][j]+temp[i][j]))-(s2*(j+1)*(M[i][j+1]+temp[i][j+1])));
						dM_k3[i][j]=dM_k3[i][j]-dm;
						dM_k3[i][j+1]=dM_k3[i][j+1]+dm;
					}
					
					//Death from all points on membrane
					dm=h*d*(M[i][j]+temp[i][j]);
					D_k3=D_k3+dm;
					dM_k3[i][j]=dM_k3[i][j]-dm;
					
					//Entry process given sufficient CD81
					if (i==dim) {
						dm=h*de*(M[i][j]+temp[i][j]);
						E_k3=E_k3+dm; //In
						dM_k3[i][j]=dM_k3[i][j]-dm; //Out
					}
				}
			}
		}
		
		//Step four
		D_k4=0;
		E_k4=0;
		temp=dM_k3;
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
				if (domask==0||mask[i][j]>0) {
					double dm=0;
					double k1=0;
					if (Ne-i-j>=0) {
						k1=rates[Ne-i-j];
					}
					double k2=rates[j];;
					if (i<M.size()-1) { //SR-B1 independent binding
						if (j<Ne) {
							dm=h*pec1*k1*(M[i][j]+temp[i][j]);
							dM_k4[i][j]=dM_k4[i][j]-dm;
							dM_k4[i+1][j]=dM_k4[i+1][j]+dm;
						}
						if (j>0) { //SR-B1 dependent binding
							dm=h*pec2*k2*(M[i][j]+temp[i][j]);
							dM_k4[i][j]=dM_k4[i][j]-dm;
							dM_k4[i+1][j-1]=dM_k4[i+1][j-1]+dm;
						}
					}
					if (j<M[i].size()-1) {  //Binding of SR-B1
						dm=h*((pss1*k1*(M[i][j]+temp[i][j]))-(s2*(j+1)*(M[i][j+1]+temp[i][j+1])));
						dM_k4[i][j]=dM_k4[i][j]-dm;
						dM_k4[i][j+1]=dM_k4[i][j+1]+dm;
					}
					
					//Death from all points on membrane
					dm=h*d*(M[i][j]+temp[i][j]);
					D_k4=D_k4+dm;
					dM_k4[i][j]=dM_k4[i][j]-dm;
					
					//Entry process given sufficient CD81
					if (i==dim) {
						dm=h*de*(M[i][j]+temp[i][j]);
						E_k4=E_k4+dm; //In
						dM_k4[i][j]=dM_k4[i][j]-dm; //Out
					}
				}
			}
		}
		
		
		//Final step
		D=D+(1./6.)*(D_k1+(2*D_k2)+(2*D_k3)+D_k4);
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
				M[i][j]=M[i][j]+(1./6.)*(dM_k1[i][j]+(2*dM_k2[i][j])+(2*dM_k3[i][j])+dM_k4[i][j]);
			}
		}
		//Calculate entry - all viruses that have the required number of CD81 receptors.
		E=E+(1./6.)*(E_k1+(2*E_k2)+(2*E_k3)+E_k4);
		
		//Calculate proportion of viruses that have gained entry
		
		//Proportion of viruses attached to the membrane
		B = 0;
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
				B=B+M[i][j];
			}
		}
		//cout << "B " << B << "\n";
		
		int err=0;
		//Correct for machine precision errors
		double T=D+B+E;
		//cout << "T " << T << "\n";
		D=D/T;
		E=E/T;
		B=B/T;
		for (int i=rs;i<M.size();i++) {
			for (int j=cs[i];j<M[i].size();j++) {
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
		
		if (verb==1) {
			cout << D << " " << E << " " << B << " " << D+E << " ";
			for (int i=0;i<M.size();i++) {
				for (int j=0;j<M[i].size();j++) {
					cout << M[i][j] << " ";
				}
			}
			cout << "\n";
		}
		//cout << "D " << D << "\n";
		//cout << "E " << E << "\n";
		
		chk=D+E;
		
		int change=0;
		//Clear up the matrix M.  Remove things less than 1e-20.
		CleanupMatrix (change,cut,Ne,pec1,pec2,pss1,cs,M);
		
		//Update the value of rs - exclude columns once they are cleaned up
		UpdateColumnsMat (change,rs,cs,M);
	}
	//cout << "E = " << E << "\n";
	return E;
	
}

//Need to change this to account for the model....
void CleanupMatrix (int& change, double cut, double Ne, double pec1, double pec2, double pss1, vector<int>& cs, vector< vector<double> >& M) {
    //cout << "Cleanup matrix\n";
    for (int i=0;i<M.size();i++) {
        if (cs[i]==M[i].size()-1||M[i][cs[i]+1]>M[i][cs[i]]) {
            //Test that the column has been activated i.e. flow has occured through SRB1
            for (int j=cs[i];j<M[i].size();j++) {
                if (M[i][j]>0) {
                    int stop=1;
                    if (M[i][j]<cut) {
                        //If an element of the matrix is less than the cutoff and will not receive input from other
                        //elements in future, set this to zero.  Skip this element in future calculations.
                        if (i==0||cs[i-1]>j) { //Here check that we are not receiving input from M[i-1][j+1]
                            stop=0;
                            //cout << "Remove element " << i << " " << j << "\n";
							
                            //Deterministic removal of M[i][j] into adjacent cells
                            double k1=(Ne-i+0.)/(Ne+0.);
                            double k2=(Ne-j+0.)/(Ne+0.);
                            double k3=(j+0.)/(Ne+0.);
							
                            double out_cd81_1=pec1*k1*k2*M[i][j];
                            double out_cd81_2=pec2*k1*k3*M[i][j];
                            if (i==M.size()-1) {
                                out_cd81_1=0;
                                out_cd81_2=0;
                            }
                            if (j==0) {
                                out_cd81_2=0;
                            }
                            double out_srb1=pss1*k2*M[i][j];
                            double out_tot=out_cd81_1+out_cd81_2+out_srb1;
                            double prop1=out_cd81_1/out_tot;
                            double prop2=out_cd81_2/out_tot;
                            double prop3=out_srb1/out_tot;
							
                            //Do the clearance
                            if (i<M.size()-1) {
                                M[i+1][j]=M[i+1][j]+(prop1*M[i][j]);
                                if (j>0) {
                                    M[i+1][j-1]=M[i+1][j-1]+(prop2*M[i][j]);
                                }
                            }
                            M[i][j+1]=M[i][j+1]+(prop3*M[i][j]);
                            M[i][j]=0;
							
                            //Update column start position
                            cs[i]=j;
                            change=1;
                        }
                    }
                    if (stop==1) {
                        break;
                    }
                }
            }
        }
    }
}

void UpdateColumnsMat (int& change, int& rs, vector<int>& cs, vector< vector<double> >& M) {
    //Here rs indicates the first non-zero column.  Update this when the final element of the first column is removed
    if (change==1) {
        //cout << "Update columns\n";
        for (int i=rs;i<cs.size();i++) {
            //cout << i << " " << cs[i] << "\n";
            if (cs[i]==M[i].size()-i-1) {
                //cout << "Add to rs\n";
                rs++;
            }
        }
    }
}

void MakeMask (double cut, vector< vector<double> >& M, vector< vector<int> >& mask) {
    for (int i=0;i<M.size();i++) {
        for (int j=0;j<M[i].size();j++){
            mask[i][j]=0;
            //Elements of mask based on movement from previous states
            if (j>0) {
                if (mask[i][j-1]>0) {
                    mask[i][j]=mask[i][j-1]-1;
                }
            }
            if (i>0) {
                if (mask[i-1][j]>mask[i][j]+1) {
                    mask[i][j]=mask[i-1][j]-1;
                }
                if (j<mask[i].size()-1) {
                    if (mask[i-1][j+1]>mask[i][j]+1) {
                        mask[i][j]=mask[i-1][j+1]-1;
                    }
                }
            }
            //Elements of mask based on frequency threshold
            if (M[i][j]>cut) {
                mask[i][j]=5;
            }
        }
    }
}


void UpdateColumns (int& rs, double cut, vector<double>& M) {
    if (M[rs]<1e-20) {
        M[rs+1]=M[rs+1]+M[rs];
        rs++;
    }
}

void MatrixSetup (int dim, vector<double>& M, vector<double>& dM, vector<double>& zero) {
	for (int i=0;i<=dim;i++) {
        M.push_back(0);
	}
	dM=M;
	zero=M;
}

void MatrixSetup6 (int dim, int Ne, vector< vector<double> >& M, vector< vector<double> >& dM, vector< vector<double> >& zero, vector< vector<int> >& mask) {
    for (int i=0;i<=dim;i++) {
        //cout << "i " << i << " " << dim << "\n";
        vector<double> mm;
        for (int j=0;j<=Ne;j++) {
            mm.push_back(0);
        }
        M.push_back(mm);
    }
    dM=M;
    zero=M;
    for (int i=0;i<=dim;i++) {
        vector<int> mm;
        for (int j=0;j<=Ne;j++) {
            mm.push_back(0);
        }
        mask.push_back(mm);
    }
}



void HalfTemp(vector<double>& temp) {
	for (int i=0;i<temp.size();i++) {
			temp[i]=temp[i]/2;
	}
}

void HalfTempMat(vector< vector<double> >& temp) {
    for (int i=0;i<temp.size();i++) {
        for (int j=0;j<temp[i].size();j++) {
            temp[i][j]=temp[i][j]/2;
        }
    }
}


void FinalCalc (int& rs, double& d, double& de, double& D, double& E, vector<double>& M) {
	//Calculation for when only the last element of the vector is occupied.  No more binding to CD81 implies no need to track the full system
	double Mtot=M[rs];

    //cout << "E " << E << " D " << D << " M " << Mtot << "\n";
	double e_prop=de/(d+de);
	double d_prop=d/(d+de);
	//cout << "Props " << e_prop << " " << d_prop << "\n";
	
	double rem=(0.999-(E+D))/Mtot;
	//cout << "Remainder " << rem << "\n";
	E=E+(rem*e_prop*Mtot);
	D=D+(rem*d_prop*Mtot);
}

void FinalCalcMat (int& rs, vector<int>& cs, double& d, double& de, double& D, double& E, vector< vector<double> >& M) {
    //Calculation for when only the last column is occupied.  No more binding to CD81 implies no need to track the full system
    double Mtot=0;
    for (int i=rs;i<M.size();i++) {
        for (int j=cs[i];j<M[i].size();j++) {
            Mtot=Mtot+M[i][j];
        }
    }
    //cout << "E " << E << " D " << D << " M " << Mtot << "\n";
    double e_prop=de/(d+de);
    double d_prop=d/(d+de);
    //cout << "Props " << e_prop << " " << d_prop << "\n";
    
    double rem=(0.999-(E+D))/Mtot;
    //cout << "Remainder " << rem << "\n";
    E=E+(rem*e_prop*Mtot);
    D=D+(rem*d_prop*Mtot);
}



void FindLikelihoodRK (run_params p, double cons, double t, double k, double& h, double theta, double& lL, int maxv, int verb, int verb_model, int dim, vector<double> params, vector<double> cd81_val, vector<double> srb_val, vector< vector<dat> >& all_dat, vector<double> fact_store) {
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

	if (verb_model==1) { //Print progress only for the default case
        double e=-10;
		if (p.model==0) {
			e=BasicModelRK0(dim,verb_model,1,h,1,1,d,c1,Ne,de);
        } else if (p.model==6) {
            e=BasicModelRK6(dim,verb_model,1,h,1,1,d,Ne,c1,c2,s1,de);
		} else if (p.model==21) {
			e=BasicModelRK21(dim,verb_model,1,h,1,1,d,Ne,c1,c2,c3,s1,de);
		} else if (p.model==22) {
			e=BasicModelRK22(dim,verb_model,1,h,1,1,d,Ne,c1,c2,s1,s2,de);
		}
	}
    double h_store=h;
	
	//Calculate values across range of CD81, keeping SRB1 fixed
	vector<double> evals_cd81;
	for (int i=0;i<cd81_val.size();i++) {
		double e=-10;
		h=h*2;
		while (e==-10) {
			e=-10;
			if (p.model==0) {
				e=BasicModelRK0(dim,0,1,h,cd81_val[i],1,d,c1,Ne,de);
            } else if (p.model==6) {
                e=BasicModelRK6(dim,verb_model,1,h,cd81_val[i],1,d,Ne,c1,c2,s1,de);
			} else if (p.model==21) {
				e=BasicModelRK21(dim,verb_model,1,h,cd81_val[i],1,d,Ne,c1,c2,c3,s1,de);
			} else if (p.model==22) {
				e=BasicModelRK22(dim,verb_model,1,h,cd81_val[i],1,d,Ne,c1,c2,s1,s2,de);
			}
            h=h/2;
		}
        //cout << i << " " << cd81_val[i] << " " << e << "\n";

		evals_cd81.push_back(e); //Gives proportion of viruses that get into a cell
	}
	vector<double> evals_srb1;

    //Calculate values across range of SRB1, keeping CD81 fixed
    for (int i=0;i<srb_val.size();i++) {
		double e=-10;
		h=h*2;
		while (e==-10) {
			e=-10;
			if (p.model==0) {
				e=BasicModelRK0(dim,0,1,h,1,srb_val[i],d,c1,Ne,de);
            } else if (p.model==6) {
                e=BasicModelRK6(dim,0,1,h,1,srb_val[i],d,Ne,c1,c2,s1,de);
			} else if (p.model==21) {
				e=BasicModelRK21(dim,0,1,h,1,srb_val[i],d,Ne,c1,c2,c3,s1,de);
			} else if (p.model==22) {
				e=BasicModelRK22(dim,0,1,h,1,srb_val[i],d,Ne,c1,c2,s1,s2,de);
			}
            h=h/2;
		}
		evals_srb1.push_back(e);
        //cout << i << " " << e << "\n";
	}

    h=h_store;
    
    //Find probabilities of entry given k input viruses
	vector< vector<double> > pvals_cd81;
	vector< vector<double> > pvals_srb1;
	ConstructPVals(maxv,evals_cd81,evals_srb1,pvals_cd81,pvals_srb1);
	lL=CalcLike(verb,theta,pvals_cd81,pvals_srb1,all_dat,fact_store);
    if (p.usegamma==1) {
		rat=s1/c1;
        double gl=CalcGamma(cons,t,k,rat);
        //cout << "Ratio " << rat << " " << gl << "\n";
        lL=lL+gl;
    }
}

void ConstructPVals (int maxv, vector<double> evals_cd81, vector<double> evals_srb1, vector< vector<double> >& pvals_cd81, vector< vector<double> >& pvals_srb1) {
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
	//		cout << "p " << p << "\n";
			if (p==0) {
				p=1e-10;
			}
            p=p*all_dat[i][j].count;
			double val=CalcDoublePoisson(p,theta,all_dat[i][j].foci,fact_store);
			lL=lL+val;
		}
	}
	//cout << "Log likelihood " << lL << "\n";
	if (verb==2) {
		OutputData(all_dat,cd81_x,srb_x);
	}
	if (verb>=1) {
		OutputFreqs(cd81_indices,srb_indices,cd81_x,srb_x);
	}
	//cout << "L " << lL << "\n";
	return lL;
}

