using namespace std;
#include "basicmodel.h"
#include "io.h"

void GetOptions (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.dim=1;
	p.model=1;
	p.max_its=1000000;
	p.N=100;
	p.r=0.01;
	p.c1=0.01;
	p.c2=0.02;
    p.c3=0.01;
	p.de=0.01;
	p.s1=0.01;
    p.s2=0.01;
    p.Nc=100;
    p.Ns=100;
    p.Ne=100;
	p.gp1=0;
	p.gp2=0;
    p.usegamma=0;
    p.precision=0;
    p.smallchange=0;
    p.threshold=1e-5;
	p.grid=0;
    p.nfix=0;
	int x=1;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--dim")==0) {
			x++;
			p.dim=atoi(argv[x]);
		} else if (p_switch.compare("--model")==0) {
			x++;
			p.model=atoi(argv[x]);
		} else if (p_switch.compare("--its")==0) {
			x++;
			p.max_its=atoi(argv[x]);
		} else if (p_switch.compare("--s1")==0) {
			x++;
			p.s1=atof(argv[x]);
        } else if (p_switch.compare("--s2")==0) {
            x++;
            p.s2=atof(argv[x]);
		} else if (p_switch.compare("--r")==0) {
			x++;
			p.r=atof(argv[x]);
		} else if (p_switch.compare("--N")==0) {
			x++;
			p.N=atof(argv[x]);
        } else if (p_switch.compare("--Nc")==0) {
            x++;
            p.Nc=atof(argv[x]);
        } else if (p_switch.compare("--Ns")==0) {
            x++;
            p.Ns=atof(argv[x]);
        } else if (p_switch.compare("--Ne")==0) {
            x++;
            p.Ne=atoi(argv[x]);
		} else if (p_switch.compare("--c1")==0) {
			x++;
			p.c1=atof(argv[x]);
		} else if (p_switch.compare("--c2")==0) {
			x++;
			p.c2=atof(argv[x]);
        } else if (p_switch.compare("--c3")==0) {
            x++;
            p.c3=atof(argv[x]);
		} else if (p_switch.compare("--de")==0) {
			x++;
			p.de=atof(argv[x]);
        } else if (p_switch.compare("--precision")==0) {
            x++;
            p.precision=atoi(argv[x]);
        } else if (p_switch.compare("--threshold")==0) {
            x++;
            p.threshold=atof(argv[x]);
        } else if (p_switch.compare("--nfix")==0) {
            x++;
            p.nfix=atoi(argv[x]);
        } else if (p_switch.compare("--usegamma")==0) {
            x++;
            p.usegamma=atoi(argv[x]);
      } else if (p_switch.compare("--small")==0) {
            x++;
            p.smallchange=atoi(argv[x]);
		} else if (p_switch.compare("--makegrid")==0) {
			x++;
			p.grid=atoi(argv[x]);
		} else if (p_switch.compare("--gp1")==0) {
			x++;
			p.gp1=atof(argv[x]);
		} else if (p_switch.compare("--gp2")==0) {
			x++;
			p.gp2=atof(argv[x]);
        } else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void GetAllDat (vector< vector<dat> >& all_dat) {
	vector<dat> d1;
	GetData("../Data/Set3_CD81_final.dat",1,d1);
	all_dat.push_back(d1);
	d1.clear();
	GetData("../Data/Set3_SRB1_final.dat",0,d1);
	all_dat.push_back(d1);
	d1.clear();
	GetData("../Data/Set5_CD81_final.dat",1,d1);
	all_dat.push_back(d1);
	d1.clear();
	GetData("../Data/Set5_SRB1_final.dat",0,d1);
	all_dat.push_back(d1);
	d1.clear();
	GetData("../Data/Set6_CD81_final.dat",1,d1);
	all_dat.push_back(d1);
	d1.clear();
	GetData("../Data/Set6_SRB1_final.dat",0,d1);
	all_dat.push_back(d1);
	d1.clear();
	GetData("../Data/Set7_CD81_final.dat",1,d1);
	all_dat.push_back(d1);
	d1.clear();
	GetData("../Data/Set7_SRB1_final.dat",0,d1);
	all_dat.push_back(d1);
	d1.clear();
}

void GetData (string name, int cd, vector<dat>& dat1) {
	ifstream dat_file1;
	dat_file1.open(name.c_str());
	double x;
	int n;
	char c;
	for (int i=0;i<1000;i++) {
		dat d;
		d.cd81=cd;
		if (!(dat_file1 >> c)) break;
		if (!(dat_file1 >> n)) break;
		if (!(dat_file1 >> x)) break;
		d.a=x;
		if (!(dat_file1 >> x)) break;
		d.sd=x;
		if (!(dat_file1 >> n)) break;
		if (!(dat_file1 >> x)) break;
		d.growth=x;
		if (!(dat_file1 >> x)) break;
		d.count=x;
		if (!(dat_file1 >> n)) break;
		d.foci=n;
		dat1.push_back(d);
	}
}

void ConstructOutputVectors(vector< vector<dat> >& all_dat, vector<double>& cd81_indices, vector<double>& srb_indices, vector<double>& cd81_x, vector<double>& srb_x) {
	int max_c=0;
	int max_s=0;
	for (int i=0;i<all_dat.size();i++) {
		for (int j=0;j<all_dat[i].size();j++) {
			if (all_dat[i][j].cd81==1) {
				if (all_dat[i][j].index>max_c) {
					max_c=all_dat[i][j].index;
				}
			} else {
				if (all_dat[i][j].index>max_s) {
					max_s=all_dat[i][j].index;
				}
			}
		}
	}
	for (int i=0;i<=max_c;i++) {
		cd81_indices.push_back(-1);
		cd81_x.push_back(-1);
	}
	for (int i=0;i<=max_s;i++) {
		srb_indices.push_back(-1);
		srb_x.push_back(-1);
	}
}

void ModifyOutputVectors (int i, int j, vector< vector<dat> >& all_dat, double p, vector<double>& cd81_indices, vector<double>& srb_indices, vector<double>& cd81_x, vector<double>& srb_x) {
	if (all_dat[i][j].cd81==1) {
		if (cd81_indices[all_dat[i][j].index]<0) {
			cd81_indices[all_dat[i][j].index]=p;
			cd81_x[all_dat[i][j].index]=all_dat[i][j].a;
		}
	} else {
		if (srb_indices[all_dat[i][j].index]<0) {
			srb_indices[all_dat[i][j].index]=p;
			srb_x[all_dat[i][j].index]=all_dat[i][j].a;
		}
	}
}

void OutputData (vector< vector<dat> >& all_dat, vector<double>& cd81_x, vector<double>& srb_x) {
	ofstream data_file;
	data_file.open("Data.out");
	data_file << "CD81\n";
	for (int i=0;i<cd81_x.size();i++) {
		data_file << cd81_x[i] << " ";
		for (int j=0;j<all_dat.size();j++) {
			for (int k=0;k<all_dat[j].size();k++) {
				if (all_dat[j][k].cd81==1) {
					if (all_dat[j][k].index==i) {
						data_file << (all_dat[j][k].foci+0.)/all_dat[j][k].count << " ";
					}
				}
			}
		}
		data_file << "\n";
	}
	data_file << "SRB1\n";
	for (int i=0;i<srb_x.size();i++) {
		data_file << srb_x[i] << " ";
		for (int j=0;j<all_dat.size();j++) {
			for (int k=0;k<all_dat[j].size();k++) {
				if (all_dat[j][k].cd81==0) {
					if (all_dat[j][k].index==i) {
						data_file << (all_dat[j][k].foci+0.)/all_dat[j][k].count << " ";
					}
				}
			}
		}
		data_file << "\n";
	}
}

void OutputFreqs (vector<double>& cd81_indices, vector<double>& srb_indices, vector<double>& cd81_x, vector<double>& srb_x) {
	ofstream out_file;
	out_file.open("Freqs.out");
	out_file << "CD81\n";
	for (int i=0;i<cd81_indices.size();i++) {
		out_file << cd81_x[i] << " " << cd81_indices[i] << "\n";
	}
	out_file << "SRB1\n";
	for (int i=0;i<srb_indices.size();i++) {
		out_file << srb_x[i] << " " << srb_indices[i] << "\n";
	}
	out_file.close();
}

void PrintParameters (double lL, vector<double> params, ofstream& rec_file) {
    cout << "Better lL = " << lL << " ";
    rec_file << "Better lL = " << lL << " ";
    for (int i=1;i<params.size();i++) {
        cout << params[i] << " ";
        rec_file << params[i] << " ";
    }
    cout << "\n";
    rec_file << "\n";
}
