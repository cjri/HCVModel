using namespace std;
#include "basicmodel.h"
#include "io.h"

void GetOptions (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.dim=1;
	p.model=1;
	p.rk=0;
	p.over=0;
	p.max_its=100000;
	p.e0=0.01;
	p.s0=0.01;
	p.p1=0.01;
	p.p2=0.02;
	p.de=0.01;
	p.gp1=0;
	p.gp2=0;
	p.lambda=0;
    p.precision=0;
    p.smallchange=0;
	p.show_all=0;
	p.grid=0;
	int x=1;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--dim")==0) {
			x++;
			p.dim=atoi(argv[x]);
		} else if (p_switch.compare("--model")==0) {
			x++;
			p.model=atoi(argv[x]);
		} else if (p_switch.compare("--rk")==0) {
			x++;
			p.rk=atoi(argv[x]);
		} else if (p_switch.compare("--over")==0) {
			x++;
			p.over=atoi(argv[x]);
		} else if (p_switch.compare("--its")==0) {
			x++;
			p.max_its=atoi(argv[x]);
		} else if (p_switch.compare("--e0")==0) {
			x++;
			p.e0=atof(argv[x]);
		} else if (p_switch.compare("--s")==0) {
			x++;
			p.s0=atof(argv[x]);
		} else if (p_switch.compare("--c1")==0) {
			x++;
			p.p1=atof(argv[x]);
		} else if (p_switch.compare("--c2")==0) {
			x++;
			p.p2=atof(argv[x]);
		} else if (p_switch.compare("--de")==0) {
			x++;
			p.de=atof(argv[x]);
        } else if (p_switch.compare("--precision")==0) {
            x++;
            p.precision=atoi(argv[x]);
        } else if (p_switch.compare("--small")==0) {
            x++;
            p.smallchange=atoi(argv[x]);
		} else if (p_switch.compare("--show")==0) {
			x++;
			p.show_all=atoi(argv[x]);
		} else if (p_switch.compare("--makegrid")==0) {
			x++;
			p.grid=atoi(argv[x]);
		} else if (p_switch.compare("--gp1")==0) {
			x++;
			p.gp1=atof(argv[x]);
		} else if (p_switch.compare("--gp2")==0) {
			x++;
			p.gp2=atof(argv[x]);
		} else if (p_switch.compare("--lambda")==0) {
			x++;
			p.lambda=atof(argv[x]);
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

void PrintParameters (int model, double lL, double s0, double p1, double p2, double de, ofstream& rec_file) {
	if (model==1) {
		cout << "Better lL = " << lL << " " << s0 << " " << p1 << " " << p2 << "\n";
		rec_file << "Better lL = " << lL << " " << s0 << " " << p1 << " " << p2 << "\n";
	}
	if (model==2) {
		cout << "Better lL = " << lL << " " << s0 << " " << p1 << " " << p2 << " " << de << "\n";
		rec_file << "Better lL = " << lL << " " << s0 << " " << p1 << " " << p2 << " " << de << "\n";
	}
}
