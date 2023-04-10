#define _USE_MATH_DEFINES

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <ctime>
#include <fstream>
#include <iostream>
#include <chrono>
#include <random>
#include <set>
#include <cstring>
#include <algorithm>

#define chain_size 50
#define trials_per_temp 100000

void init_uniform(Measurements main_measurements, Interactions main_interactions, Site *spin, int Ns, int L);
void init_random(Measurements main_measurements, Interactions main_interactions, Site *spin, int Ns, int L);
inline int index(int x, int y, int L);
void set_coordinates(Site *spin, int Ns, int L);
void set_nn(Site *spin, int Ns, int L);
int mod(int x, int m);
//void Wolff_Cluster_Sim(double beta, Site *spin, int *Ns, int *L, int pick_file);
//double Wolff_Sweep(double beta, Site *spin, int *Ns, int *L);
void randomize(Site *spin, int Ns, int L);
double rand1();
void clear_files(const string L_name);
void print(const string L_name, double E1, double E2, double M1, double M2, double M4, double beta, int Ns);
double theta(int p1, int p2);
int update_site(Measurements main_measurements, Interactions main_interactions, Properties main_properties, int k, double beta, Site *spin, int Ns, int L);
double MC_sweep(Measurements main_measurements, Interactions main_interactions, Properties main_properties, double beta, Site *spin, int Ns, int L);
void Monte_Carlo_Sim(Interactions main_interactions, Measurements main_measurements, Properties main_properties, double beta, Site *spin, int Ns, int L, const string L_name);

using namespace std;

const int N_nn1 = 6;    //Number of Nearest neighbors
const double Jnn = -1;  // nearest-neighbor ferromagnetic coupling
const int q = 6; //number of potts states
const double Dp = 0.14;

class Measurements {
    public:
    double Ed_curr;

    double clock_energy(Properties main_properties, Site *spin, int Ns, int L) {
        double sum = 0;
        for(int i=0; i<Ns; i++) {
           for(int k=0; k<N_nn1; k++) {
               sum += main_properties.V_clock[spin[i].potts][spin[i].nn1[k]->potts];
           }
     }
        return 0.5 * Jnn * sum;
    }

    double dipolar_energy(Interactions main_interactions, Site *spin, int Ns, int L) {
    
        double sum = 0;
        for(int i=0; i<Ns; i++) {
    
            for(int j=0; j<Ns; j++) {
                sum += main_interactions.Ud(spin, i, j, L)*spin[i].Sz*spin[j].Sz;
            }
        }
    
        return 0.5*Dp*sum;
    }

    inline double E_tot(Interactions main_interactions, Properties main_properties, Site *spin, int Ns, int L) {
    
        return this->clock_energy(main_properties, spin, Ns, L) + this->dipolar_energy(main_interactions, spin, Ns, L);
    }

    double potts_magnetization(Site *spin, int Ns, int L) {
    
        double xsum = 0, ysum = 0, sum = 0;

        for (int i = 0; i < Ns; i++) {
            xsum += cos(theta(spin[i].potts, 0));
            ysum += sin(theta(spin[i].potts, 0));
        }
        xsum /= (double) Ns;
        ysum /= (double) Ns;

        sum = sqrt((xsum*xsum) + (ysum*ysum));

        return sum;
    }

    double ising_magnetization(Site *spin, int Ns) {
        int sum = 0;
    
        for(int i=0; i<Ns; i++) {
        
            sum += spin[i].Sz;
        }
    
        return ((double) sum)/((double) Ns);
}
    
};

class Interactions {
    public:
    double *Vd;

    Interactions(int L) {
        double* Vd = new double[L*L];
    }

    void compute_Vd(int L, int N_max) {
        double dr0[2], dr[2];
    
        for(int x=0; x<L; x++)
            for(int y=0; y<L; y++) {
        
            if(x == 0 && y == 0) this->Vd[x*L + y] = 0;
            else {
        
                this->Vd[x*L + y] = 0;
            
                dr0[0] = x + 0.5*y;
                dr0[1] = 0.5*sqrt(3.)*y;
            
                for(int m=-N_max; m<=N_max; m++)
                for(int n=-N_max; n<=N_max; n++) {
                
                    dr[0] = dr0[0] + L*(m + 0.5*n);
                    dr[1] = dr0[1] + L*(0.5*sqrt(3.)*n);
                
                    double r2 = pow(dr[0], 2) + pow(dr[1], 2);
                
                    this->Vd[x*L + y] += 1./pow(r2, 1.5);
                
                
                }
            }
        }
    }

    inline double Ud(Site *spin, int i, int j, int L) {
        int dx = spin[j].lattice_pt[0] - spin[i].lattice_pt[0];
        int dy = spin[j].lattice_pt[1] - spin[i].lattice_pt[1];
    
        return 0.5*(this->Vd[mod(dx, L)*L + mod(dy, L)] + this->Vd[mod(-dx, L)*L + mod(-dy, L)]);
    }

    ~Interactions() {
        delete[] Vd;
    }
};

class Properties {
    public:
    int idx_m[q][q];
    double V_clock[q][q];
    double pb[q][q][q];

    Properties() {

        for (int i = 0; i < q; i++) {
            for (int j = 0; j < q; j++) {
                V_clock[i][j] = cos((i - j)*2.*M_PI/((double)q));
            }
        }

    idx_m[0][0] = 3;
    idx_m[0][1] = 2;
    idx_m[0][2] = 1;
    idx_m[0][3] = 0;
    idx_m[0][4] = 5;
    idx_m[0][5] = 4;
    
    idx_m[3][0] = 3;
    idx_m[3][1] = 2;
    idx_m[3][2] = 1;
    idx_m[3][3] = 0;
    idx_m[3][4] = 5;
    idx_m[3][5] = 4;
    
    idx_m[1][0] = 5;
    idx_m[1][1] = 4;
    idx_m[1][2] = 3;
    idx_m[1][3] = 2;
    idx_m[1][4] = 1;
    idx_m[1][5] = 0;
    
    idx_m[4][0] = 5;
    idx_m[4][1] = 4;
    idx_m[4][2] = 3;
    idx_m[4][3] = 2;
    idx_m[4][4] = 1;
    idx_m[4][5] = 0;
    
    idx_m[2][0] = 1;
    idx_m[2][1] = 0;
    idx_m[2][2] = 5;
    idx_m[2][3] = 4;
    idx_m[2][4] = 3;
    idx_m[2][5] = 2;
    
    idx_m[5][0] = 1;
    idx_m[5][1] = 0;
    idx_m[5][2] = 5;
    idx_m[5][3] = 4;
    idx_m[5][4] = 3;
    idx_m[5][5] = 2;
    }

    void set_pb(double beta) {
        for(int i=0; i<q; i++) {
            for(int j=0; j<q; j++) {
                for(int k=0; k<q; k++) {
    
            double f1 = cos((j - i)*2.*M_PI/((double) q));
            double f2 = cos((k - i)*2.*M_PI/((double) q));
            
            pb[i][j][k] = 1. - exp(2.*beta*Jnn*f1*f2);
            
                }
            }
        }
    }
};

class Site {
public:
    int idx;    // index 0, 1, 2, ... Ns-1
    int x, y;   // coordinates
    int Sz;     // Ising spins
    int potts;
    int Ns;
    int L;
    int lattice_pt[2];
    // idx of its four neighbors
    Site *nn1[N_nn1];

};

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

// a simple data structure for square-lattice spin

class Cluster {
public:
    vector<int> sites;
    double size;
};

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage: ./potts_parallel <Lattice Length>" << endl;
        return 1;
    }

    const int L = atof(argv[1]);
    const int Ns = L * L;

    Site* spin = new Site[Ns];
    Properties main_properties;
    Interactions main_interactions(L);
    Measurements main_measurements;

    const string L_size(argv[1]);
    const string L_name = "Potts" + L_size;
    
    clear_files(L_name);
    
    init_uniform(main_measurements, main_interactions, spin, Ns, L);

    set_coordinates(spin, Ns, L);
    set_nn(spin, Ns, L);
    
    randomize(spin, Ns, L);

    double del_T = 0.1;

    cout << "Thank you for choosing the Potts Model :)" << endl;
            
    for(double T = 1.5; T>0; T -= del_T) {

        cout << "Lattice Length = " << L << endl;
        cout << "T = " << T << endl;

        main_properties.set_pb(1./T);

        Monte_Carlo_Sim(main_interactions, main_measurements, main_properties, 1./T, spin, Ns, L, L_name);
	
	if (T <= 1.1 && T >= 0.65) {
	  del_T = 0.05;
	} else {
	  del_T = 0.1;
	}

    }


    delete[] spin;

    return 0;
} 
    
void init_uniform(Measurements main_measurements, Interactions main_interactions, Site *spin, int Ns, int L) {
    for (int i=0; i < Ns; i++) {
        spin[i].Sz = 1;
        spin[i].potts = 0;
    }
    main_measurements.Ed_curr = main_measurements.dipolar_energy(main_interactions, spin, Ns, L);
}

void init_random(Measurements main_measurements, Interactions main_interactions, Site *spin, int Ns, int L) {
    for (int i = 0; i < Ns; i++) {
        spin[i].potts = (int) (rand1() * 6);
        spin[i].Sz = (spin[i].potts % 2 == 0) ? 1 : -1;
    }
    main_measurements.Ed_curr = main_measurements.dipolar_energy(main_interactions, spin, Ns, L);
}

inline int index(int x, int y, int L) {
    return x*L+y;
}

void set_coordinates(Site *spin, int Ns, int L) {
    int y = 0;
    int x = 0;
    for(int i=0; i<Ns; i++) {
        spin[i].lattice_pt[0] = x;
        spin[i].lattice_pt[1] = y;
        spin[i].idx = i;
        spin[i].x = x + 0.5*y;
        spin[i].y = 0.5*sqrt(3)*y;
        x++;
        if (x%L == 0) {
            y++;
            x = 0;
        }
    }
}

void randomize(Site *spin, int Ns, int L) { //sets all site spins randomly to any of the q possible states
    for(int i=0; i<Ns; i++) {
        double r = rand1();
        spin[i].potts = (int)(r*6);
    } 
}

void set_nn(Site *spin, int Ns, int L) { //sets the 2*D nearest neighbors for each site in the system
    int k;
    int x, y;
    for(int i=0; i<Ns; i++) { 
        x = spin[i].lattice_pt[0];
		y = spin[i].lattice_pt[1];

        k = index(mod(x+1, L), y, L);
        spin[i].nn1[0] = &spin[k];
        
        k = index(x, mod(y+1, L), L);
        spin[i].nn1[1] = &spin[k];
        
        k = index(mod(x-1, L), mod(y+1, L), L);
        spin[i].nn1[2] = &spin[k];
        
        k = index(mod(x-1, L), y, L);
        spin[i].nn1[3] = &spin[k];
        
        k = index(x, mod(y-1, L), L);
        spin[i].nn1[4] = &spin[k];
        
        k = index(mod(x+1, L), mod(y-1, L), L);
        spin[i].nn1[5] = &spin[k];

    }
}

int mod(int x, int m) {
	if (x>=0 && x<m)
		return x;
	else if (x<0)
		return m-1-mod(-1-x,m);
	else
		return x%m;
}

double theta(int p1, int p2) {
    double theta;

    theta = (2*M_PI*((double)p1-(double)p2))/(double)q;
    return theta;
}

//test a new state, accept or reject the state based on parameters. 
int update_site(Measurements main_measurements, Interactions main_interactions, Properties main_properties, int k, double beta, Site *spin, int Ns, int L) {

    int p_new = (int) (rand1() * (double)q);
    int Sz_new = (p_new % 2 == 0) ? 1 : -1;

    double delPSpin = 0;
    double delISpin = 0;
    double delE = 0;

    for(int l=0; l<N_nn1; l++) {
        delPSpin += (main_properties.V_clock[p_new][spin[k].nn1[l]->potts] - main_properties.V_clock[spin[k].potts][spin[k].nn1[l]->potts]);
    }

    if (Sz_new = spin[k].Sz) {
        delISpin = 0;
    } else {
        for (int m = 0; m < Ns; m++) {
            delISpin += (Sz_new - spin[k].Sz)*spin[k].Sz*main_interactions.Ud(spin, k, m, L);
        }
    }

    //change in energy when spin is flipped
    delE = (-Jnn * Jnn * delPSpin) + (Dp * delISpin);

    double r = rand1();

    if (delE == 0) {
        //if there is no change in the energy, coin flip to see if accepted
        if (r < 0.5) {
            spin[k].potts = p_new;
            spin[k].Sz = Sz_new;
            main_measurements.Ed_curr += Dp * delISpin;

            return 1;
        } else return 0;
    } else {
        if (r < exp(-delE * beta)) {
            spin[k].potts = p_new;
            spin[k].Sz = Sz_new;
            main_measurements.Ed_curr += Dp * delISpin;

            return 1;
        } else {
        return 0;
        }
    }
}

double MC_sweep(Measurements main_measurements, Interactions main_interactions, Properties main_properties, double beta, Site *spin, int Ns, int L) {
    int hits = 0;
    for(int i=0; i<Ns; i++) {
        hits += update_site(main_measurements, main_interactions, main_properties, i, beta, spin, Ns, L);
    }
    return ((double) hits)/((double) Ns);      // success rate
}

void Monte_Carlo_Sim(Interactions main_interactions, Measurements main_measurements, Properties main_properties, double beta, Site *spin, int Ns, int L, const string L_name) {

    int thermalize = 20000;
    int nsweep = 50;
    int ndata = 50000000;

    //Run 5000 sweeps of the system to achieve equilibrium
    double accepted = 0;
    
    for (int i = 0; i < thermalize; i++) {
        accepted += MC_sweep(main_measurements, main_interactions, main_properties, beta, spin, Ns, L);
    }
    cout << "spin update rate = " << accepted/((double) thermalize) << endl;

    //now start to collect data
    double E1 = 0, E2 = 0;
    double PM1 = 0, PM2 = 0, PM4 = 0;
    double IM1 = 0, IM2 = 0, IM4 = 0;
    double avg_accept = 0;

    int system_check = 100000;

    for (int n = 0; n < ndata; n++) {

        if(n % system_check == 0) cout << "n = " << n << endl;

        accepted = 0;
        for (int r = 0; r < nsweep; r++) {

            accepted += MC_sweep(main_measurements, main_interactions, main_properties, beta, spin, Ns, L);

        }
        accepted /= ((double) nsweep);
        avg_accept = (n * avg_accept + accepted) / (n + 1.);

        double e = main_measurements.clock_energy(main_properties, spin, Ns, L) + main_measurements.Ed_curr;

        E1 = (n * E1 + e) / (n + 1.);
        E2 = (n * E2 + e*e) / (n + 1.);

        double potts_mag = main_measurements.potts_magnetization(spin, Ns, L);

        PM1 = (n * PM1 + potts_mag) / (n + 1.);
        PM2 = (n * PM2 + pow(potts_mag, 2)) / (n + 1.);
        PM4 = (n * PM4 + pow(potts_mag, 4)) / (n + 1.);

        double ising_mag = main_measurements.ising_magnetization(spin, Ns);

        IM1 = (n * IM1 + ising_mag) / (n + 1.);
        IM2 = (n * IM2 + pow(ising_mag, 2)) / (n + 1.);
        IM4 = (n * IM4 + pow(ising_mag, 4)) / (n + 1.);
    }

    print(L_name, E1, E2, PM1, PM2, PM4, beta, Ns);
}

double rand1() {
    
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    double rand1 = dis(rng);
    
    return rand1;
}

void clear_files(const string L_name) {
    ofstream clear1, clear2, clear3, clear4, clear5;
    
    clear1.open(L_name + "_energy.dat", fstream::trunc);
    clear2.open(L_name + "_heat.dat", fstream::trunc);
    clear3.open(L_name + "_m.dat", fstream::trunc);
    clear4.open(L_name + "_susceptibility.dat", fstream::trunc);
    clear5.open(L_name + "_binder4.dat", fstream::trunc);

    clear1.close();
    clear2.close();
    clear3.close();
    clear4.close();
    clear5.close();

}

void print(const string L_name, double E1, double E2, double M1, double M2, double M4, double beta, int Ns) {
    
        ofstream energy;
        energy.open(L_name + "_energy.dat", fstream::app);
        ofstream heat;
        heat.open(L_name + "_heat.dat", fstream::app);
        ofstream m;
        m.open(L_name + "_m.dat", fstream::app);
        ofstream susceptibility;
        susceptibility.open(L_name + "_susceptibility.dat", fstream::app);
        ofstream binder4;
        binder4.open(L_name + "_binder4.dat", fstream::app);

        energy << 1./beta << ", " << E1/((double) Ns) << endl;
        heat << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) Ns) << endl;
        m << 1./beta << ", " << M1 << endl;
        susceptibility << 1./beta << ", " << beta*(M2 - M1*M1) << endl;
        binder4 << 1./beta << ", " << 1. - (M4/(3 * pow(M2, 2))) << endl;

        m.close();
        susceptibility.close();
        binder4.close();
        energy.close();
        heat.close();

}
