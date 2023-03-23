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

using namespace std;

const int N_nn1 = 4;    //Number of Nearest neighbors

class Site {
public:
    int idx;    // index 0, 1, 2, ... Ns-1
    int x, y;   // coordinates
    int Sz;     // Ising spins
    int potts;
    int Ns;
    int L;
    // idx of its four neighbors
    Site *nn1[N_nn1];

};

double get_energy(Site *spin, int Ns, int L);
double get_magnetization(Site *spin, int Ns, int L);
void init_chain(Site *spin, int Ns, int L);
void set_coordinates(Site *spin, int Ns, int L);
void set_nn(Site *spin, int Ns, int L);
int mod_LR(int x, int m);
int mod_UD(int x, int m);
//void Wolff_Cluster_Sim(double beta, Site *spin, int *Ns, int *L, int pick_file);
//double Wolff_Sweep(double beta, Site *spin, int *Ns, int *L);
void randomize(Site *spin, int Ns, int L);
double rand1();
void clear_files(const string L_name);
void print(const string L_name, double E1, double E2, double M1, double M2, double M4, double beta, int Ns);
double theta(int p1, int p2);
void init_v_clock();
int update_site(int k, double beta, Site *spin, int Ns, int L);
double MC_sweep(double beta, Site *spin, int Ns, int L);
void Monte_Carlo_Sim(double beta, Site *spin, int Ns, int L, const string L_name);

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());


const int D = 2;        // D-dimensional lattice

const double Jnn = -1;  // nearest-neighbor ferromagnetic coupling
const double q = 6; //number of potts states

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

    const string L_size(argv[1]);
    const string L_name = "Potts" + L_size;
    
    clear_files(L_name);
    init_v_clock();
    
    init_chain(spin, Ns, L);

    set_coordinates(spin, Ns, L);
    set_nn(spin, Ns, L);
    
    randomize(spin, Ns, L);

    double del_T = 0.1;

    cout << "Thank you for choosing the Potts Model :)" << endl;
            
    for(double T = 1.5; T>0; T -= del_T) {

        cout << "Lattice Length = " << L << endl;
        cout << "T = " << T << endl;
            
        Monte_Carlo_Sim(1./T, spin, Ns, L, L_name);
    }


    delete[] spin;

    return 0;
} 

double V_clock[6][6];

void init_v_clock() {

    for (int i = 0; i < q; i++) {
        for (int j = 0; j < q; j++) {
            
            V_clock[i][j] = cos((i-j)*2.*M_PI/((double) q));
        }
    }
}
    
void init_chain(Site *spin, int Ns, int L) {
    for (int i=0; i < Ns; i++) {
        spin[i].Sz = 1;
        spin[i].potts = 1;
    }
}

void set_coordinates(Site *spin, int Ns, int L) {
    int y = 0;
    int x = 0;
    for(int i=0; i<Ns; i++) {
        spin[i].idx = i;
        spin[i].x = x;
        spin[i].y = y;
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
    for(int i=0; i<Ns; i++) {
        
        k = mod_LR(spin[i].x + 1, L);
        k += (L * spin[i].y);
        spin[i].nn1[0] = &spin[k];
        
        k = mod_LR(spin[i].x - 1, L);
        k += (L * spin[i].y);
        spin[i].nn1[1] = &spin[k];
        
        k = mod_UD(i-L, Ns);
        spin[i].nn1[2] = &spin[k];

        k = mod_UD(i+L, Ns);
        spin[i].nn1[3] = &spin[k];

    }
}

int mod_LR(int x, int m) { //implements periodic boundary condition in x-direction
    if (x>=0 && x<m) {
        return x;
    } else if (x<0) {
        return m-1-mod_LR(-1-x,m);
    } else {
        return x%m;
    }
}

int mod_UD(int x, int m) { //implements periodic boundary condition in y-direction
    if (x>=0 && x<m)
        return x;
    else if (x<0)
        return m-1-mod_UD(-1-x,m);
    else
        return x%m;
}

double get_magnetization(Site *spin, int Ns, int L) {
    
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

double get_energy(Site *spin, int Ns, int L) {
    double sum = 0;
    for(int i=0; i<Ns; i++) {
        for(int k=0; k<N_nn1; k++) {
            sum += V_clock[spin[i].potts][spin[i].nn1[k]->potts];
        }
    }
    return 0.5 * Jnn * sum;
}

double theta(int p1, int p2) {
    double theta;

    theta = (2*M_PI*((double)p1-(double)p2))/q;
    return theta;
}

//test a new state, accept or reject the state based on parameters. 
int update_site(int k, double beta, Site *spin, int Ns, int L) {

    int p_new = (int) (rand1() * q);
    double delS = 0;
    double delE = 0;

    for(int l=0; l<N_nn1; l++) {
        delS += (V_clock[p_new][spin[k].nn1[l]->potts] - V_clock[spin[k].potts][spin[k].nn1[l]->potts]);
    }

    double Hamiltonian = Jnn * delS;

    //change in energy when spin is flipped (should there be a negative on spin[k].Sz?)
    delE = -Jnn * Hamiltonian;

    double r = rand1();

    if (delE == 0) {
        //if there is no change in the energy, coin flip to see if accepted
        if (r < 0.5) {
            spin[k].potts = p_new;
            return 1;
        } else return 0;
    } else {
        if (r < exp(-delE * beta)) {
            spin[k].potts = p_new;
            return 1;
        } else {
        return 0;
        }
    }
}

double MC_sweep(double beta, Site *spin, int Ns, int L) {
    int hits = 0;
    for(int i=0; i<Ns; i++) {
        hits += update_site(i, beta, spin, Ns, L);
    }
    return ((double) hits)/((double) Ns);      // success rate
}

void Monte_Carlo_Sim(double beta, Site *spin, int Ns, int L, const string L_name) {

    int thermalize = 20000;
    int nsweep = 50;
    int ndata = 500000;

    //Run 5000 sweeps of the system to achieve equilibrium
    double accepted = 0;
    
    for (int i = 0; i < thermalize; i++) {
        accepted += MC_sweep(beta, spin, Ns, L);
    }
    cout << "spin update rate = " << accepted/((double) thermalize) << endl;

    //now start to collect data
    double E1 = 0, E2 = 0;
    double M1 = 0, M2 = 0, M4 = 0;
    double avg_accept = 0;

    int system_check = 100000;

    for (int n = 0; n < ndata; n++) {

        if(n % system_check == 0) cout << "n = " << n << endl;

        accepted = 0;
        for (int r = 0; r < nsweep; r++) {

            accepted += MC_sweep(beta, spin, Ns, L);

        }
        accepted /= ((double) nsweep);
        avg_accept = (n * avg_accept + accepted) / (n + 1.);

        double e = get_energy(spin, Ns, L);

        E1 = (n * E1 + e) / (n + 1.);
        E2 = (n * E2 + e*e) / (n + 1.);

        double mag = get_magnetization(spin, Ns, L);

        M1 = (n * M1 + mag) / (n + 1.);
        M2 = (n * M2 + pow(mag, 2)) / (n + 1.);
        M4 = (n * M4 + pow(mag, 4)) / (n + 1.);
    }

    print(L_name, E1, E2, M1, M2, M4, beta, Ns);
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