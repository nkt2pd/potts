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

const int N_nn1 = 4;

class Site {
public:
    int idx;    // index 0, 1, 2, ... Ns-1
    int x, y;   // coordinates
    int Sz;     // Ising spins
    int Ns;
    int L;
    // idx of its four neighbors
    Site *nn1[N_nn1];

};

double get_energy(Site *spin, int *Ns, int *L);
double get_magnetization(Site *spin, int *Ns, int *L);
void init_chain(Site *spin, int *Ns, int *L);
void set_coordinates(Site *spin, int *Ns, int *L);
void set_nn(Site *spin, int *Ns, int *L);
int mod_LR(int x, int m);
int mod_UD(int x, int m);
void Wolff_Cluster_Sim(double beta, Site *spin, int *Ns, int *L, int pick_file);
double Wolff_Sweep(double beta, Site *spin, int *Ns, int *L);
void randomize(Site *spin, int *Ns, int *L);
double rand1();
void clear_files();
void print(int i, double E1, double E2, double M1, double M2, double M4, double beta, int *Ns);
double theta(int p1, int p2);
void init_v_clock();

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

const int L1 = 4;       // chain size
const int L2 = 8;
const int L3 = 12;
const int L4 = 16;
const int L5 = 20;
const int L6 = 32;
const int L7 = 48;
const int L8 = 72;
const int D = 2;        // D-dimensional lattice
const int Ns1 = L1 * L1;      // number of spins
const int Ns2 = L2 * L2;
const int Ns3 = L3 * L3;
const int Ns4 = L4 * L4;
const int Ns5 = L5 * L5;
const int Ns6 = L6 * L6;
const int Ns7 = L7 * L7;
const int Ns8 = L8 * L8;
const string L1_name = "Potts4";
const string L2_name = "Potts8";
const string L3_name = "Potts12";
const string L4_name = "Potts16";
const string L5_name = "Potts20";
const string L6_name = "Potts32";
const string L7_name = "Potts48";
const string L8_name = "Potts72";


const double Jnn = -1;  // nearest-neighbor ferromagnetic coupling
const double q = 6; //number of potts states

// a simple data structure for square-lattice spin

class Cluster {
public:
    vector<int> sites;
    double size;
};

Site spin1[Ns1];
Site spin2[Ns2];
Site spin3[Ns3];
Site spin4[Ns4];
Site spin5[Ns5];
Site spin6[Ns6];
Site spin7[Ns7];
Site spin8[Ns8];

int main() {

    spin1[0].Ns = Ns1;
    spin1[0].L = L1;
    spin2[0].Ns = Ns2;
    spin2[0].L = L2;
    spin3[0].Ns = Ns3;
    spin3[0].L = L3;
    spin4[0].Ns = Ns4;
    spin4[0].L = L4;
    spin5[0].Ns = Ns5;
    spin5[0].L = L5;
    spin6[0].Ns = Ns6;
    spin6[0].L = L6;
    spin7[0].Ns = Ns7;
    spin7[0].L = L7;
    spin8[0].Ns = Ns8;
    spin8[0].L = L8;
        
    Site *spin;
    int *Ns;
    int *L;

    clear_files();
        
    for (int i = 0; i < 8; i++) {
        if (i == 0) {
            spin = spin1;
            Ns = &spin1[0].Ns;
            L = &spin1[0].L;
        } else if (i == 1) {
            spin = spin2;
            Ns = &spin2[0].Ns;
            L = &spin2[0].L;
        } else if (i == 2) {
            spin = spin3;
            Ns = &spin3[0].Ns;
            L = &spin3[0].L;
        } else if (i == 3) {
            spin = spin4;
            Ns = &spin4[0].Ns;
            L = &spin4[0].L;
        } else if (i == 4) {
            spin = spin5;
            Ns = &spin5[0].Ns;
            L = &spin5[0].L;
        } else if (i == 5) {
            spin = spin6;
            Ns = &spin6[0].Ns;
            L = &spin6[0].L;
        } else if (i == 6) {
            spin = spin7;
            Ns = &spin7[0].Ns;
            L = &spin7[0].L;
        } else {
            spin = spin8;
            Ns = &spin8[0].Ns;
            L = &spin8[0].L;
        }
            
        init_chain(spin, Ns, L);
            
        set_coordinates(spin, Ns, L);
        set_nn(spin, Ns, L);

        randomize(spin, Ns, L);

        double del_T = 0.1;

        cout << "Thank you for choosing the Potts Model :)" << endl;
            
        for(double T = 1.5; T>0; T -= del_T) {

            cout << "Lattice Length = " << *L << endl;
            cout << "T = " << T << endl;
            
            Wolff_Cluster_Sim(1./T, spin, Ns, L, i);
        }
    }
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
    
void init_chain(Site *spin, int *Ns, int *L) {
    for (int i=0; i < *Ns; i++) {
        spin[i].Sz = 1;
    }
}

void set_coordinates(Site *spin, int *Ns, int *L) {
    int y = 0;
    int x = 0;
    for(int i=0; i<spin[0].Ns; i++) {
        spin[i].idx = i;
        spin[i].x = x;
        spin[i].y = y;
        x++;
        if (x%*L == 0) {
            y++;
            x = 0;
        }
    }
}

void randomize(Site *spin, int *Ns, int *L) { //sets all site spins randomly to any of the q possible states
    for(int i=0; i<*Ns; i++) {
        double r = rand1();
        spin[i].Sz = (int)(r*6);
    } 
}

void set_nn(Site *spin, int *Ns, int *L) { //sets the 2*D nearest neighbors for each site in the system
    int k;
    for(int i=0; i<spin[0].Ns; i++) {
        
        k = mod_LR(spin[i].x + 1, *L);
        k += (*L * spin[i].y);
        spin[i].nn1[0] = &spin[k];
        
        k = mod_LR(spin[i].x - 1, *L);
        k += (*L * spin[i].y);
        spin[i].nn1[1] = &spin[k];
        
        k = mod_UD(i-*L, *Ns);
        spin[i].nn1[2] = &spin[k];

        k = mod_UD(i+*L, *Ns);
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

double get_magnetization(Site *spin, int *Ns, int *L) {
    
    double xsum = 0, ysum = 0, sum = 0;

    for (int i = 0; i < *Ns; i++) {
        xsum += cos(theta(spin[i].Sz, 0));
        ysum += sin(theta(spin[i].Sz, 0));
    }
    xsum /= (double) *Ns;
    ysum /= (double) *Ns;

    sum = sqrt((xsum*xsum) + (ysum*ysum));

    return sum;
}

double get_energy(Site *spin, int *Ns, int *L) {
    double sum = 0;
    for(int i=0; i<*Ns; i++) {
        for(int k=0; k<N_nn1; k++) {
            sum += V_clock[spin[i].Sz][spin[i].nn1[k]->Sz];
        }
    }
    return 0.5 * Jnn * sum;
}

double theta(int p1, int p2) {
    double theta;

    theta = (2*M_PI*((double)p1-(double)p2))/q;
    return theta;
}

double Wolff_Sweep(double beta, Site *spin, int *Ns, int *L) {
    double prob = 1 - exp(2*beta*Jnn);
    double r = rand1();
    double R = 0;
    r *= *L* (*L);

    if (r == 900.) r--;

    Cluster C;
    Cluster F_old;
    Cluster F_new;

    F_old.sites.clear();
    C.sites.clear();
    C.size = 0;
    C.sites.push_back((int)r);
    F_old.sites.push_back((int)r);

    vector<int> in_cluster(*Ns);

    in_cluster[(int)r] = 1;

    while(!F_old.sites.empty()) {
        F_new.sites.clear();

        for (int i = 0; i < F_old.sites.size(); i++) {

            for (int j = 0; j < N_nn1; j++) {

                int k = spin[F_old.sites[i]].nn1[j]->idx;

                if (spin[F_old.sites[i]].nn1[j]->Sz == spin[F_old.sites[i]].Sz && in_cluster[k] == 0) {

                    if (rand1() < prob) {

                        F_new.sites.push_back(spin[k].idx);
                        C.sites.push_back(spin[k].idx);
                        in_cluster[k] = 1;
                        C.size += 1.;

                    }
                }
            }
        }
            F_old.sites = F_new.sites;
    }
    
    //Once cluster is made, flip all the spins in cluster randomly
    R = rand1();

    for (int l = 0; l < C.sites.size(); l++) {

        spin[C.sites[l]].Sz = (int) (R*6);
    }

    return C.size;
}

void Wolff_Cluster_Sim(double beta, Site *spin, int *Ns, int *L, int pick_file) {

    int thermalize = 5000;
    int nsweep = 20;
    int ndata = 600000;

    //Run 5000 sweeps of the system to achieve equilibrium
    double total_members = 0;
    for (int i = 0; i < thermalize; i++) {
        total_members += Wolff_Sweep(beta, spin, Ns, L);
    }
    cout << "Average Cluster Size = " << total_members/((double) thermalize) << endl;

    //now start to collect data
    double E1 = 0, E2 = 0;
    double M1 = 0, M2 = 0, M4 = 0, abs_M = 0;
    double avg_size = 0;

    int system_check = 100000;

    for (int n = 0; n < ndata; n++) {

        if(n % system_check == 0) cout << "n = " << n << endl;

        total_members = 0;

        for (int r = 0; r < nsweep; r++) {
            total_members += Wolff_Sweep(beta, spin, Ns, L);
        }

        total_members /= ((double) nsweep);
        avg_size = (n * avg_size + total_members) / (n + 1.);

        double e = get_energy(spin, Ns, L);
        E1 = (n * E1 + e) / (n + 1.);
        E2 = (n * E2 + e*e) / (n + 1.);

        double mag = get_magnetization(spin, Ns, L);
        M1 = (n * M1 + mag) / (n + 1.);
        M2 = (n * M2 + pow(mag, 2)) / (n + 1.);
        M4 = (n * M4 + pow(mag, 4)) / (n + 1.);

    }

    print(pick_file, E1, E2, M1, M2, M4, beta, Ns);

}

double rand1() {
    
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    double rand1 = dis(rng);
    
    return rand1;
}

void clear_files() {
    ofstream clear1, clear2, clear3, clear4, clear5,
             clear6, clear7, clear8, clear9, clear10,
             clear11, clear12, clear13, clear14, clear15,
             clear16, clear17, clear18, clear19, clear20,
             clear21, clear22, clear23, clear24, clear25,
             clear26, clear27, clear28, clear29, clear30,
             clear31, clear32, clear33, clear34, clear35,
             clear36, clear37, clear38, clear39, clear40;
    
    clear1.open("Potts4_energy.dat", fstream::trunc);
    clear2.open("Potts4_heat.dat", fstream::trunc);
    clear3.open("Potts4_binder4.dat", fstream::trunc);
    clear4.open("Potts4_m.dat", fstream::trunc);
    clear5.open("Potts4_susceptibility.dat", fstream::trunc);
    clear6.open("Potts8_energy.dat", fstream::trunc);
    clear7.open("Potts8_heat.dat", fstream::trunc);
    clear8.open("Potts8_binder4.dat", fstream::trunc);
    clear9.open("Potts8_m.dat", fstream::trunc);
    clear10.open("Potts8_susceptibility.dat", fstream::trunc);
    clear11.open("Potts12_energy.dat", fstream::trunc);
    clear12.open("Potts12_heat.dat", fstream::trunc);
    clear13.open("Potts12_binder4.dat", fstream::trunc);
    clear14.open("Potts12_m.dat", fstream::trunc);
    clear15.open("Potts12_susceptibility.dat", fstream::trunc);
    clear16.open("Potts16_energy.dat", fstream::trunc);
    clear17.open("Potts16_heat.dat", fstream::trunc);
    clear18.open("Potts16_binder4.dat", fstream::trunc);
    clear19.open("Potts16_m.dat", fstream::trunc);
    clear20.open("Potts16_susceptibility.dat", fstream::trunc);
    clear21.open("Potts20_energy.dat", fstream::trunc);
    clear22.open("Potts20_heat.dat", fstream::trunc);
    clear23.open("Potts20_binder4.dat", fstream::trunc);
    clear24.open("Potts20_m.dat", fstream::trunc);
    clear25.open("Potts20_susceptibility.dat", fstream::trunc);
    clear26.open("Potts32_energy.dat", fstream::trunc);
    clear27.open("Potts32_heat.dat", fstream::trunc);
    clear28.open("Potts32_binder4.dat", fstream::trunc);
    clear29.open("Potts32_m.dat", fstream::trunc);
    clear30.open("Potts32_susceptibility.dat", fstream::trunc);
    clear31.open("Potts48_energy.dat", fstream::trunc);
    clear32.open("Potts48_heat.dat", fstream::trunc);
    clear33.open("Potts48_binder4.dat", fstream::trunc);
    clear34.open("Potts48_m.dat", fstream::trunc);
    clear35.open("Potts48_susceptibility.dat", fstream::trunc);
    clear36.open("Potts72_energy.dat", fstream::trunc);
    clear37.open("Potts72_heat.dat", fstream::trunc);
    clear38.open("Potts72_binder4.dat", fstream::trunc);
    clear39.open("Potts72_m.dat", fstream::trunc);
    clear40.open("Potts72_susceptibility.dat", fstream::trunc);

    clear1.close();
    clear2.close();
    clear3.close();
    clear4.close();
    clear5.close();
    clear6.close();
    clear7.close();
    clear8.close();
    clear9.close();
    clear10.close();
    clear11.close();
    clear12.close();
    clear13.close();
    clear14.close();
    clear15.close();
    clear16.close();
    clear17.close();
    clear18.close();
    clear19.close();
    clear20.close();
    clear21.close();
    clear22.close();
    clear23.close();
    clear24.close();
    clear25.close();
    clear26.close();
    clear27.close();
    clear28.close();
    clear29.close();
    clear30.close();
    clear31.close();
    clear32.close();
    clear33.close();
    clear34.close();
    clear35.close();
    clear36.close();
    clear37.close();
    clear38.close();
    clear39.close();
    clear40.close();
}

void print(int i, double E1, double E2, double M1, double M2, double M4, double beta, int *Ns) {
    if (i == 0) {
        ofstream energy;
        energy.open(L1_name + "_energy.dat", fstream::app);
        ofstream heat;
        heat.open(L1_name + "_heat.dat", fstream::app);
        ofstream m;
        m.open(L1_name + "_m.dat", fstream::app);
        ofstream susceptibility;
        susceptibility.open(L1_name + "_susceptibility.dat", fstream::app);
        ofstream binder4;
        binder4.open(L1_name + "_binder4.dat", fstream::app);

        energy << 1./beta << ", " << E1/((double) *Ns) << endl;
        heat << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) *Ns) << endl;
        m << 1./beta << ", " << M1 << endl;
        susceptibility << 1./beta << ", " << beta*(M2 - M1*M1)/((double) *Ns) << endl;
        binder4 << 1./beta << ", " << 1. - (M4/(3 * pow(M2, 2))) << endl;

        m.close();
        susceptibility.close();
        binder4.close();
        energy.close();
        heat.close();

    } else if (i == 1) {
        ofstream energy;
        energy.open(L2_name + "_energy.dat", fstream::app);
        ofstream heat;
        heat.open(L2_name + "_heat.dat", fstream::app);
        ofstream m;
        m.open(L2_name + "_m.dat", fstream::app);
        ofstream susceptibility;
        susceptibility.open(L2_name + "_susceptibility.dat", fstream::app);
        ofstream binder4;
        binder4.open(L2_name + "_binder4.dat", fstream::app);

        energy << 1./beta << ", " << E1/((double) *Ns) << endl;
        heat << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) *Ns) << endl;
        m << 1./beta << ", " << M1/((double) *Ns) << endl;
        susceptibility << 1./beta << ", " << beta*(M2 - M1*M1)/((double) *Ns) << endl;
        binder4 << 1./beta << ", " << 1. - (M4/(3 * pow(M2, 2))) << endl;

        m.close();
        susceptibility.close();
        binder4.close();
        energy.close();
        heat.close();

    } else if (i == 2) {
        ofstream energy;
        energy.open(L3_name + "_energy.dat", fstream::app);
        ofstream heat;
        heat.open(L3_name + "_heat.dat", fstream::app);
        ofstream m;
        m.open(L3_name + "_m.dat", fstream::app);
        ofstream susceptibility;
        susceptibility.open(L3_name + "_susceptibility.dat", fstream::app);
        ofstream binder4;
        binder4.open(L3_name + "_binder4.dat", fstream::app);

        energy << 1./beta << ", " << E1/((double) *Ns) << endl;
        heat << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) *Ns) << endl;
        m << 1./beta << ", " << M1/((double) *Ns) << endl;
        susceptibility << 1./beta << ", " << beta*(M2 - M1*M1)/((double) *Ns) << endl;
        binder4 << 1./beta << ", " << 1. - (M4/(3 * pow(M2, 2))) << endl;

        m.close();
        susceptibility.close();
        binder4.close();
        energy.close();
        heat.close();

    } else if (i == 3) {
        ofstream energy;
        energy.open(L4_name + "_energy.dat", fstream::app);
        ofstream heat;
        heat.open(L4_name + "_heat.dat", fstream::app);
        ofstream m;
        m.open(L4_name + "_m.dat", fstream::app);
        ofstream susceptibility;
        susceptibility.open(L4_name + "_susceptibility.dat", fstream::app);
        ofstream binder4;
        binder4.open(L4_name + "_binder4.dat", fstream::app);

        energy << 1./beta << ", " << E1/((double) *Ns) << endl;
        heat << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) *Ns) << endl;
        m << 1./beta << ", " << M1/((double) *Ns) << endl;
        susceptibility << 1./beta << ", " << beta*(M2 - M1*M1)/((double) *Ns) << endl;
        binder4 << 1./beta << ", " << 1. - (M4/(3 * pow(M2, 2))) << endl;

        m.close();
        susceptibility.close();
        binder4.close();
        energy.close();
        heat.close();

    } else if (i == 4) {
        ofstream energy;
        energy.open(L5_name + "_energy.dat", fstream::app);
        ofstream heat;
        heat.open(L5_name + "_heat.dat", fstream::app);
        ofstream m;
        m.open(L5_name + "_m.dat", fstream::app);
        ofstream susceptibility;
        susceptibility.open(L5_name + "_susceptibility.dat", fstream::app);
        ofstream binder4;
        binder4.open(L5_name + "_binder4.dat", fstream::app);

        energy << 1./beta << ", " << E1/((double) *Ns) << endl;
        heat << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) *Ns) << endl;
        m << 1./beta << ", " << M1/((double) *Ns) << endl;
        susceptibility << 1./beta << ", " << beta*(M2 - M1*M1)/((double) *Ns) << endl;
        binder4 << 1./beta << ", " << 1. - (M4/(3 * pow(M2, 2))) << endl;

        m.close();
        susceptibility.close();
        binder4.close();
        energy.close();
        heat.close();

    } else if (i == 5 ) {
        ofstream energy;
        energy.open(L6_name + "_energy.dat", fstream::app);
        ofstream heat;
        heat.open(L6_name + "_heat.dat", fstream::app);
        ofstream m;
        m.open(L6_name + "_m.dat", fstream::app);
        ofstream susceptibility;
        susceptibility.open(L6_name + "_susceptibility.dat", fstream::app);
        ofstream binder4;
        binder4.open(L6_name + "_binder4.dat", fstream::app);

        energy << 1./beta << ", " << E1/((double) *Ns) << endl;
        heat << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) *Ns) << endl;
        m << 1./beta << ", " << M1/((double) *Ns) << endl;
        susceptibility << 1./beta << ", " << beta*(M2 - M1*M1)/((double) *Ns) << endl;
        binder4 << 1./beta << ", " << 1. - (M4/(3 * pow(M2, 2))) << endl;

        m.close();
        susceptibility.close();
        binder4.close();
        energy.close();
        heat.close();

    } else if (i == 6) {
        ofstream energy;
        energy.open(L7_name + "_energy.dat", fstream::app);
        ofstream heat;
        heat.open(L7_name + "_heat.dat", fstream::app);
        ofstream m;
        m.open(L7_name + "_m.dat", fstream::app);
        ofstream susceptibility;
        susceptibility.open(L7_name + "_susceptibility.dat", fstream::app);
        ofstream binder4;
        binder4.open(L7_name + "_binder4.dat", fstream::app);

        energy << 1./beta << ", " << E1/((double) *Ns) << endl;
        heat << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) *Ns) << endl;
        m << 1./beta << ", " << M1/((double) *Ns) << endl;
        susceptibility << 1./beta << ", " << beta*(M2 - M1*M1)/((double) *Ns) << endl;
        binder4 << 1./beta << ", " << 1. - (M4/(3 * pow(M2, 2))) << endl;

        m.close();
        susceptibility.close();
        binder4.close();
        energy.close();
        heat.close();

    } else {
        ofstream energy;
        energy.open(L8_name + "_energy.dat", fstream::app);
        ofstream heat;
        heat.open(L8_name + "_heat.dat", fstream::app);
        ofstream m;
        m.open(L8_name + "_m.dat", fstream::app);
        ofstream susceptibility;
        susceptibility.open(L8_name + "_susceptibility.dat", fstream::app);
        ofstream binder4;
        binder4.open(L8_name + "_binder4.dat", fstream::app);

        energy << 1./beta << ", " << E1/((double) *Ns) << endl;
        heat << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) *Ns) << endl;
        m << 1./beta << ", " << M1/((double) *Ns) << endl;
        susceptibility << 1./beta << ", " << beta*(M2 - M1*M1)/((double) *Ns) << endl;
        binder4 << 1./beta << ", " << 1. - (M4/(3 * pow(M2, 2))) << endl;

        m.close();
        susceptibility.close();
        binder4.close();
        energy.close();
        heat.close();

    }
}