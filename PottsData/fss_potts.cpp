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
#include <string>
#include <algorithm>

/* 
Plan: File for each temperature, scan one line at a time from each lattice size file
      into the corresponding temperature file for each temperature
*/

int how_many_temps();
void make_files(int *num_temps);

int main() {

    int num_temps = 0;
    num_temps = how_many_temps();

    make_files(&num_temps);

    return 0;
}

int how_many_temps() {
    int num_temps = 0;

    std::ifstream arbFile;
    std::string line;
    arbFile.open("Potts4_m.dat");

    while (std::getline(arbFile, line)) {
        num_temps++;
    }

    return num_temps;

}

void make_files(int *num_temps) {
    std::ifstream infile1, infile2, infile3, infile4, infile5, infile6;

    std::ofstream ofile1, ofile2, ofile3, ofile4, ofile5, ofile6, ofile7, ofile8,
                  ofile9, ofile10, ofile11, ofile12, ofile13, ofile14, ofile15;

    infile1.open("Potts4_m.dat");
    infile2.open("Potts8_m.dat");
    infile3.open("Potts12_m.dat");
    infile4.open("Potts16_m.dat");
    infile5.open("Potts20_m.dat");
    infile6.open("Potts32_m.dat");

    ofile1.open("t=1.5_m.vs.L_fss.dat");
    ofile2.open("t=1.4_m.vs.L_fss.dat");
    ofile3.open("t=1.3_m.vs.L_fss.dat");
    ofile4.open("t=1.2_m.vs.L_fss.dat");
    ofile5.open("t=1.1_m.vs.L_fss.dat");
    ofile6.open("t=1.0_m.vs.L_fss.dat");
    ofile7.open("t=0.9_m.vs.L_fss.dat");
    ofile8.open("t=0.8_m.vs.L_fss.dat");
    ofile9.open("t=0.7_m.vs.L_fss.dat");
    ofile10.open("t=0.6_m.vs.L_fss.dat");
    ofile11.open("t=0.5_m.vs.L_fss.dat");
    ofile12.open("t=0.4_m.vs.L_fss.dat");
    ofile13.open("t=0.3_m.vs.L_fss.dat");
    ofile14.open("t=0.2_m.vs.L_fss.dat");
    ofile15.open("t=0.1_m.vs.L_fss.dat");

    std::string buffer;
    double vals = 0;

    std::vector<std::ofstream> streams;
    streams.push_back(std::move(ofile1));
    streams.push_back(std::move(ofile2));
    streams.push_back(std::move(ofile3));
    streams.push_back(std::move(ofile4));
    streams.push_back(std::move(ofile5));
    streams.push_back(std::move(ofile6));
    streams.push_back(std::move(ofile7));
    streams.push_back(std::move(ofile8));
    streams.push_back(std::move(ofile9));
    streams.push_back(std::move(ofile10));
    streams.push_back(std::move(ofile11));
    streams.push_back(std::move(ofile12));
    streams.push_back(std::move(ofile13));
    streams.push_back(std::move(ofile14));
    streams.push_back(std::move(ofile15));


    for (int i = 0; i < *num_temps; i++) {
        infile1 >> buffer;
        infile2 >> buffer;
        infile3 >> buffer;
        infile4 >> buffer;
        infile5 >> buffer;
        infile6 >> buffer;

        infile1 >> buffer;
        vals = stod(buffer);
        streams[i] << "4, " << vals << std::endl;

        infile2 >> buffer;
        vals = stod(buffer);
        streams[i] << "8, " << vals << std::endl;

        infile3 >> buffer;
        vals = stod(buffer);
        streams[i] << "12, " << vals << std::endl;

        infile4 >> buffer;
        vals = stod(buffer);
        streams[i] << "16, " << vals << std::endl;

        infile5 >> buffer;
        vals = stod(buffer);
        streams[i] << "20, " << vals << std::endl;

        infile6 >> buffer;
        vals = stod(buffer);
        streams[i] << "32, " << vals << std::endl;
    }

}

/*
File format

temp, L, m

in each file, want one m value for each value of L. Each file represents one temperature.
*/