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
    std::ifstream infile1, infile2, infile3, infile4, infile5, infile6, infile7, infile8;

    std::ofstream ofile1, ofile2, ofile3, ofile4, ofile5, ofile6, ofile7, ofile8,
                  ofile9, ofile10, ofile11, ofile12, ofile13, ofile14, ofile15,
                  ofile16, ofile17, ofile18, ofile19, ofile20;

    infile1.open("Potts4_m.dat");
    infile2.open("Potts8_m.dat");
    infile3.open("Potts12_m.dat");
    infile4.open("Potts16_m.dat");
    infile5.open("Potts20_m.dat");
    infile6.open("Potts32_m.dat");
    infile7.open("Potts48_m.dat");
    infile8.open("Potts72_m.dat");

    ofile1.open("t=1.5_m.vs.L_fss.dat");
    ofile2.open("t=1.4_m.vs.L_fss.dat");
    ofile3.open("t=1.3_m.vs.L_fss.dat");
    ofile4.open("t=1.2_m.vs.L_fss.dat");
    ofile5.open("t=1.1_m.vs.L_fss.dat");
    ofile6.open("t=1.05_m.vs.L_fss.dat");
    ofile7.open("t=1.0_m.vs.L_fss.dat");
    ofile8.open("t=0.95_m.vs.L_fss.dat");
    ofile9.open("t=0.9_m.vs.L_fss.dat");
    ofile10.open("t=0.85_m.vs.L_fss.dat");
    ofile11.open("t=0.8_m.vs.L_fss.dat");
    ofile12.open("t=0.75_m.vs.L_fss.dat");
    ofile13.open("t=0.7_m.vs.L_fss.dat");
    ofile14.open("t=0.65_m.vs.L_fss.dat");
    ofile15.open("t=0.55_m.vs.L_fss.dat");
    ofile16.open("t=0.45_m.vs.L_fss.dat");
    ofile17.open("t=0.35_m.vs.L_fss.dat");
    ofile18.open("t=0.25_m.vs.L_fss.dat");
    ofile19.open("t=0.15_m.vs.L_fss.dat");
    ofile20.open("t=0.05_m.vs.L_fss.dat");

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
    streams.push_back(std::move(ofile16));
    streams.push_back(std::move(ofile17));
    streams.push_back(std::move(ofile18));
    streams.push_back(std::move(ofile19));
    streams.push_back(std::move(ofile20));


    for (int i = 0; i < *num_temps; i++) {
        infile1 >> buffer;
        infile2 >> buffer;
        infile3 >> buffer;
        infile4 >> buffer;
        infile5 >> buffer;
        infile6 >> buffer;
        infile7 >> buffer;
        infile8 >> buffer;

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

        infile7 >> buffer;
        vals = stod(buffer);
        streams[i] << "48, " << vals << std::endl;

        infile8 >> buffer;
        vals = stod(buffer);
        streams[i] << "72, " << vals << std::endl;
    }

    infile1.close();
    infile2.close();
    infile3.close();
    infile4.close();
    infile5.close();
    infile6.close();
    infile7.close();
    infile8.close();

    ofile1.close();
    ofile2.close();
    ofile3.close();
    ofile4.close();
    ofile5.close();
    ofile6.close();
    ofile7.close();
    ofile8.close();
    ofile9.close();
    ofile10.close();
    ofile11.close();
    ofile12.close();
    ofile13.close();
    ofile14.close();
    ofile15.close();
    ofile16.close();
    ofile17.close();
    ofile18.close();
    ofile19.close();
    ofile20.close();



}

/*
File format

temp, m

in each file, want one m value for each value of L. Each file represents one temperature.
*/