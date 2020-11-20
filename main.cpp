#include "IsingModel.hpp"
#include <iostream>
#include <armadillo>
#include <fstream>
#include <omp.h>
#include <string>

using namespace std;
using namespace arma;

// output file
ofstream ofile;

// Main program begins here

int main(int argc, char* argv[])
{
  string filename, order;
  int NSpins, MCcycles;
  double InitialTemp, FinalTemp, TempStep;
  if (argc <= 5) {
    cout << "Bad Usage: " << argv[0] <<
      " read output file, Number of spins, MC cycles, initial and final temperature and tempurate step" << endl;
    exit(1);
  }
  if (argc > 1) {
    filename=argv[1];
    NSpins = atoi(argv[2]);
    MCcycles = atoi(argv[3]);
    InitialTemp = atof(argv[4]);
    FinalTemp = atof(argv[5]);
    TempStep = atof(argv[6]);
    order = argv[7];
  }
  // Declare new file name and add lattice size to file name
  string fileout = "./textfiles/"+filename;
  string argument = to_string(NSpins);
  fileout.append(argument);
  ofile.open(fileout);



  if(filename=="Lattice"){
  // Start Monte Carlo sampling by looping over the selcted Temperatures
  #pragma omp_parallel for
  for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature+=TempStep){

    IsingModel *current_temp_run;
    current_temp_run = new IsingModel(Temperature,NSpins,MCcycles, order);
    // Start Monte Carlo computation and get expectation values
    current_temp_run->MetropolisSampling();
    current_temp_run->WriteResultstoFile(ofile);




  }
  }
  if(filename=="Energy"){
    double Temperature = InitialTemp;

    IsingModel *energy_run;
    energy_run = new IsingModel(Temperature,NSpins,MCcycles, order);
    // Start Monte Carlo computation and get expectation values
    energy_run->Get_Energy();
    vec energy = energy_run->Energy_vector;
    vec magneticMoment = energy_run->MagneticMoment_vector;
    vec number_of_flips = energy_run->number_of_flips_vector;



    ofile << "MC_cycles=" << MCcycles
    << " T=" << Temperature << " number_of_spins="<< NSpins*NSpins <<
    " order=" << order << endl;

    //ofile << "variance=" << energy_run->variance << endl;

    ofile <<" Energy-vector: MagnetcMoment-vector: Number of flips:" << endl;
    for(int i=0;i<=MCcycles;i++){
      ofile << energy(i) << " " << magneticMoment(i) << " " << number_of_flips(i) << endl;
    }
  }


  ofile.close();  // close output file
  return 0;
}
