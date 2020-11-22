#include "IsingModel.hpp"
#include <iostream>
#include <armadillo>
#include <fstream>
#include <omp.h>
#include <string>
#include <time.h>

using namespace std;
using namespace arma;

// output file
ofstream ofile;

//Based on argv arguments one of these will run, explained in readme.
void parallelization_test(int NSpins, int MCcycles, double InitialTemp, double FinalTemp, double TempStep, string order, ofstream &ofile);
void Lattice(int NSpins, int MCcycles, double InitialTemp, double FinalTemp, double TempStep, string order, ofstream &ofile);
void Energy(int NSpins, int MCcycles, double InitialTemp, double FinalTemp, double TempStep, string order, ofstream &ofile);

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
    //The order can be Unordered or Ordered, if Unordered the SpinMatrix will be given random -1 or 1,
    //if Ordered all spins will be up
    order = argv[7];
  }
  // Declare new file name and add lattice size to file name
  string fileout = "./textfiles/"+filename;
  string argument = to_string(NSpins);
  fileout.append(argument);
  ofile.open(fileout);

  //Running functions based on filenames
  if(filename=="Lattice"){
    Lattice(NSpins,MCcycles,InitialTemp,FinalTemp,TempStep,order,ofile);
  }
  if(filename=="Energy"){
    Energy(NSpins,MCcycles,InitialTemp,FinalTemp,TempStep,order,ofile);
  }

  if(filename=="parallelization_test"){
    parallelization_test(NSpins,MCcycles,InitialTemp,FinalTemp,TempStep, order, ofile);
}

  ofile.close();  // close output file
  return 0;
}

//This function saves the energy, magneticMoment and #of flips as function MCcycle
//It will run with the Initial temperature, the fifth argv  (argv[4]).
void Energy(int NSpins, int MCcycles, double InitialTemp, double FinalTemp, double TempStep, string order, ofstream &ofile)
{
  double Temperature = InitialTemp;

  //Making class pointer and running for one temperature
  IsingModel *energy_run;
  energy_run = new IsingModel(Temperature,NSpins,MCcycles, order);
  //Running the MC simulation
  energy_run->Get_Energy();

  //extracting data from class
  vec energy = energy_run->Energy_vector;
  vec magneticMoment = energy_run->MagneticMoment_vector;
  vec number_of_flips = energy_run->number_of_flips_vector;

  //Writing the data to class
  ofile << "MC_cycles=" << MCcycles
  << " T=" << Temperature << " number_of_spins="<< NSpins*NSpins <<
  " order=" << order << endl;

  ofile <<" Energy-vector: MagnetcMoment-vector: Number of flips:" << endl;
  for(int i=0;i<=MCcycles;i++){
    ofile << energy(i) << " " << magneticMoment(i) << " " << number_of_flips(i) << endl;
  }
}

//This function runs a MC simulation between an initial temp (argv[4]) and a final temp (argv[5]) with
//a temperature step length (argv[6])
void Lattice(int NSpins, int MCcycles, double InitialTemp, double FinalTemp, double TempStep, string order, ofstream &ofile)
{
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


//This function is for testing the parallelization times. It runs 10 times unparallelized and 10 times paralleizied,
//and saves the data to file.
//It runs the same way as the Lattice function except not saving the expectation values ect. to file.
void parallelization_test(int NSpins, int MCcycles, double InitialTemp, double FinalTemp, double TempStep, string order, ofstream &ofile){
  double start, finish, runtime;
  ofile << "MC_cycles=" << MCcycles << " L=" << NSpins << endl;
  ofile << "t-unparellized: t-paralleized:" << endl;

  for(int i=0;i<10;i++){
    start = clock();
    for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature+=TempStep){
      IsingModel *current_temp_run;
      current_temp_run = new IsingModel(Temperature,NSpins,MCcycles, order);
      current_temp_run->MetropolisSampling();
    }
    finish = clock();
    runtime = (finish -start)/CLOCKS_PER_SEC;
    ofile << runtime << " ";

  start = clock();
  #pragma omp_parallel for
  for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature+=TempStep){
    IsingModel *current_temp_run;
    current_temp_run = new IsingModel(Temperature,NSpins,MCcycles, order);
    // Start Monte Carlo computation and get expectation values
    current_temp_run->MetropolisSampling();
  }
  finish = clock();
  runtime = 0.5*(finish -start)/CLOCKS_PER_SEC;
  ofile << runtime << endl;
}
}
