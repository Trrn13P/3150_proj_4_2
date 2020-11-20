#include <armadillo>
#include <fstream>
#include <iostream>

using namespace arma;
using namespace std;

#ifndef ISINGMODEL_HPP
#define ISINGMODEL_HPP

class IsingModel {
private:
  /* data */

public:
  void MetropolisSampling();
  void InitializeLattice();
  void WriteResultstoFile(ofstream &ofile);

  void Get_Energy();

  float random_between_zero_and_one();
  int minus_one_or_one();


  double Temperature;
  int NSpins, MCcycles;

  int number_of_flips;

  mat SpinMatrix;
  vec Energy_vector;
  vec MagneticMoment_vector;
  vec number_of_flips_vector;
  vec EnergyDifference = zeros<mat>(17);
  vec ExpectationValues = zeros<mat>(5);

  double Energy, MagneticMoment;

  string order;

  IsingModel(double Temperature_, int NSpins_, int MCcycles_, string order_){
    order = order_;
    NSpins = NSpins_;
    MCcycles = MCcycles_;
    Temperature = Temperature_;

    // setup array for possible energy changes
    for( int de =-8; de <= 8; de+=4) EnergyDifference(de+8) = exp(-de/Temperature);

    // Initialize the lattice spin values
    SpinMatrix = zeros<mat>(NSpins,NSpins);
    Energy_vector = zeros<vec>(MCcycles+1);
    MagneticMoment_vector = zeros<vec>(MCcycles+1);
    number_of_flips_vector = zeros<vec>(MCcycles+1);

    //    initialize energy and magnetization
    double Energy = 0.;     double MagneticMoment = 0.;

  }
};
#endif
