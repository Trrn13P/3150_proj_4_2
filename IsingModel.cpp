#include "IsingModel.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
using namespace  std;
using namespace arma;


// inline function for PeriodicBoundary boundary conditions
inline int PeriodicBoundary(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}



/*
returns a random number between 0 and 1.
*/
float IsingModel::random_between_zero_and_one(){
  return rand()*1./RAND_MAX;
}

/*
This function returns a -1 or 1 randomly
*/
int IsingModel::minus_one_or_one(){
  float number = random_between_zero_and_one();
  return 2*round(number)-1;
}


// The Monte Carlo part with the Metropolis algo with sweeps over the lattice
void IsingModel::MetropolisSampling()
{
  // Initialize the seed and call the Mersienne algo
  std::random_device rd;
  std::mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in [[0, 1]
  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);


  // initialize SpinMatrix
  InitializeLattice();
  // Start Monte Carlo cycles
  for (int cycles = 1; cycles <= MCcycles; cycles++){
    // The sweep over the lattice, looping over all spin sites
    for(int x =0; x < NSpins; x++) {
      for (int y= 0; y < NSpins; y++){
	int ix = (int) (RandomNumberGenerator(gen)*(double)NSpins);
	int iy = (int) (RandomNumberGenerator(gen)*(double)NSpins);
	int deltaE =  2*SpinMatrix(ix,iy)*
	  (SpinMatrix(ix,PeriodicBoundary(iy,NSpins,-1))+
	   SpinMatrix(PeriodicBoundary(ix,NSpins,-1),iy) +
	   SpinMatrix(ix,PeriodicBoundary(iy,NSpins,1)) +
	   SpinMatrix(PeriodicBoundary(ix,NSpins,1),iy));
     //Metropolis test
	if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) ) {
	  SpinMatrix(ix,iy) *= -1.0;  // flip one spin and accept new spin config
	  MagneticMoment += (double) 2*SpinMatrix(ix,iy);
	  Energy += (double) deltaE;
	}
      }
    }
    // update expectation values  for local node
    ExpectationValues(0) += Energy;    ExpectationValues(1) += Energy*Energy;
    ExpectationValues(2) += MagneticMoment;
    ExpectationValues(3) += MagneticMoment*MagneticMoment;
    ExpectationValues(4) += fabs(MagneticMoment);
  }
} // end of Metropolis sampling over spins

// function to initialise energy, spin matrix and magnetization
void IsingModel::InitializeLattice()
{
  // setup spin matrix and initial magnetization
  //order is argv[7]
  if(order=="Ordered"){
    for(int x =0; x < NSpins; x++) {
      for (int y= 0; y < NSpins; y++){
        SpinMatrix(x,y) = 1.0; // all spins up
        MagneticMoment +=  (double) SpinMatrix(x,y);
      }
    }
  }
  else if(order=="Unordered"){
    for(int x =0; x < NSpins; x++) {
      for (int y= 0; y < NSpins; y++){
        SpinMatrix(x,y) = minus_one_or_one(); // random spin orientation
        MagneticMoment +=  (double) SpinMatrix(x,y);
      }
    }
  }
  else{
    //Prints error message in terminal if provided bad argv[7]
    cout << "Wrong Ordered/Unordered" << endl;
    exit(1);
  }

  // setup initial energy
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      Energy -=  (double) SpinMatrix(x,y)*
	(SpinMatrix(PeriodicBoundary(x,NSpins,-1),y) +
	 SpinMatrix(x,PeriodicBoundary(y,NSpins,-1)));
    }
  }
}// end function initialise


//Not changed from Morten's program, writes expectation values to file.
void IsingModel::WriteResultstoFile(ofstream &ofile)
{
  double norm = 1.0/((double) (MCcycles));  // divided by  number of cycles
  double E_ExpectationValues = ExpectationValues(0)*norm;
  double E2_ExpectationValues = ExpectationValues(1)*norm;
  double M_ExpectationValues = ExpectationValues(2)*norm;
  double M2_ExpectationValues = ExpectationValues(3)*norm;
  double Mabs_ExpectationValues = ExpectationValues(4)*norm;
  // all expectation values are per spin, divide by 1/NSpins/NSpins
  double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/NSpins/NSpins;
  double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/NSpins/NSpins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << Temperature;
  ofile << setw(15) << setprecision(8) << E_ExpectationValues/NSpins/NSpins;
  ofile << setw(15) << setprecision(8) << Evariance/Temperature/Temperature;
  ofile << setw(15) << setprecision(8) << M_ExpectationValues/NSpins/NSpins;
  ofile << setw(15) << setprecision(8) << Mvariance/Temperature;
  ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues/NSpins/NSpins << endl;
} // end output function


//This function is almost identical to the MetropolisSampling function except
//it doesn't save expectation values, but instead saves MagneticMoment, # of flips and Energy
//for each MC-cycle.
void IsingModel::Get_Energy()
{
  // Initialize the seed and call the Mersienne algo
  std::random_device rd;
  std::mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in [[0, 1]
  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);


  // initialize SpinMatrix
  InitializeLattice();
  Energy_vector(0) = Energy;
  MagneticMoment_vector(0) = MagneticMoment;
  number_of_flips = 0;

  // Start Monte Carlo cycles
  for (int cycles = 1; cycles <= MCcycles; cycles++){
    // The sweep over the lattice, looping over all spin sites
    for(int x =0; x < NSpins; x++) {
      for (int y= 0; y < NSpins; y++){
	int ix = (int) (RandomNumberGenerator(gen)*(double)NSpins);
	int iy = (int) (RandomNumberGenerator(gen)*(double)NSpins);
	int deltaE =  2*SpinMatrix(ix,iy)*
	  (SpinMatrix(ix,PeriodicBoundary(iy,NSpins,-1))+
	   SpinMatrix(PeriodicBoundary(ix,NSpins,-1),iy) +
	   SpinMatrix(ix,PeriodicBoundary(iy,NSpins,1)) +
	   SpinMatrix(PeriodicBoundary(ix,NSpins,1),iy));
	if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) ) {
	  SpinMatrix(ix,iy) *= -1.0;  // flip one spin and accept new spin config
	  MagneticMoment += (double) 2*SpinMatrix(ix,iy);
	  Energy += (double) deltaE;
    number_of_flips+=1;

	}
  MagneticMoment_vector(cycles) = MagneticMoment;
  Energy_vector(cycles) = Energy;
  number_of_flips_vector(cycles) = number_of_flips;
      }
    }
  }
} // end of Metropolis sampling over spins
