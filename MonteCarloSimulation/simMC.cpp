// File:  simMC.cpp
// -----------------
// This file contains the implementation of the Simulation class,
// which defines all methods and data for a Molecular Simulation,
// using approaches such as Molecular Dynamics, Monte Carlo etc.

#include "simMC.h"

// Method:  readSimParameters
// Usage:   readSimParameters();
// -----------------------------
// Reads in parameters for NVE ensemble

void Simulation::readSimParameters()
{
  int simulation;
  // open the simMC.dat file

  std::ifstream in;
  in.open("simMC.dat");

  if(in.fail()){
	  std::cout << "Cannot open simMC.dat!\n";
	 return;
  }

  in >> simulation;
  in.close();

  switch(simulation)
  {
    case 2:  //Monte Carlo
        mc = new MonteCarlo();
        mc->run(); 		
	break;
    default:
		std::cout << "Invalid simulation selected!  Aborting" << std::endl;
	break;
  }
  return;
}

