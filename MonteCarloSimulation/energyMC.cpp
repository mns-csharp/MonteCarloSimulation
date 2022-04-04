// File: energyMC.cpp
// ----------------
// File containing functions to implement the abstract
// Energy class.
 

#include "energyMC.h"
 
// constructor
// Method: Energy  
// Usage:   n = Energy(atoms);
// --------------------------
// Used to construct the force component of any derived energy classes.

Energy::Energy(Atom **theAtoms)
{
  atoms = theAtoms;
}

// destructor
Energy::~Energy(){}
