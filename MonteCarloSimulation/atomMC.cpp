// File: atomMC.cpp
// --------------
// File containing the method bodies to implement the
// atom class.
 
#include "atomMC.h"

// Method:  setType
// Usage:   setType(int atomType);
//  ------------------------------
// Used to identify the atom as belonging to a certain
// _type_ of component within the ensemble

void Atom::setType(int t)
{
  type = t;
}

// Method: setrCutOff
// Usage:  setrCutOff(double **atomic_rCutOffValue);
// -------------------------------------------------
// Used to assign the values for the rCutOffs for the various
// atom types that make up the ensemble.

void Atom::setrCutOff(double **newrCut)
{
  rCutOff = newrCut;
}


// Method:  setMass
// Usage:   setMass(double atomicMass);
//  -----------------------------------
// Assign the mass of the atom

void Atom::setMass(double newMass)
{
  mass = newMass;
}

// Method:  setPosition
// Usage:   setPosition(double *atomicPosition);
// ---------------------------------------------
// Sets the pointer within the atom aboject to point to a
// new 1D array of doubles, representing the position
// of the atom in its various dimensions.  Assuming a 3D
// simulation, the required array would have (only) three
// elements, where:
// atomicPosition[0] = position in X dimension
// atomicPosition[1] = position in Y dimension
// atomicPosition[2] = position in Z dimension
// Note:  this method requires an _array_ to be passed to it
// which has already had memory allocated to it.  It does not
// take on the values in the array, but the address of the
// array itself.

void Atom::setPosition(double *newPos)
{
  position = newPos;
}

// Method setTrialPosition
void Atom::setTrialPosition(double *newPos)
{
  trialPosition[0] = newPos[0];
  trialPosition[1] = newPos[1];
  trialPosition[2] = newPos[2];
}

void Atom::reSetPosition()
{
  position[0] = trialPosition[0];
  position[1] = trialPosition[1];
  position[2] = trialPosition[2];
}
// Access (Get) Methods
// Usage:   n = getXxxx();
//  --------------------
// These methods return the _value_ stored for the mass or
// type of the atom.

// Method: getType
// Usage: p = getType();
// --------------------- 
// Gets the type associated with each atom.

int Atom::getType()
{
  return type;
}

// Method: getMass
// Usage: n = getMass();
// --------------------- 
// Get values of the atomic masses.

double Atom::getMass()
{
  return mass;
}



// Access (Get) Methods
// Usage:   n = getXxxxx();
// ------------------------
// These methods return a reference to the arrays which hold
// the values for rCutOff, acceleration, position, velocity,
// time step, and force, as described above their respective
// set method.

// Method: getrCutOff
// Usage: n = getrCutOff();
// ------------------------
//  returns the cut off distances for atom pairs

double **Atom::getrCutOff()
{
  return rCutOff;
}


// Method: getPosition
// Usage: n = getPosition();
// ------------------------- 
//  return position of atom

double *Atom::getPosition()
{
  return position;
}


// Method: getTrialPosition
double *Atom::getTrialPosition()
{
  return trialPosition;
}
// Constructor:Atom
// Usage: Atom atom;
// ----------------- 
// Builds the Atom class

Atom::Atom(int theType, double theMass, int dimensions)
{
  //int i;
  mass = theMass;
  type = theType;

  // allocate memory for the arrays
  position = new double[dimensions];
  trialPosition = new double[dimensions];
}

// 2nd constructor for array reference construction
Atom::Atom(){}
