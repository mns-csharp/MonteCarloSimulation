// File: ljatomMC.cpp
// ------------------ 
// File containing the functions to implement the
// LJatom class.

#include "LJatom.h"

// Methods:  setSigma and setEpsilon
// Usage:    setSigma(atomicSigma);  (or)
//	     setEpsilon(atomicEpsilon);
// ---------------------------------------
// These methods accept a _reference_ to a 2D array of doubles,
// representing the sigma and epsilon values for the interactions
// between each __atomic pair type__.
// each dimesnion will be _exactly_ the number of atom types in
// the ensemble, where
// atomicSigma[0][0] = the sigma value for interactions between an
// atom of type [0] and an atom of type 0,
// atomicSigma[0][1] = the sigma value for interactions between an
// atom of type 0 and an atom of type 1, and so on, likewise
// for atomicEpsilon
	

// Method: setEpsilon
// Usage: setEpsilon(double **)
// ---------------------------- 
// Assign values of LJ epsilon parameter

void LJatom::setEpsilon(double **newEpsilon)
{
  epsilon = newEpsilon;
}

// Method: setSigma
// Usage: setSigma(double **);
// --------------------------- 
// Set the LJ parameter Sigma

void LJatom::setSigma(double **newSigma)
{
  sigma = newSigma;
}


// Methods:  getEpsilon and getSigma
// Usage:    n = getEpsilon();  (or)
//           n = getSigma();
// ---------------------------------
// These methods return references to the arrays which store the values
// for sigma and epsilon for that atom.

// Method: getEpsilon
// Usage: n = getEpsilon();
// ------------------------ 
// Get values of the Lennard-Jones epsilon parameter.

double **LJatom::getEpsilon()
{
  return epsilon;
}

// Method: getSigma
// Usage: n = getSigma();
// ---------------------- 
// Get values of the Lennard-Jones sigma parameter.

double **LJatom::getSigma()
{
  return sigma;
}


// Constructor: LJatom
// Usage: LJatom atom;
// ------------------- 
// Builds the LJatom class

LJatom::LJatom(int theType, double theMass, double **ep,
		double **sig, double **rC, int dimensions)
       :Atom(theType, theMass, dimensions)
{
  epsilon = ep;
  sigma = sig;
  rCutOff = rC;
}

LJatom::LJatom()
{
}



