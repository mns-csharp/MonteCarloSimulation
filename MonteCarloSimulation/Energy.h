// File:  energyMC.h
// -----------------
// This file contains the definition of abstract class
// Energy, which defines all methods and data common to
// different force types (eg derived classes), as well as
// foreshadowing methods of those derived classes through
// the use of abstract methods.


#ifndef _energyMC_h
#define _energyMC_h

#include <iostream>
#include <math.h>
#include "Atom.h"

class Energy 
{
  protected:
   Atom **atoms; 	// reference to the array of atoms
  public:
   Energy(Atom **);      // constructor for Energy 
   virtual ~Energy();   // destructor for Energy 
   // Abstract Methods
   // ----------------
   //abstract methods for subclasses
   // for LJforce
   virtual void setEnergy(int, double, double *)=0;
   virtual double getTrialPotE(int, double, Atom **)=0;
   virtual void lrc(int, int *, double,
		double *) = 0;  // set long range corrections
};

#endif

