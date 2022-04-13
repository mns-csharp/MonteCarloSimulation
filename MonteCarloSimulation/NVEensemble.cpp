// File: nveMC.cpp
// ---------------- 
// File containing functions to implement the
// NVEensemble class.
 
#include <iostream>
#include <math.h>
#include <string.h>
#include "NVEensemble.h"

// Method: getEnergyLRC
// Usage: n = getEnergyLRC();
// -------------------------- 
// Gets the long range energy correction.

double NVEensemble::getEnergyLRC()
{
  return energyLRC;
}

// Method: initialCoord
// Usage: initialCoord();
// ---------------------- 
// Places atoms upon a lattice to get initial coordinates

void NVEensemble::initialCoord()
{
  int i, j, x, y, z, offset;
  double cells, dcells, cellL, halfCellL, *tempPos, *position;
 
  // Determine the number of unit cells in each coordinate
  // direction
  dcells = pow(0.25 * (double)numAtom, 1.0/3.0);
  cells =  (int) nearestInt(dcells, 1.0);
 
  // check if numAtom is an non-fcc number of molecules
  // and increase the number of cells if necessary
 
  while((4 * cells * cells * cells) < numAtom)
                 cells = cells + 1;
 
  // Determine length of the unit cell
 
  cellL = boxLen/ (double) cells;
  halfCellL = 0.5 * cellL;

  // Construct the unit cell
  // point to atoms position
  position = atoms[0]->getPosition();
  position[0] = 0.0;
  position[1] = 0.0;
  position[2]= 0.0;

  position = atoms[1]->getPosition();
  position[0] = halfCellL;
  position[1] = halfCellL;
  position[2]= 0.0;

  position = atoms[2]->getPosition();
  position[0] = 0.0;
  position[1] = halfCellL;
  position[2]= halfCellL;

  position = atoms[3]->getPosition();
  position[0] = halfCellL;
  position[1] = 0.0;
  position[2]= halfCellL;
 
  for(i = 4; i < numAtom; i++)
  {
    position = atoms[i]->getPosition();
    position[0] = 0.0;
    position[1] = 0.0;
    position[2] = 0.0;
  } //init all other atoms to 0
 
  // Build the lattice from the unit cell by
  // repeatly translating the four vectors of
  // the unit cell through a distance cellL in
  // the x, y and z directions
 
  offset = 0;
 
  for (z = 1; z <= cells; z++)
    for (y = 1; y <= cells; y++)
      for (x = 1; x <= cells; x++){
        for (i = 0; i < 4; i++){
          j = i + offset;
          if(j < numAtom){
	    tempPos = atoms[j]->getPosition();
	    position = atoms[i]->getPosition();
            tempPos[0] = position[0] + cellL * (x-1);
            tempPos[1] = position[1] + cellL * (y-1);
            tempPos[2] = position[2] + cellL * (z-1);
          }
        }
 
      offset = offset + 4;
  }

  // Shift centre of box to the origin.

  for (i = 0; i < numAtom; i++){
      tempPos = atoms[i]->getPosition();
	 tempPos[0] -= halfCellL;
	 tempPos[1] -= halfCellL;
	 tempPos[2] -= halfCellL;
  }
}

// Method: reSetEnergy
// Usage:  reSetEnergy(energy);
// ----------------------------
// Re-sets the energy of the ensemble following
// an acceted move
  
void NVEensemble::reSetEnergy(double eng)
{
 potEnergy = eng;
}

// Method: setEnergy
// Usage: setEnergy();
// -------------------
// Set the energy calculations using LJenergy.
// This is an N* N calculation to set the energy
// at the start of the simulation.

void NVEensemble::setEnergy()
{
  theEnergy->setEnergy(numAtom, boxLen,
             &potEnergy);
}

// Method: getTrialPotE
// Usage: n = getTrialPotE(atoms)
// ---------------------------------
// Performs N*N calculations to update energy after
// all atoms are moved.
  
double NVEensemble::getTrialPotE(Atom **a)
{
  return theEnergy->getTrialPotE(numAtom, boxLen, a);
}  
// Method: lrc
// Usage: lrc();
//  ------------
// initiate the force calculations for long range corrections

void NVEensemble::lrc()
{
  theEnergy->lrc(numComp, comp, boxVol, &energyLRC);
}

// Method: NVEensemble
// Usage: NVEensemble;
// ------------------- 
//  instantiate the NVEensemble class

NVEensemble::NVEensemble(Atom **atoms, int numComp, int numAtom,
   double totEfixed, double density, 
   double *molFract, int *comp)
   : Ensemble(atoms, numComp, numAtom, totEfixed,
     density, molFract, comp)
{
  setVolume();
  setLength();
  setComp();
  initialCoord();
  theEnergy = new LJenergy(atoms);
}
