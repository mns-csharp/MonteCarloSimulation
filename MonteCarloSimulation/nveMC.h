//  File:  nveMC.h
//  -----------------
//  This file contains the definition of the NVEensemble class,
//  which defines all methods and data for the microcanonical ensemble.

#ifndef _nveMC_h
#define _nveMC_h

#include "ljatomMC.h"
#include "auxfuncMC.h"
#include "ljenergyMC.h"
#include "ensembleMC.h"

class NVEensemble : public Ensemble
{
  private:
    double energyLRC;    // long range energy correction
  public:
    NVEensemble(Atom **, int, int, double,
	 double, double *, int *);// ensemble constructor
    void readInNVE();                   // read in ensemble data
    virtual void initialCoord();        //place atoms on lattice

    virtual double getTrialPotE(Atom **);
    virtual void reSetEnergy(double);
    virtual void setEnergy();           // Order N*N energy calculation
    virtual double getEnergyLRC();      // get energy long range correction
    virtual void lrc();                 // trigger the long range corrections
};

#endif
