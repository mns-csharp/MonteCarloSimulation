#include "MonteCarlo.h"

// File:  mc.cpp
// Class: MonteCarlo
// -----------------
// This file contains the implementation details for the Monte Carlo
// class to perform an MC simulation

void MonteCarlo::runNVE()
{
	std::ofstream out;
  out.open("mc_nve.out", std::ios::out);

  if(out.fail()) {
	  std::cout << "could not open output file!" << std::endl;
     return;
  }
  int i, index, j = 0, n, nTotal, seed = -50, accept = 0;
  int nStep  = getnStep();
  int nEquil = getnEquil();
  int nSize  = getnSize();
  int num       = ensemble->getNumAtom();
  double length = ensemble->getLength();
  double avPotEnergy = 0.0;
  double avKinEnergy = 0.0;
  double avTotEnergy = 0.0;
  double avTemperature = 0.0;
  double erKinEnergy = 0.0;
  double erPotEnergy = 0.0;
  double erTotEnergy = 0.0;
  double erTemperature = 0.0;
  double *accumPotE, *accumKinE;
  double *accumTotE, *accumTemp;
  double *trial, *position, f;
  double eOld, eTrial, trialPotE, potE, kinE, totE, temp;
  double rMax, eFixed, kNew, kOld, wFactor;
  double tally = 0.0, lrc;
  double const MAX = 0.1;
  rMax = MAX*length; 
  nTotal = nStep - nEquil;
  n = nTotal/nSize;
  f = 3.0*num/2.0 - 1.0; 
  accumPotE = new double [n];
  accumKinE = new double [n];
  accumTotE = new double [n];
  accumTemp = new double [n];
  position = new double [3];
  trial = new double [3];

  for (i = 0; i < n; i++)
  {
    accumPotE[i] = 0.0;
    accumKinE[i] = 0.0;
    accumTotE[i] = 0.0;
    accumTemp[i] = 0.0;
  }
  atom = ensemble->getAtoms();
  ensemble->lrc(); //perform long range corrections
  lrc = ensemble->getEnergyLRC();
  ensemble->setEnergy();  //set initial energy 
  eFixed = ensemble->gettotEfixed();
  
  std::cout << "Results being output to file \"mc_nve.out\"" << std::endl;

  for(i = 0; i < nStep; i++)
  {
     potE = ensemble->getPotEnergy(); 
     kOld = eFixed - (potE + lrc); 
     for (index = 0; index < num; ++index)
     {
       position  = atom[index]->getPosition();
       // find trial position 
       trial[0] =  position[0] + rMax*(2.0*random(&seed)-1.0);
       trial[1] =  position[1] + rMax*(2.0*random(&seed)-1.0);
       trial[2] =  position[2] + rMax*(2.0*random(&seed)-1.0);
       // apply periodic boundary conditions
       trial[0] -= length * nearestInt(trial[0],length);
       trial[1] -= length * nearestInt(trial[1],length);
       trial[2] -= length * nearestInt(trial[2],length);
       atom[index]->setTrialPosition(trial);
      }
      // calculate energy of trial ensemble
      trialPotE = ensemble->getTrialPotE(atom); 
      kNew = eFixed - (trialPotE +lrc);
      if (kNew > 0)
      {
        wFactor = pow(kNew/kOld, f);
        if (wFactor > 1 || wFactor > random(&seed))  //accept move
        { 
          for (index = 0; index < num; ++index)   
            atom[index]->reSetPosition();
          ensemble->setKineticE(kNew);
          ensemble->reSetEnergy(trialPotE);
          accept++;
        }

       if(i >= nEquil)
       {
         potE = ensemble->getPotEnergy() + lrc;
         kinE = ensemble->getKineticE();
         totE = potE + kinE;
         temp = 2.0*kinE/(3.0*num);
         avPotEnergy += potE;
	 avKinEnergy += kinE;
         avTotEnergy += totE;
         avTemperature += temp;

	 if((i != nEquil) && ((i % nSize) == 0))
         {
	   accumPotE[j] /= nSize;
	   accumKinE[j] /= nSize;
           accumTotE[j] /= nSize;
           accumTemp[j] /= nSize;
	   j++;
	 }

	 accumPotE[j] += potE;
	 accumKinE[j] += kinE;
         accumTotE[j] += totE;
         accumTemp[j] += temp;
        }
     }
 
     tally = accept/(double) (i + 1);
     if (tally < 0.5)
        rMax *=0.95;
     else
        rMax *= 1.05; 
  }
  avPotEnergy /= nTotal;
  avKinEnergy /= nTotal;
  avTotEnergy /= nTotal;
  avTemperature /= nTotal;

  for(i = 0; i < n; i++)
  {
    erPotEnergy += pow((accumPotE[i] - avPotEnergy),2);
    erKinEnergy += pow((accumKinE[i] - avKinEnergy),2);
    erTotEnergy += pow((accumTotE[i] - avTotEnergy),2);
    erTemperature += pow((accumTemp[i] - avTemperature),2); 
  }
//  out <<"Average Potential Energy:\t" << avPotEnergy
//		<<"\t+/-\t"<<sqrt(erPotEnergy)/nTotal<<endl;
//  out <<"Average Kinetic Energy:\t" << avKinEnergy
//		<<"\t+/-\t"<<sqrt(erKinEnergy)/nTotal<<endl;
  out <<"Average Potential Energy:\t" << avPotEnergy
       <<"  +/-  "<<sqrt(erPotEnergy)/nTotal<< std::endl;
  out <<"Average Kinetic Energy:\t\t" << avKinEnergy
       <<"  +/-  "<<sqrt(erKinEnergy)/nTotal<< std::endl;
  out <<"Average Total Energy:\t\t" << avTotEnergy
      <<"  +/-  "<<sqrt(erTotEnergy)/nTotal<< std::endl;
  out <<"Average Temperature:\t\t" << avTemperature
      <<"  +/-  "<<sqrt(erTemperature)/nTotal<< std::endl;
  out <<"Acceptance Rate:\t\t"<<100*tally <<" %"<< std::endl;
  out.close();
}

// Method: readInNVE
// Usage: readInNVE();
// ------------------- 
// ReadIn reads-in the NVE ensemble settings from the
// NVEfileMC.dat data file.

void MonteCarlo::readInNVE(int interPotential)
{
  int i, j;
  double numType, *molFract, **epsilon, **sigma, 
         **rCut,  potE, density;
  int dimensions = 3, numComp, numAtom, *comp;
  Atom *newAtom, **atoms;

  // open the NVEfileMC.dat file

  std::ifstream in;
  in.open("NVEfileMC.dat");

  if(in.fail()){
	  std::cout << "Cannot open NVEFileMC.dat!\n";
	 return;
  }


  in >> numComp;

  if(numComp <= 0){
	  std::cout << "Number of components must be > 0!\n";
	 return;
  }

  comp = new int[numComp];

  in >> numAtom;

  if(numAtom < 4){
	  std::cout << "At least 4 atoms are required!\n";
    return;
  }

  if(!(molFract = new double [numComp]))
  {
	  std::cout << "Cannot allocate memory to molFract" << std::endl;
	 return;
  }
  for(i = 0; i < numComp; i++)
      in >> molFract[i];

  double *mass = new double[numComp];

  for(i = 0; i < numComp; i++)
	 in >> mass[i];

  in >> potE;
  in >> density;
  in.close();

  // construct array of atoms
  if(interPotential == 1) // Lennard-Jones
  {
	in.open("paramLJ.dat");
	if(in.fail()){
		std::cout << "Cannot open paramLJ.dat!\n";
	  return;
	}

	// assign memory to arrays
	if(!(atoms = new Atom * [numAtom]))
	{
		std::cout << "Cannot allocate memory to atoms!\n";
	 return;
	}

	epsilon = getMemory(numComp);
	sigma = getMemory(numComp);
	rCut = getMemory(numComp);

	for (i = 0; i < numComp; i++)
	{
	  for (j = 0; j < numComp; j++)
	  {
	    in >> epsilon[i][j] >> sigma[i][j] >> rCut[i][j];

	    //adjust for indistinguishable pairs
	    if (j!=i)
	    {
	     epsilon[j][i] = epsilon[i][j];
	     sigma[j][i] = sigma[i][j];
	     rCut[j][i] = rCut[i][j];
	    }
	  }
	}

	// close input file
	in.close();

	int *atnum = new int[numComp + 1];
	atnum[0] = 0;

	// calculate the numebr of atoms of each type from the molFract
	for(i = 0; i < numComp; i++)
	{
	  numType = numAtom * molFract[i];
	  int rem = (int) numType;

	  if ((numType - rem) >= .5)
	     atnum[i+1] = (int) numType + 1;
	  else
	     atnum[i+1] = (int) numType;
	}

	int sum = 0;
	for (i = 0; i <= numComp; i++)
	sum += atnum[i];

	while (atnum[numComp] < numAtom)
	       atnum[numComp]++;

	for(i = 0; i < numComp; i++)
		atnum[i+1] += atnum[i];

	// loop through the number of components
	for (i = 0; i < numComp; i++)
	{
	  //loop through the number of atoms of that type
	  for (int j = atnum[i]; j < atnum[i+1]; j++)
	  {
	    //create atom(s) and assign mass, type, sigma, epsilon, and rCut
	    //(derivatives - 2) because acceleration and velocity are stored
	    //separately from the higher derivatives of time
	    newAtom = new LJatom(i, mass[i], epsilon, 
				 sigma, rCut, dimensions);

	    //store reference to atom in array
	    atoms[j] = newAtom;
  	  }
	}
  } // end of Lennard-Jones atom array creation

  ensemble = new NVEensemble(atoms, numComp, numAtom, potE,
                             density, molFract, comp);

  return;
}

// Method: getnSize()
// Usage: n = getSize();
// --------------------
// Return the size of the blocks for averaging
// ensemble properties.
 	
int MonteCarlo::getnSize()
{
   return nSize;
}

// Method: getnEquil()
// Usage: n = getnEquil();
// -----------------------
// Return the number of steps prior to
// equilibrium.

int MonteCarlo::getnEquil()
{
  return nEquil;
}

// Method: getnStep()
// Usage: n = getnStep();
// ----------------------
// Return the total number of simulation
// steps.

int MonteCarlo::getnStep()
{
   return nStep;
}
 
// Method:  run
// Usage:  run();
// --------------
// Runs the appropriate simularion run method 
// determined by the choice of ensemble.

void MonteCarlo::run()
{
  switch(theEnsemble)
  {
    case 1: // NVE ensemble
      runNVE();
      break;
      // other ensembles can be inserted here
  }
}

// Constructor
// -----------
// Accesses parameter file mc.dat, which identifies the choice of
// the choice of intermolecular potential, and ensemble, and
// constructs the MonteCarlo object.
 
MonteCarlo::MonteCarlo()
{
  int potential;
  std::ifstream in;

  in.open("mc.dat");

  if(in.fail())
  {
	  std::cout << "Unable to open mc.dat for MC parameters" << std::endl;
    return;
  }

  in >> nStep;
  in >> nEquil;
  in >> nSize;
  in >> theEnsemble >> potential;
  in.close();

  switch(theEnsemble)
  {
    case 1: //NVE ensemble
      readInNVE(potential);
      break;
     //other ensembles here
  }

}

