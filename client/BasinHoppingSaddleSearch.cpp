#include "BasinHoppingSaddleSearch.h"
#include "Log.h"
#include <stdio.h>
#include "NudgedElasticBand.h"
#include "MinModeSaddleSearch.h"
#include "LowestEigenmode.h"
#include "Dimer.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"

BasinHoppingSaddleSearch::BasinHoppingSaddleSearch(Matter *reactantPassed, Matter *displacementPassed,
                                           Parameters *parametersPassed)
{
    reactant = new Matter(parameters);
    *reactant = *reactantPassed;
    parameters = parametersPassed;

    saddle = displacementPassed;

    eigenvector.resize(reactant->numberOfAtoms(), 3);
    eigenvector.setZero();
}

BasinHoppingSaddleSearch::~BasinHoppingSaddleSearch()
{
    delete reactant;
    delete product;
}

int BasinHoppingSaddleSearch::run(void)
{
                                                      //minimize "saddle"
    saddle -> relax(false, true, false, "displacementmin");
    product = new Matter(parameters);
    *product = *saddle;
                                                     //accept or reject based on boltzman exp(-de/(kB*parameters->temperature))
    double eproduct, ereactant, de;
    eproduct = product->getPotentialEnergy();
    ereactant = reactant -> getPotentialEnergy();
    de = eproduct - ereactant;
    double kB = parameters->kB;
    double Temperature = parameters->temperature;
    double arg = -de/(kB*Temperature);
    double p=exp(arg);
    double r=helper_functions::random(); 
    if(ereactant < eproduct){
	if(r>p){                                      //reject
	    return 1;
	}
    }
                                                     // NEB reactant to minimized "saddle"
    NudgedElasticBand neb(reactant, product, parameters);
    neb.neb_images[0].matter2con("neb_initial_band.con", false);
    for(size_t j{1}; j<neb.nimages; j++){
	neb.neb_images[j].matter2con("neb_initial_band", true);
    }
    neb.compute();
                                                     // pick the maximum energy image along the band
    double Emax = -1e100;
    int HighestImage = 0;
    
    for (size_t i{1}; i<neb.nimages; i++) {
	double Etest = neb.neb_images[i].getPotentialEnergy();
	printf ("i: %li Etest: %f \n", i, Etest);
	if(Etest>Emax){
	    Emax = Etest;
	    HighestImage = i;
	}
    }
                                                      //do dimer
    //Calculate initial direction
    AtomMatrix r_1 = neb.neb_images[HighestImage-1].getPositions();
    AtomMatrix r_2 = neb.neb_images[HighestImage].getPositions();
    AtomMatrix r_3 = neb.neb_images[HighestImage+1].getPositions();
    AtomMatrix direction = (r_3-r_1)/2;
    MinModeSaddleSearch dim(&neb.neb_images[HighestImage],direction.normalized(),ereactant, parameters);
    dim.run();
    *saddle = neb.neb_images[HighestImage];
    eigenvalue = dim.getEigenvalue();
    eigenvector = dim.getEigenvector();
    return 0;


}

double BasinHoppingSaddleSearch::getEigenvalue()
{
    return eigenvalue;
}

AtomMatrix BasinHoppingSaddleSearch::getEigenvector()
{
    return eigenvector;
}
