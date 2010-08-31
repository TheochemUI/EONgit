/*
 *===============================================
 *  Created by Andreas Pedersen on 10/4/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#include "Potentials.h"


using namespace constants;

// which potential to use is decided at preprocessor level
Potentials::Potentials(Parameters *parameters){
    parameters_ = parameters;
//_______________________    
// To use a new potential.
// An interface should be created in the file NewPotential_interface.cpp. 
// Code will stop at runtime if used and no new potential has been defined!
    if(parameters_->getPotentialTag() == getPotentialNewPotential()){
        //interface_ = new NewPotential();
        //interface_->initialize();
		printf("The new potential must be commented in Potentials.cpp.\n");
        std::exit(1);
    }
//_______________________  
    else if(parameters_->getPotentialTag() == getPotentialLJ()){
        interface_ = new LJ();
        interface_->initialize();
    }
    else if(parameters_->getPotentialTag() == getPotentialMorse()){
        interface_ = new Morse();
        interface_->initialize();
    }
    else if(parameters_->getPotentialTag() == getPotentialEMT()){
        interface_ = new EffectiveMediumTheory();
        interface_->initialize();
    }
    else if(parameters_->getPotentialTag() == getPotentialEAM()){
        interface_ = new EAM();
        interface_->initialize();
    }
	else if(parameters_->getPotentialTag() == getPotentialQSC()){
		interface_ = new QSC();
		interface_->initialize();
    }
    else if(parameters_->getPotentialTag() == getPotentialZpIce()){
        interface_ = new ZpIce();
        interface_->initialize();
    }
    else if(parameters_->getPotentialTag() == getPotentialTip4p()){
        interface_ = new Tip4p();
        interface_->initialize();
    }
#ifndef NO_FORTRAN
#error
    else if(parameters_->getPotentialTag() == getPotentialAluminum())
    {
        interface_ = new Aluminum();
        interface_->initialize();
    }
    else if(parameters_->getPotentialTag() == getPotentialLenosky()){
        interface_ = new Lenosky();
        interface_->initialize();
    }
    else if(parameters_->getPotentialTag() == getPotentialSW()){
        interface_ = new SW();
        interface_->initialize();
    }
    else if(parameters_->getPotentialTag() == getPotentialTersoff()){
        interface_ = new Tersoff();
        interface_->initialize();
    }
    else if(parameters_->getPotentialTag() == getPotentialEDIP()){
        interface_ = new EDIP();
        interface_->initialize();
    }
#endif
	else if(parameters_->getPotentialTag() == getPotentialVASP()){
        printf("VASP potential not implemented yet. Please use different one.\n");
        std::exit(1);
        //interface_ = new vasp();
        //interface_->initialize();
    }		
	else{
		printf("Potential tag not recognized: %ld\n", parameters_->getPotentialTag());
		std::exit(1);
	}	
};

Potentials::~Potentials()
{
    interface_->cleanMemory();
};

// An alike function should be provided by the force calculator.
void Potentials::force(long nAtoms, const double *positions, const long *atomicNrs, 
                          double *forces, double *energy, const double *box){
    // The call is passed to the specified force calculator.
    interface_->force(nAtoms, positions, atomicNrs, forces, energy, box);
    
    if(parameters_->getPotentialNoTranslation()){
        double tempForceX = 0;
        double tempForceY = 0;
        double tempForceZ = 0;
        
        for(long int i=0; i<nAtoms; i++) {
            tempForceX = tempForceX+forces[ 3*i ];
            tempForceY = tempForceY+forces[3*i+1];
            tempForceZ = tempForceZ+forces[3*i+2];
        }
        tempForceX = tempForceX/nAtoms;
        tempForceY = tempForceY/nAtoms;
        tempForceZ = tempForceZ/nAtoms;
        
        for(long int i=0; i<nAtoms; i++) {
            forces[ 3*i ] = forces[ 3*i ]-tempForceX;
            forces[3*i+1] = forces[3*i+1]-tempForceY;
            forces[3*i+2] = forces[3*i+2]-tempForceZ;
        }
    }
    return;
};
