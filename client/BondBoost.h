//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef BONDBOOST_H
#define BONDBOOST_H

#include "Matter.h"
#include "HelperFunctions.h"
#include "Parameters.h"

#include "Eigen.h"

/** Functionality relying on the conjugate gradients algorithm. The object is capable of minimizing an Matter object or modified forces being passed in.*/
class BondBoost {

public:
    /** Constructor to be used when a structure is minimized.
    @param[in]   *matter        Pointer to the Matter object to be relaxed.
    @param[in]   *parameters    Pointer to the Parameter object containing the runtime parameters.*/
    BondBoost(Matter *matt, Parameters *params);
    ~BondBoost();///< Destructor.

    void initial();
    double boost();  

private:
    Matrix<double, Eigen::Dynamic, 1> Rmdsteps();
    long BondSelect();
    double Booststeps();
    long nAtoms;///< Number of free coordinates.
    Matter *matter;///< Pointer to atom object \b outside the scope of the class.    
    Parameters *parameters;///< Pointer to a structure outside the scope of the class containing runtime parameters. 
    long  *BAList;
    long  *RAList;
    long  *TABAList;
    long  *BBAList;
    double  *Epsr_Q;
    Matrix<double, Eigen::Dynamic, 1> TABLList;//EquilibriumTaggedAtomInvolvedBondLengthList;
    Matrix<double, Eigen::Dynamic, 1> EBBLList;//EquilibriumBoostBondLengthList
    Matrix<double, Eigen::Dynamic, 1> CBBLList;//CurrentBoostBondLengthList
    long  nBAs;
    long  nRAs;
    long  nTABs;
    long  nReg;
    long  nBBs;
    double SDtime;
    double SPtime;
    double SDtime_B;
};

class Hyperdynamics {
        public:
            enum{
                NONE,
                BOND_BOOST
            };
};

#endif
