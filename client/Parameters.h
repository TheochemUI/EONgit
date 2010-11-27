//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//
//-----------------------------------------------------------------------------------
/*
 *===============================================
 *  EON Parameters.h
 *===============================================
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <ctype.h>

using namespace std;

#include "Constants.h"

/** Contains all runtime parameters and results. No functionality just bookkeeping.*/
class Parameters {
public:
    Parameters();
    ~Parameters();
    int load(string filename);
    int load(FILE *file);
/** All input parameters. If one would use values different from the default values these should be specified in the file with the name set in (Constants::PARMS_FILE_NAME).*/
    enum JobType {PROCESS_SEARCH, SADDLE_SEARCH, MINIMIZATION, PARALLEL_REPLICA, REPLICA_EXCHANGE, HESSIAN, DIMER_DR, DIMER_ROTATION, DISPLACEMENT_SAMPLING, TEST};
    JobType jobType;
    long randomSeed; // seed for random generator
    long reactantStateTag; // tag to describe to which reactant state the saddle point connects
    long potentialTag; // tag to describe which potential to use. Compare with values in Constants.cpp
    long potentialNoTranslation; // translation will be removed, handled in the potential class   
    long minimizeOnly; // only perform minimization, not saddle search
    long minimizeBox; // also minimize the box dimensions if minimize_only_ is true
    long getPrefactorsTag; // tag to describe if the prefactors should be determined. 
    double convergedRelax; // converge criterion during relaxation [eV/A]
    long maximumIterations; // max iterations for saddle point searches and minimization

    double cgCurvatureStep; // finite difference step size used in conjugate gradients
    double cgMaxMoveFullRelax; // maximum displacement vector for a step during minimization
    double qmTimeStep; // time step size used in Quickmin

    double maxDifferencePos; // The distance criterion for comparing geometries
//    double neighborCutoff; // radius used in the local atomic structure analysis
//GH: this needs to be a parameter instead of a #define in Epicenters

    long processSearchMinimizeFirst;

    bool saddleRefine; // refine saddle point
    long saddleMaxJumpAttempts; // how many times the initial displacement should try to bring the system directly to convex region. If 0 a search is started after the displacement no matter what
    long saddleMaxIterations; // max iterations for saddle point searches and minimization [GH: fix comment]
    long saddleMaxIterationsConcave; // max iterations for saddle point searches and minimization [GH: fix comment]
    long saddleLowestEigenmodeDetermination; // the algorithm to be used for lowest eigenmode determination; compare with values in Constants.cpp (now in SaddlePoint.h)
    long saddleTypePerturbation; // displacement type to use; compare with values in Constants.cpp (now in SaddlePoint.h)
    double saddleConverged; // converge criterion during saddle point search [eV/A]
    double saddleMaxStepSize; // Max length of the norm of the displacement when positive eigenvalue [A]
    double saddleMaxEnergy; // Energy above product state that will cause termination of the saddle point search [eV]
    double saddleNormPerturbation; // The norm or the perturbation vector [A]
    double saddleMaxSinglePerturbation; // max value of displacement in x, y and z direction for atoms being perturbated [A]
    double saddleWithinRadiusPerturbated; // Atoms within this radius this of the one defining the center of the displacement are also being dispalced with the value sizePerturbation_SP_ [A]
    double saddlePerpendicularForceRatio; // proportion to keep of the perpendicular force when the lowest eigenvalue is positive 

    long hessianMaxSize; // if specified, the size of the hessian determined will be equal or smaller than this value
    int hessianKind;
    double hessianMinDisplacement; // atomic displacement between min1 and the saddle point or min2 and the saddle point causing the atom to be accounted for in the Hessian [A]
    double hessianWithinRadiusDisplaced; // atoms within this radius of one the atom considered displace are also accounted for in the Hessian [A]
    double hessianPrefactorMax; // max prefactor allowed
    double hessianPrefactorMin; // min prefactor allowed

    long dimerRotationsHigh;
    long dimerRotationsLow;
    double dimerWindowHigh;
    double dimerWindowLow;
    //long dimerRotationsNewSearch; // number of iteration before starting a new saddle point search used in dimer
    double dimerSeparation; // distance between the two dimer images
    double dimerRotationAngle; // finite difference rotation angle
    double dimerMaxIterations;
    
    long   displaceNSamples;
    long   displaceIterMax;             // The maximum number of rotations to perform on the dimer.
    double displaceTorqueConvergence;   // The convergence criteria of the dimer rotation.
    double displaceMaxCurvature;        // The maximum curvature for which a sample is considered good. Used to avoid shallow but negative curvatures.
    double displaceMaxDE;               // The maximum dE for which a sample is considered good.
    string displaceCutoffs;
    string displaceMagnitudes;
    
    double lanczosConvergence; // Difference between the lowest eignevalues of two successive iterations.
    int lanczosIteration; // maximum number of iteration

    double mdTimeStep;
    double mdTemperature;
    double PRD_MaxMovedDist;
    bool   mdRefine;
    bool   mdAutoStop;
    long   RefineAccuracy;
    long   mdSteps;
    long   DephaseSteps;
    long   DephaseConstrain;
    long   DH_CheckType;
    long   CheckFreq;
    long   NewRelaxSteps;
    bool   BondBoost;
    long   BBRMDS;
    double BBDVMAX;
    double BBQRR; 
    double BBPRR; 
    double BBQcut;
 
    double ThermoType;
    double Andersen_Alpha;
    double Andersen_Tcol;
  
private:
    string toLowerCase(string s);
};
#endif
