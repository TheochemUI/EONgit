//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "Dimer.h"

using namespace helper_functions;

Dimer::Dimer(Matter const *matter, Parameters *params)
{
    parameters    = params;
    matterCenter  = new Matter(parameters);
    matterDimer   = new Matter(parameters);
    *matterCenter = *matter;
    *matterDimer  = *matter;
    nAtoms = matter->numberOfAtoms();
    
    direction.resize(nAtoms, 3);
    rotationalPlane.resize(nAtoms, 3);
    direction.setZero();
    rotationalPlane.setZero();
    totalForceCalls = 0;
}

Dimer::~Dimer()
{
    delete matterCenter;
    delete matterDimer;
}

// was estimateLowestEigenmode. rename to compute
void Dimer::compute(Matter const *matter, AtomMatrix initialDirection)
{
    long rotations = 0;
    long forceCallsCenter;
    long forceCallsDimer;
    double rotationalForce1;
    double rotationalForce2;
    double curvature, rotationalForceChange, forceDimer, rotationAngle;
    double lengthRotationalForceOld;
    double torque = 0;
    bool doneRotating = false;

    *matterCenter = *matter; 
    rotationalForceChange = forceDimer = rotationAngle = curvature = 0;
    rotationalForce1 = 0;
    rotationalForce2 = 0;
    AtomMatrix rotationalForce(nAtoms,3);
    AtomMatrix rotationalForceOld(nAtoms, 3);
    AtomMatrix rotationalPlaneOld(nAtoms, 3);
    rotationalForce.setZero();
    rotationalForceOld.setZero();
    rotationalPlaneOld.setZero();
    initialDirection.normalize();
    direction = initialDirection;

    statsAngle = 0;
    lengthRotationalForceOld = 0;
    forceCallsCenter = matterCenter->getForceCalls();
    forceCallsDimer = matterDimer->getForceCalls();

    // uses two force calls per rotation
    while(!doneRotating)
    {
        // calculate the rotational force and curvature
        curvature = calcRotationalForce(rotationalForce);  

        // determine the new rotational plane
        determineRotationalPlane(rotationalForce, rotationalForceOld, rotationalPlaneOld, &lengthRotationalForceOld);

        // calculate the torque on the dimer
        torque = rotationalForce.squaredNorm();
        assert(std::isnormal(torque));

        // convergence scheme
        if((torque > parameters->dimerTorqueMax && rotations >= parameters->dimerRotationsMax) ||
           (torque < parameters->dimerTorqueMax && torque >= parameters->dimerTorqueMin 
            && rotations >= parameters->dimerRotationsMin) ||
           (torque < parameters->dimerTorqueMin))
        {
            doneRotating = true;
        }

        // rotational force along the rotational planes normal
        rotationalForce1 = (rotationalForce.cwise()*rotationalPlane).sum();

        rotate(parameters->dimerRotationAngle);
        
        if(!doneRotating)
        {
            // rotated dimer
            curvature = calcRotationalForce(rotationalForce); // <-- bad notation

            rotationalForce2 = (rotationalForce.cwise()*rotationalPlane).sum();

            rotationalForceChange = ((rotationalForce1 - rotationalForce2) / 
                                     parameters->dimerRotationAngle);

            forceDimer = (rotationalForce1 + rotationalForce2)/2.0;

            rotationAngle = (atan(2.0*forceDimer/rotationalForceChange)/2.0 - parameters->dimerRotationAngle/2.0);

            if(rotationalForceChange < 0)
            {
                rotationAngle = rotationAngle + M_PI/2.0;
            }

            rotate(rotationAngle);
            rotationalPlaneOld = rotationalPlane; //XXX: Is this copying correctly???
            rotations++;
        }

        #ifndef NDEBUG
            printf("DIMERROT   -----   ---------  % 9.3e   ---------  % 9.3e  % 9.3e  %9ld   ---------\n",
            torque, curvature, rotationAngle*(180.0/M_PI), rotations);
        #endif
    }

    statsTorque = torque;
    statsCurvature = curvature;
    direction.normalize();
    statsAngle = acos((direction.cwise() * initialDirection).sum());
    statsAngle *= (180.0/M_PI);
    statsRotations = rotations;

    eigenvalue = curvature;

    forceCallsCenter = matterCenter->getForceCalls()-forceCallsCenter;
    forceCallsDimer = matterDimer->getForceCalls()-forceCallsDimer;

    totalForceCalls += forceCallsCenter+forceCallsDimer;

    return;
}

double Dimer::getEigenvalue()
{
    return eigenvalue;
}

AtomMatrix Dimer::getEigenvector()
{
      return direction;
}

double Dimer::calcRotationalForce(AtomMatrix &rotationalForce)
{
    double projectedForceA, projectedForceB;
    AtomMatrix posCenter(nAtoms,3);
    AtomMatrix posDimer(nAtoms,3);
    AtomMatrix forceCenter(nAtoms,3);
    AtomMatrix forceA(nAtoms,3);
    AtomMatrix forceB(nAtoms,3);

    posCenter = matterCenter->getPositions();

    // displace to get the dimer configuration A
    posDimer = posCenter + direction*parameters->dimerSeparation;

    // obtain the force for the dimer configuration
    matterDimer->setPositions(posDimer);
    forceA = matterDimer->getForces();

    // use forward difference to obtain the force for configuration B
    forceCenter = matterCenter->getForces();
    forceB = 2.0*forceCenter - forceA;

    projectedForceA = (direction.cwise() * forceA).sum();
    projectedForceB = (direction.cwise() * forceB).sum();

    // remove force component parallel to dimer
    forceA = makeOrthogonal(forceA, direction);
    forceB = makeOrthogonal(forceB, direction);

    // determine difference in force orthogonal to dimer
    rotationalForce = (forceA - forceB)/parameters->dimerSeparation;

    // curvature along the dimer
    return (projectedForceB-projectedForceA)/(2*parameters->dimerSeparation);
}

void Dimer::determineRotationalPlane(AtomMatrix rotationalForce,
                                     AtomMatrix &rotationalForceOld,
                                     AtomMatrix rotationalPlaneOld,
                                     double* lengthRotationalForceOld)
{
    double a, b, gamma = 0.0;

    a = fabs((rotationalForce.cwise() * rotationalForceOld).sum());
    b = rotationalForceOld.squaredNorm();
    if(a < 0.5*b)
    {
        gamma = (rotationalForce.cwise() * (rotationalForce - rotationalForceOld)).sum()/b;
    }
    else
        gamma = 0.0;

    // new rotational plane based on the current rotational force and the previous rotational plane force
    rotationalPlane = rotationalForce + rotationalPlaneOld * (*(lengthRotationalForceOld)) * gamma;

    // plane normal is made orthogonal to the dimer direction and normalized
    *lengthRotationalForceOld = rotationalPlane.norm();
    rotationalPlane = makeOrthogonal(rotationalPlane, direction);
    rotationalPlane.normalize();

    rotationalForceOld = rotationalForce;

    return;
}

void Dimer::rotate(double rotationAngle)
{
    double cosAngle, sinAngle;

    statsAngle += rotationAngle;

    cosAngle = cos(rotationAngle);
    sinAngle = sin(rotationAngle);

    direction = direction * cosAngle + rotationalPlane * sinAngle;
    rotationalPlane = rotationalPlane * cosAngle - direction * sinAngle;
    direction.normalize();
    rotationalPlane.normalize();

    // remove component from rotationalPlane parallel to direction
    rotationalPlane = makeOrthogonal(rotationalPlane, direction);
    rotationalPlane.normalize();

    return;
}
