/*
 *===============================================
 *  Created by Rye Terrell on 6/7/10.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#include "QMBox.h"

using namespace helper_functions;
using namespace constants;

QMBox::QMBox(Matter *matter, Parameters *parameters)
{
    matter_ = matter;    
    parameters_ = parameters;
    nFreeCoord_ = 3 * matter->numberOfFreeAtoms();
    tempListDouble_ = new double[nFreeCoord_];
    boxforce_ = new double[3];
    boxv_ = new double[3];
    boxv_[0] = 0.0; boxv_[1] = 0.0; boxv_[2] = 0.0; 
    dR = 0.001;
    dT = 0.01;
    qmBox_ = new Quickmin(matter_, parameters);
};

QMBox::~QMBox(){

    // matter_ should not be deleted
    // parameters_ should not be deleted
    // forces_ should not be deleted
    // Are pointers to objects outside the scope
    
    delete [] tempListDouble_;
    delete [] boxforce_;
    delete [] boxv_;
    return;
};


void QMBox::increment_velocity()
{
    double bx = matter_->getBoundary(0);
    double by = matter_->getBoundary(1);
    double bz = matter_->getBoundary(2);
    double E = matter_->potentialEnergy();
    matter_->setBoundary(0, bx + dR);
    double dxE = matter_->potentialEnergy();
    matter_->setBoundary(0, bx);
    matter_->setBoundary(1, by + dR);
    double dyE = matter_->potentialEnergy();
    matter_->setBoundary(1, by);
    matter_->setBoundary(2, bz + dR);
    double dzE = matter_->potentialEnergy();
    matter_->setBoundary(2, bz);
    boxforce_[0] = (E - dxE) / dR;
    boxforce_[1] = (E - dyE) / dR;
    boxforce_[2] = (E - dzE) / dR;
    for(int i = 0; i < 3; i++)
    {
        if(boxforce_[i] * boxv_[i] < 0.0)
        {
            boxv_[i] *= 0;//boxforce_[i] * dT;
            dT *= 0.999;
        }
        else
        {
            boxv_[i] += boxforce_[i] * dT;
        }
    }
    return;
}

void QMBox::oneStep()
{
    increment_velocity();
    double bx = matter_->getBoundary(0);
    double by = matter_->getBoundary(1);
    double bz = matter_->getBoundary(2);
    //printf("%12.8f   %12.8f   %12.8f   %12.8f\n", bx, by, bz, matter_->potentialEnergy());
    matter_->setBoundary(0, bx + boxv_[0] * dT);
    matter_->setBoundary(1, by + boxv_[1] * dT);
    matter_->setBoundary(2, bz + boxv_[2] * dT);
    double scalex = (bx + boxv_[0] * dT) / bx;
    double scaley = (by + boxv_[1] * dT) / by;
    double scalez = (bz + boxv_[2] * dT) / bz;
    matter_->getFreePositions(tempListDouble_);  
    for(int i = 0; i < nFreeCoord_ / 3; i++)
    {
        tempListDouble_[i * 3 + 0] *= scalex;
        tempListDouble_[i * 3 + 1] *= scaley;
        tempListDouble_[i * 3 + 2] *= scalez;
    }
    matter_->setFreePositions(tempListDouble_);  
    qmBox_->oneStep();
    return;
};

    

void QMBox::fullRelax()
{
    bool converged = false;
    long forceCallsTemp;
    forceCallsTemp = matter_->getForceCalls();  
    while(!converged)
    {
        oneStep();
        printf("%12.8f   %12.8f   %12.8f   %12.8f   %12.8f\n", boxforce_[0], boxforce_[1], boxforce_[2], matter_->potentialEnergy(), dT);
        converged = isItConverged(parameters_->getConverged_Relax());
//        std::cout<<matter_->potentialEnergy()<<"\n";
    }
    forceCallsTemp = matter_->getForceCalls() - forceCallsTemp;
    parameters_->addForceCalls(forceCallsTemp);
    return;
};


bool QMBox::isItConverged(double convergeCriterion)
{
    double diff;
    if(!qmBox_->isItConverged(convergeCriterion))
    {
        return false;
    }
    for(int i = 0; i < 3; i++)
    {
        diff = fabs(boxforce_[i]);
        if(convergeCriterion < diff)
        {
            return false;
        }
    }
    return true;
};





















































//void QMBox::oneStepPart2(double *freeForces)
//{
//    double dotVelocityForces;
//    double dotForcesForces;
//    double *velocity;
//    velocity = new double[nFreeCoord_];

//    // Keep a copy of the pointer to the force array, is being used when it is 
//    // decided if the calculation is converged
//    forces_ = freeForces;
//    
//    matter_->getFreeVelocities(velocity); 
//    multiplyScalar(tempListDouble_, freeForces, 0.5*getTimeStep(), nFreeCoord_);
//    add(velocity, tempListDouble_, velocity, nFreeCoord_);
//    
//    dotVelocityForces = dot(velocity, freeForces, nFreeCoord_);

//    // Zeroing all velocities if they are not orthogonal to the forces
//    if(dotVelocityForces < 0)
//        multiplyScalar(velocity, velocity, 0., nFreeCoord_);
//    else{
//        dotForcesForces = dot(freeForces, freeForces, nFreeCoord_);
//        multiplyScalar(velocity, freeForces, 
//                       dotVelocityForces/dotForcesForces, nFreeCoord_);
//    }
//    matter_->setFreeVelocities(velocity);      

//    delete [] velocity;
//    return;
//};

//void QMBox::oneStepPart1(double *freeForces){
//    double *positions;
//    double *velocity;
//    positions = new double[nFreeCoord_];

//    velocity = new double[nFreeCoord_];

//    matter_->getFreeVelocities(velocity);  
//    
//    multiplyScalar(tempListDouble_, freeForces, 0.5 * getTimeStep(), nFreeCoord_);
//    add(velocity, tempListDouble_, velocity, nFreeCoord_);
//    matter_->setFreeVelocities(velocity);
//    
//    matter_->getFreePositions(positions);       
//    multiplyScalar(tempListDouble_, velocity, getTimeStep(), nFreeCoord_);
//    add(positions, tempListDouble_, positions, nFreeCoord_);
//    matter_->setFreePositions(positions);  

//    delete [] positions;
//    delete [] velocity;
//    return;
//};


