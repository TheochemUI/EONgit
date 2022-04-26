
#include "Quickmin.h"
#include "HelperFunctions.h"
#include <cmath>

Quickmin::Quickmin(ObjectiveFunction *objfPassed, Parameters *parametersPassed)
{
    objf = objfPassed;
    parameters = parametersPassed;
    dt = parametersPassed->optTimeStep;
    velocity.resize(objf->degreesOfFreedom());
    velocity.setZero();
    iteration = 0;
}

Quickmin::~Quickmin()
{
    return;
}

int Quickmin::step(double maxMove)
{
    VectorXd force = -objf->getGradient();
    if (parameters->optQMSteepestDecent) {
        velocity.setZero();
    }
    else {
        if (velocity.dot(force) < 0) {
            velocity.setZero();
        }
        else {
            VectorXd f_unit = force/force.norm();
            velocity = velocity.dot(f_unit) * f_unit;
        }
    }
    
    velocity += force * dt;
    VectorXd dr = helper_functions::maxAtomMotionAppliedV(velocity * dt, parameters->optMaxMove);
    objf->setPositions(objf->getPositions() + dr);  
    iteration++;
//    return objf->isConverged();
    if(objf->isConverged()) return 1;
    return 0;
}

int Quickmin::run(int maxSteps, double maxMove)
{
    while(!objf->isConverged() && iteration < maxSteps) {
        step(maxMove);
    }
//    return objf->isConverged();
    if(objf->isConverged()) return 1;
    return 0;
}

// Previous path must be the full one, including endpoints
int Quickmin::step(const double maxMove, const std::vector<Matter> prevPath, bool& notStoppedEarly){
    int stepval = step(maxMove);
    size_t nfree = prevPath.front().numberOfFreeAtoms();
    VectorXd curPath = objf->getPositions();
    AtomMatrix curPoint = AtomMatrix::Constant(nfree, 3, 0);
    AtomMatrix prevPoint = AtomMatrix::Constant(nfree, 3, 0);
    double diff{-1};
    for (size_t idx{1}; idx < prevPath.size()-1; idx++){
        curPoint = AtomMatrix::Map(curPath.segment(3*nfree*(idx-1), 3*nfree).data(), nfree, 3);
        prevPoint = prevPath[idx].getPositionsFree();
        auto diff = (prevPoint - curPoint).norm();
        if (diff <= 5){
            notStoppedEarly = true;
        } else {
            notStoppedEarly = false;
        }
    }
    return stepval;
}
