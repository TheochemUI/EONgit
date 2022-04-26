
//Based on the SteepestDescent minimizer written in ASE.

#include "SteepestDescent.h"
#include "Log.h"
#include <cassert>
#include <cmath>

SteepestDescent::SteepestDescent(ObjectiveFunction *objfPassed, Parameters *parametersPassed)
{
    objf = objfPassed;
    parameters = parametersPassed;

    iteration = 0;
}

SteepestDescent::~SteepestDescent()
{
    return;
}

int SteepestDescent::step(double maxMove)
{
    VectorXd r = objf->getPositions();
    VectorXd f = -objf->getGradient();

    VectorXd dr;
    double alpha = parameters->optSDAlpha;
    if (parameters->optSDTwoPoint == true && iteration > 0) {
        VectorXd dx = r-rPrev;
        VectorXd dg = -f+fPrev;
        alpha = dx.dot(dx)/dx.dot(dg);
        if (alpha < 0) {
            alpha = parameters->optSDAlpha;
        }
        log_file("[SD] alpha: %.4e\n", alpha);
    }

    dr = alpha*f;
    dr = helper_functions::maxAtomMotionAppliedV(dr, maxMove);

    objf->setPositions(r+dr);

    rPrev = r;
    fPrev = f;

    iteration++;

//    return objf->isConverged();
    if(objf->isConverged()) return 1;
    return 0;
}


int SteepestDescent::run(int maxSteps, double maxMove)
{
    while(!objf->isConverged() && iteration < maxSteps) {
        step(maxMove);
    }
//    return objf->isConverged();
    if(objf->isConverged()) return 1;
    return 0;
}

// Previous path must be the full one, including endpoints
int SteepestDescent::step(const double maxMove, const std::vector<Matter> prevPath, bool& notStoppedEarly){
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
