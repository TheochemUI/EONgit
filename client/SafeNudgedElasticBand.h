#ifndef SAFENUDGEDELASTICBAND_H_
#define SAFENUDGEDELASTICBAND_H_

#include <math.h>
#include <cmath>

#include "Eigen.h"

#include "Matter.h"
#include "HelperFunctions.h"
#include "Parameters.h"

class Matter;
class Parameters;

// NEB method for determining a minimum energy path between two matter objects
class SafeNudgedElasticBand {

public:

    enum{
        STATUS_GOOD, //0
        STATUS_INIT, //1
        STATUS_BAD_MAX_ITERATIONS, //2
    };

    SafeNudgedElasticBand(Matter initialPassed, Matter finalPassed, Parameters* parametersPassed);
    ~SafeNudgedElasticBand();

    void clean(void);
    int compute(void);
    void updateForces(void);
    double convergenceForce(void);
    void findExtrema(void);
    void printImageData(bool writeToFile=false);

    size_t natoms, nimages, totImages;
    size_t maxEnergyImage, climbingImage, numExtrema;
    std::vector<Matter> imageArray;
    std::vector<AtomMatrix> tangentArray;
    std::vector<AtomMatrix> projectedForceArray;
    bool movedAfterForceCall;
    std::vector<double> extremumEnergies;
    std::vector<double> extremumPositions;
    std::vector<double> extremumCurvatures;

private:

    Parameters *parameters;

};

#endif // SAFENUDGEDELASTICBAND_H_
