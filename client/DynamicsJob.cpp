#include <stdio.h>
#include <string>

#include "Dynamics.h"
#include "DynamicsJob.h"
#include "Parameters.h"
#include "Potential.h"
#include "HelperFunctions.h"

using namespace std;
using namespace helper_functions;

DynamicsJob::DynamicsJob(Parameters *params)
{
    parameters = params;
}

DynamicsJob::~DynamicsJob() {}

std::vector<std::string> DynamicsJob::run(void)
{
    Matter *R = new Matter(parameters);
    Matter *F = new Matter(parameters);
    R->con2matter("pos.con");
    *F = *R;

    Dynamics *d = new Dynamics(R, parameters);
    d->run();

    *F = *R;
    std::string productFilename("final.con");
    returnFiles.push_back(productFilename);
    F->matter2con(productFilename);

    delete R;
    delete F;
    delete d;

    std::vector<std::string> returnFiles;
    return returnFiles;
}

