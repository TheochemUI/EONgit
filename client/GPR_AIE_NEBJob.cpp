#include "GPR_AIE_NEBJob.h"
#include "ConjugateGradients.h"
#include "Potential.h"
#include "Log.h"

#include <stdio.h>
#include <string>

GPR_AIE_NEBJob::GPR_AIE_NEBJob(Parameters *parametersPassed) : fCallsNEB{0}, fCallsGPR{0}
{
    gprfunc = std::make_unique<gpr::GaussianProcessRegression>();
    parameters=parametersPassed;
}

GPR_AIE_NEBJob::~GPR_AIE_NEBJob()
{}

std::vector<std::string> GPR_AIE_NEBJob::run(void)
{
    long status;
    int f1;

    // TODO: Deal with tsInterpolate, see SafeNudgedElasticBandJob
    string reactantFilename = helper_functions::getRelevantFile("reactant.con");
    string productFilename = helper_functions::getRelevantFile("product.con");

    auto eonp = this->parameters;
    Matter matterInit(eonp), matterFinal(eonp);
    matterInit.con2matter(reactantFilename);
    matterFinal.con2matter(productFilename);

    // Prepare data for GPR
    auto config_data = helper_functions::eon_matter_to_frozen_conf_info(&matterInit,
                                                                        this->parameters->gprPotActiveRadius);
    auto atoms_config = std::get<gpr::AtomsConfiguration>(config_data);
    // Setup the runs
    auto imgArray = helper_functions::prepInitialPath(eonp);
    auto obspath = helper_functions::prepInitialObs(imgArray);
    // Setup GPR
    auto eondat = std::make_pair(*eonp, matterInit);
    *this->gprfunc = helper_functions::initializeGPR(*this->gprfunc, atoms_config, obspath, eondat);
    this->gprfunc->setHyperparameters(obspath, atoms_config, true);
    this->gprfunc->optimize(obspath);
    // Prepare GPR potential
    GPRPotential gprpot{eonp};
    // gprpot.registerGPRObject(this->gprfunc.get());
    // Set Potentials
    LJ ljpot;
    ljpot.setParams(eonp);
    matterInit.setPotential(&ljpot);
    matterFinal.setPotential(&ljpot);

    auto nebInit = new SafeNudgedElasticBand(matterInit, matterFinal, eonp);
    nebInit->compute();
    this->fCallsGPR += 1;

    bool mustUpdate = helper_functions::maybeUpdateObs(*nebInit, obspath, *eonp);

    f1 = Potential::fcalls;

    while(mustUpdate){
        this->gprfunc->setHyperparameters(obspath, atoms_config, false);
        this->gprfunc->optimize(obspath);
        // gprpot.registerGPRObject(this->gprfunc.get());
        auto nebTwo = new SafeNudgedElasticBand(matterInit, matterFinal, eonp);
        nebTwo->compute();
        this->fCallsGPR += 1;
        mustUpdate = helper_functions::maybeUpdateObs(*nebTwo, obspath, *eonp);
        delete nebTwo;
    };
    // If there is only one round and no updates, do not redo NEB
    if (this->fCallsGPR > 1){
        // Final round
        auto nebFin = new SafeNudgedElasticBand(matterInit, matterFinal, eonp);
        status = nebFin->compute();
        this->fCallsGPR += 1;

        fCallsNEB += Potential::fcalls - f1;

        if (status == SafeNudgedElasticBand::STATUS_INIT) {
            status = SafeNudgedElasticBand::STATUS_GOOD;
        }

        printEndState(status);
        saveData(status, nebFin);
        delete nebFin;
        delete nebInit;
        return returnFiles;
    } else {
        // TODO: Be more elegant when one round has occurred
        fCallsNEB += Potential::fcalls - f1;

        if (status == SafeNudgedElasticBand::STATUS_INIT) {
            status = SafeNudgedElasticBand::STATUS_GOOD;
        }

        printEndState(status);
        saveData(status, nebInit);
        delete nebInit;

        return returnFiles;
    }
}

void GPR_AIE_NEBJob::saveData(int status, SafeNudgedElasticBand *neb)
{
    FILE *fileResults, *fileNEB;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");

    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcalls);
    fprintf(fileResults, "%d force_calls_neb\n", fCallsNEB);
    fprintf(fileResults, "%d gpr_created\n", fCallsGPR);
    fprintf(fileResults, "%f energy_reference\n", neb->imageArray.front().getPotentialEnergy());
    fprintf(fileResults, "%li number_of_images\n", neb->nimages);
    for(long i=0; i<=neb->nimages+1; i++) {
        fprintf(fileResults, "%f image%li_energy\n", neb->imageArray[i].getPotentialEnergy()-neb->imageArray.front().getPotentialEnergy(), i);
        fprintf(fileResults, "%f image%li_force\n", neb->imageArray[i].getForces().norm(), i);
        fprintf(fileResults, "%f image%li_projected_force\n", neb->projectedForceArray[i].norm(), i);
    }
    fprintf(fileResults, "%li number_of_extrema\n", neb->numExtrema);
    for(long i=0; i<neb->numExtrema; i++) {
        fprintf(fileResults, "%f extremum%li_position\n", neb->extremumPositions[i], i);
        fprintf(fileResults, "%f extremum%li_energy\n", neb->extremumEnergies[i], i);
    }

    fclose(fileResults);

    std::string nebFilename("neb.con");
    returnFiles.push_back(nebFilename);
    fileNEB = fopen(nebFilename.c_str(), "wb");
    for(long i=0; i<=neb->nimages+1; i++) {
        neb->imageArray[i].matter2con(fileNEB);
    }
    fclose(fileNEB);

    returnFiles.push_back("neb.dat");
    neb->printImageData(true);
}

void GPR_AIE_NEBJob::printEndState(int status)
{
    log("\nFinal state: ");
    if(status == SafeNudgedElasticBand::STATUS_GOOD)
        log("GPR-AIE Nudged elastic band, successful.\n");
    else if(status == SafeNudgedElasticBand::STATUS_BAD_MAX_ITERATIONS)
        log("Nudged elastic band, too many iterations.\n");
    else
        log("Unknown status: %i!\n", status);
    return;
}
