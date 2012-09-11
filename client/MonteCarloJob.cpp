//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "MonteCarloJob.h"
#include "MonteCarlo.h"
#include "Log.h"
#include "Matter.h"
#include "HelperFunctions.h"

MonteCarloJob::MonteCarloJob(Parameters *params)
{
    parameters = params;
}

MonteCarloJob::~MonteCarloJob(){ }

std::vector<std::string> MonteCarloJob::run(void)
{
    string posInFilename("pos.con");
    string posOutFilename("out.con");

    if (parameters->checkpoint) {
        FILE *pos;
        pos = fopen("pos_cp.con", "r");
        if (pos != NULL) {
            posInFilename = "pos_cp.con";
            log("Resuming from checkpoint\n");
        }else{
            log("No checkpoint files found\n");
        }
    }

    std::vector<std::string> returnFiles;
    returnFiles.push_back(posOutFilename);

    Matter *matter = new Matter(parameters);
    matter->con2matter(posInFilename);

    //code will go
    MonteCarlo mc = MonteCarlo(matter, parameters);
    mc.run(1, 300.0, 0.2);


    FILE *fileResults;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");

    //fprintf(fileResults, "%d termination_reason\n", status);
    //fprintf(fileResults, "minimization job_type\n");
    //fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    //fprintf(fileResults, "%d total_force_calls\n", Potential::fcallsTotal);
    //if (status != STATUS_POTENTIAL_FAILED) {
    //    fprintf(fileResults, "%f potential_energy\n", pos->getPotentialEnergy());
    //}
    //fclose(fileResults);

    return returnFiles;
}
