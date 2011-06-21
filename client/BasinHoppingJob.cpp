//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <stdio.h>
#include <string>

#include "BasinHoppingJob.h"
#include "Constants.h"
#include "ConjugateGradients.h"
#include "false_boinc.h"
#include "Potential.h"
#include "HelperFunctions.h"

#ifdef BOINC
    #include <boinc/boinc_api.h>
    #include <boinc/diagnostics.h>
    #include <boinc/filesys.h>
#ifdef WIN32
    #include <boinc/boinc_win.h>
    #include <boinc/win_util.h>
#endif
#else
    #include "false_boinc.h"
#endif

using namespace std;
using namespace helper_functions;

BasinHoppingJob::BasinHoppingJob(Parameters *params)
{
    parameters = params;
    current = new Matter(parameters);
    trial = new Matter(parameters);
}

BasinHoppingJob::~BasinHoppingJob()
{
    delete current;
    delete trial;
}

std::vector<std::string> BasinHoppingJob::run(void)
{
    Matter *tmpMatter = new Matter(parameters);

    current->con2matter("reactant_passed.con");
    *trial = *current;
    *tmpMatter = *current;

    ConjugateGradients cgMin(tmpMatter, parameters);
    cgMin.setOutput(0);
    cgMin.fullRelax();

    double currentEnergy = tmpMatter->getPotentialEnergy();
    double minimumEnergy = currentEnergy;

    Matter *minimumEnergyStructure = new Matter(parameters);
    *minimumEnergyStructure = *current;

    for (int step=0; step<parameters->basinHoppingSteps; step++)
    {

        AtomMatrix displacement;
        if(parameters->basinHoppingSingleAtomDisplace)
        {
            displacement = displaceSingle();
        }
        else
        {
            displacement = displaceRandom();
        }

        trial->setPositions(current->getPositions() + displacement);

        *tmpMatter = *trial;
        ConjugateGradients cgMin(tmpMatter, parameters);
        cgMin.setOutput(0);
        cgMin.fullRelax();

        double deltaE = tmpMatter->getPotentialEnergy()-currentEnergy;
        if (abs(deltaE)>50.0) {
            printf("ERROR: huge deltaE: %f\n", deltaE);
            return returnFiles;
        }
        double p = exp(-deltaE / (parameters->temperature*8.617343e-5));

        if (randomDouble(1.0)<min(1.0, p)) 
        {
            *current = *trial;
            if(parameters->basinHoppingStayMinimized)
            {
                *current = *tmpMatter;
            }
            currentEnergy = tmpMatter->getPotentialEnergy();
            if (abs(deltaE)>parameters->structureComparisonEnergyDifference) 
            {
                *current = *tmpMatter;
                if (currentEnergy < minimumEnergy) 
                {
                    minimumEnergy = currentEnergy;
                    *minimumEnergyStructure = *current;
                }
            }
            tmpMatter->matter2xyz("movie", true);
        }
        printf("step: %6i energy: %10.4f c_energy: %10.4f de: %10.2e min_fc: %ld\n",
               step, currentEnergy, current->getPotentialEnergy(),
               deltaE, cgMin.totalForceCalls);
        boinc_fraction_done((double)(step+1)/(double)parameters->basinHoppingSteps);
    }

    /* Save Results */

    FILE *fileResults, *fileProduct;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");

    fprintf(fileResults, "%d termination_reason\n", 0);
    fprintf(fileResults, "%e minimum_energy\n", minimumEnergy);
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fclose(fileResults);

    std::string productFilename("product.con");
    returnFiles.push_back(productFilename);
    fileProduct = fopen(productFilename.c_str(), "wb");
    minimumEnergyStructure->matter2con(fileProduct);
    fclose(fileProduct);

    delete tmpMatter;
    delete minimumEnergyStructure;
  
    return returnFiles;
}

AtomMatrix BasinHoppingJob::displaceRandom()
{
    // Create a random displacement.
    AtomMatrix displacement;
    displacement.resize(trial->numberOfAtoms(), 3);
    displacement.setZero();

    for(int i = 0; i < trial->numberOfAtoms(); i++)
    {
        for(int j = 0; j < 3; j++)
        {
            if(!trial->getFixed(i))
            {
		double md = parameters->basinHoppingMaxDisplacement;
                displacement(i, j) = randomDouble(2*md) - md;
		
            }
        }
    }
    return displacement;
}

AtomMatrix BasinHoppingJob::displaceSingle()
{
    // Create a random displacement.
    AtomMatrix displacement;        
    displacement.resize(trial->numberOfAtoms(), 3);
    displacement.setZero();
    
    long ra = (long)randomDouble(trial->numberOfAtoms());
    while(trial->getFixed(ra))
    {
        ra = (long)random(trial->numberOfAtoms());
    }

    for(int j = 0; j < 3; j++)
    {
	double md = parameters->basinHoppingMaxDisplacement;
        displacement(ra,j) = randomDouble(2*md) - md;
    }

    return displacement;
}
