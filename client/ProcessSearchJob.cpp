
#include "ProcessSearchJob.h"
#include "BasinHoppingSaddleSearch.h"
#include "BiasedGradientSquaredDescent.h"
#include "DynamicsSaddleSearch.h"
#include "EpiCenters.h"
#include "Log.h"
#include "Optimizer.h"
#include "Potential.h"
#include "Prefactor.h"
#include <stdio.h>
#include <string>
// #include <cassert>

std::vector<std::string> ProcessSearchJob::run(void) {
  string reactantFilename("pos.con");
  string displacementFilename("displacement.con");
  string modeFilename("direction.dat");

  initial = std::make_shared<Matter>(pot, params);
  if (params->saddleMethod == "min_mode" ||
      params->saddleMethod == "basin_hopping" ||
      params->saddleMethod == "bgsd") {
    displacement = std::make_shared<Matter>(pot, params);
  } else if (params->saddleMethod == "dynamics") {
    displacement = NULL;
  }
  saddle = std::make_shared<Matter>(pot, params);
  min1 = std::make_shared<Matter>(pot, params);
  min2 = std::make_shared<Matter>(pot, params);

  if (!initial->con2matter(reactantFilename)) {
    printf("Stop\n");
    exit(1);
  }

  if (params->processSearchMinimizeFirst) {
    log("Minimizing initial structure\n");
    int fi = Potential::fcalls;
    initial->relax();
    // fCallsMin += Potential::fcalls - fi;
  }

  barriersValues[0] = barriersValues[1] = 0;
  prefactorsValues[0] = prefactorsValues[1] = 0;

  if (params->saddleMethod == "min_mode" ||
      params->saddleMethod == "basin_hopping" ||
      params->saddleMethod == "bgsd") {
    if (params->saddleDisplaceType == EpiCenters::DISP_LOAD) {
      // displacement was passed from the server
      if (!saddle->con2matter(displacementFilename)) {
        printf("Stop\n");
        exit(1);
      }
      *min1 = *min2 = *initial;
    } else {
      // displacement and mode will be made on the client
      // in SaddleSearch->initialize(...)
      *saddle = *min1 = *min2 = *initial;
    }
  } else {
    *saddle = *min1 = *min2 = *initial;
  }

  AtomMatrix mode;

  if (params->saddleMethod == "min_mode") {
    if (params->saddleDisplaceType == EpiCenters::DISP_LOAD) {
      // mode was passed from the server
      mode = helper_functions::loadMode(modeFilename, initial->numberOfAtoms());
    }
    saddleSearch = new MinModeSaddleSearch(
        saddle, mode, initial->getPotentialEnergy(), params, pot);
  } else if (params->saddleMethod == "basin_hopping") {
    saddleSearch = new BasinHoppingSaddleSearch(min1, saddle, pot, params);
  } else if (params->saddleMethod == "dynamics") {
    saddleSearch = new DynamicsSaddleSearch(saddle, params);
  } else if (params->saddleMethod == "bgsd") {
    saddleSearch = new BiasedGradientSquaredDescent(
        saddle, initial->getPotentialEnergy(), params);
  }

  int status = doProcessSearch();

  printEndState(status);
  saveData(status);

  // might have been forced to be equal if the structure passed to the client
  // when determining barrier and the prefactor

  return returnFiles;
}

int ProcessSearchJob::doProcessSearch(void) {
  Matter matterTemp(pot, params);
  long status;
  int f1;
  // f1 = Potential::fcalls;

  status = saddleSearch->run();
  // fCallsSaddle += Potential::fcalls - f1;

  if (status != MinModeSaddleSearch::STATUS_GOOD) {
    return status;
  }

  // relax from the saddle point

  AtomMatrix posSaddle = saddle->getPositions();
  AtomMatrix displacedPos;

  *min1 = *saddle;

  displacedPos = posSaddle - saddleSearch->getEigenvector() *
                                 params->processSearchMinimizationOffset;
  min1->setPositions(displacedPos);

  // Potential::fcalls = 0;
  log("\nStarting Minimization 1\n");
  bool converged = min1->relax(false, params->writeMovies, false, "min1");
  // fCallsMin += Potential::fcalls;

  if (!converged) {
    return MinModeSaddleSearch::STATUS_BAD_MINIMA;
  }

  *min2 = *saddle;
  displacedPos = posSaddle + saddleSearch->getEigenvector() *
                                 params->processSearchMinimizationOffset;
  min2->setPositions(displacedPos);

  // Potential::fcalls = 0;
  log("\nStarting Minimization 2\n");
  converged = min2->relax(false, params->writeMovies, false, "min2");
  // fCallsMin += Potential::fcalls;

  if (!converged) {
    return MinModeSaddleSearch::STATUS_BAD_MINIMA;
  }

  // if min2 corresponds to initial state, swap min1 && min2
  if (!(initial->compare(*min1)) && initial->compare(*min2)) {
    matterTemp = *min1;
    *min1 = *min2;
    *min2 = matterTemp;
  }

  if ((initial->compare(*min1)) == false) {
    log("initial != min1\n");
    return MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED;
  }

  if (initial->compare(*min2)) {
    // both minima are the initial state
    log("both minima are the initial state");
    return MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED;
  }

  // use the structure passed to the client when determining
  // the barrier and prefactor for the forward process
  if (!params->processSearchMinimizeFirst) {
    min1 = initial;
  }

  // Calculate the barriers
  barriersValues[0] = saddle->getPotentialEnergy() - min1->getPotentialEnergy();
  barriersValues[1] = saddle->getPotentialEnergy() - min2->getPotentialEnergy();

  if ((params->saddleMaxEnergy < barriersValues[0]) ||
      (params->saddleMaxEnergy < barriersValues[1])) {
    return MinModeSaddleSearch::STATUS_BAD_HIGH_BARRIER;
  }

  if (barriersValues[0] < 0.0 || barriersValues[1] < 0.0) {
    return MinModeSaddleSearch::STATUS_NEGATIVE_BARRIER;
  }

  // calculate the prefactor
  if (!params->prefactorDefaultValue) {
    // f1 = Potential::fcalls;

    int prefStatus;
    double pref1, pref2;
    // XXX: no get() calls
    prefStatus = Prefactor::getPrefactors(params.get(), min1.get(), saddle.get(), min2.get(),
                                          pref1, pref2);
    if (prefStatus == -1) {
      printf("Prefactor: bad calculation\n");
      return MinModeSaddleSearch::STATUS_FAILED_PREFACTOR;
    }
    // fCallsPrefactors += Potential::fcalls - f1;

    /* Check that the prefactors are in the correct range */
    if ((pref1 > params->prefactorMaxValue) ||
        (pref1 < params->prefactorMinValue)) {
      cout << "Bad reactant-to-saddle prefactor: " << pref1 << endl;
      return MinModeSaddleSearch::STATUS_BAD_PREFACTOR;
    }
    if ((pref2 > params->prefactorMaxValue) ||
        (pref2 < params->prefactorMinValue)) {
      cout << "Bad product-to-saddle prefactor: " << pref2 << endl;
      return MinModeSaddleSearch::STATUS_BAD_PREFACTOR;
    }
    prefactorsValues[0] = pref1;
    prefactorsValues[1] = pref2;

  } else {
    // use the default prefactor value specified
    prefactorsValues[0] = params->prefactorDefaultValue;
    prefactorsValues[1] = params->prefactorDefaultValue;
  }
  return MinModeSaddleSearch::STATUS_GOOD;
}

void ProcessSearchJob::saveData(int status) {
  FILE *fileResults, *fileReactant, *fileSaddle, *fileProduct, *fileMode;

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  fileResults = fopen(resultsFilename.c_str(), "wb");
  // XXX: min_fcalls isn't quite right it should get them from
  //      the minimizer. But right now the minimizers are in
  //      the SaddleSearch object. They will be taken out eventually.

  fprintf(fileResults, "%d termination_reason\n", status);
  fprintf(fileResults, "%ld random_seed\n", params->randomSeed);
  fprintf(fileResults, "%s potential_type\n",
          helper_functions::getPotentialName(params->potential).c_str());
  // fprintf(fileResults, "%d total_force_calls\n", Potential::fcallsTotal);
  // fprintf(fileResults, "%ld force_calls_minimization\n", fCallsMin);
  //    fprintf(fileResults, "%ld force_calls_minimization\n",
  //    SaddleSearch->forceCallsMinimization + fCallsMin);
  fprintf(fileResults, "%ld force_calls_saddle\n", fCallsSaddle);
  fprintf(fileResults, "%.12e potential_energy_saddle\n",
          saddle->getPotentialEnergy());
  fprintf(fileResults, "%.12e potential_energy_reactant\n",
          min1->getPotentialEnergy());
  fprintf(fileResults, "%.12e potential_energy_product\n",
          min2->getPotentialEnergy());
  fprintf(fileResults, "%.12e barrier_reactant_to_product\n",
          barriersValues[0]);
  fprintf(fileResults, "%.12e barrier_product_to_reactant\n",
          barriersValues[1]);
  if (params->saddleMethod == "min_mode") {
    fprintf(fileResults, "%.12e displacement_saddle_distance\n",
            displacement->perAtomNorm(*saddle));
  } else {
    fprintf(fileResults, "%.12e displacement_saddle_distance\n", 0.0);
  }
  if (params->saddleMethod == "dynamics") {
    DynamicsSaddleSearch *ds = (DynamicsSaddleSearch *)saddleSearch;
    fprintf(fileResults, "%.12e simulation_time\n",
            ds->time * params->timeUnit);
    fprintf(fileResults, "%.12e md_temperature\n",
            params->saddleDynamicsTemperature);
  }
  fprintf(fileResults, "%ld force_calls_prefactors\n", fCallsPrefactors);
  fprintf(fileResults, "%.12e prefactor_reactant_to_product\n",
          prefactorsValues[0]);
  fprintf(fileResults, "%.12e prefactor_product_to_reactant\n",
          prefactorsValues[1]);
  fclose(fileResults);

  std::string reactantFilename("reactant.con");
  returnFiles.push_back(reactantFilename);
  fileReactant = fopen(reactantFilename.c_str(), "wb");
  min1->matter2con(fileReactant);

  std::string modeFilename("mode.dat");
  returnFiles.push_back(modeFilename);
  fileMode = fopen(modeFilename.c_str(), "wb");
  helper_functions::saveMode(fileMode, saddle, saddleSearch->getEigenvector());
  fclose(fileMode);
  fclose(fileReactant);

  std::string saddleFilename("saddle.con");
  returnFiles.push_back(saddleFilename);
  fileSaddle = fopen(saddleFilename.c_str(), "wb");
  saddle->matter2con(fileSaddle);
  fclose(fileSaddle);

  std::string productFilename("product.con");
  returnFiles.push_back(productFilename);
  fileProduct = fopen(productFilename.c_str(), "wb");
  min2->matter2con(fileProduct);
  fclose(fileProduct);

  return;
}

void ProcessSearchJob::printEndState(int status) {
  log("[Saddle Search] Final status: ");

  if (status == MinModeSaddleSearch::STATUS_GOOD)
    log("Success\n");

  else if (status == MinModeSaddleSearch::STATUS_BAD_NO_CONVEX)
    log("Initial displacement unable to reach convex region\n");

  else if (status == MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY)
    log("Barrier too high\n");

  else if (status == MinModeSaddleSearch::STATUS_BAD_MAX_CONCAVE_ITERATIONS)
    log("Too many iterations in concave region\n");

  else if (status == MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS)
    log("Too many iterations\n");

  else if (status == MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED)
    log("Saddle is not connected to initial state\n");

  else if (status == MinModeSaddleSearch::STATUS_BAD_PREFACTOR)
    log("Prefactors not within window\n");

  else if (status == MinModeSaddleSearch::STATUS_FAILED_PREFACTOR)
    log("Hessian calculation failed\n");

  else if (status == MinModeSaddleSearch::STATUS_BAD_HIGH_BARRIER)
    log("Energy barrier not within window\n");

  else if (status == MinModeSaddleSearch::STATUS_BAD_MINIMA)
    log("Minimizations from saddle did not converge\n");

  else if (status == MinModeSaddleSearch::STATUS_NONNEGATIVE_ABORT)
    log("Nonnegative initial mode, aborting\n");

  else if (status == MinModeSaddleSearch::STATUS_NEGATIVE_BARRIER)
    log("Negative barrier detected\n");

  else if (status == MinModeSaddleSearch::STATUS_BAD_MD_TRAJECTORY_TOO_SHORT)
    log("No reaction found during MD trajectory\n");

  else if (status == MinModeSaddleSearch::STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE)
    log("Converged to stationary point with zero negative modes\n");

  else if (status == MinModeSaddleSearch::STATUS_BAD_NO_BARRIER)
    log("No forward barrier was found along minimized band\n");

  else if (status == MinModeSaddleSearch::STATUS_ZEROMODE_ABORT)
    log("Zero mode abort.\n");

  else if (status == MinModeSaddleSearch::STATUS_OPTIMIZER_ERROR)
    log("Optimizer error.\n");

  else
    log("Unknown status: %i!\n", status);

  return;
}
