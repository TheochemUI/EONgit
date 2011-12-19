//----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <errno.h>
#include <time.h>
#include "Parameters.h"
#include "INIFile.h"
#include "EpiCenters.h"
#include "HelperFunctions.h"
#include "Hessian.h"
#include "Job.h"
#include "Dynamics.h"
#include "BondBoost.h"
#include "SaddleSearch.h"
#include "ImprovedDimer.h"
#include "NudgedElasticBand.h"
#include "Potential.h"

Parameters::Parameters(){

    // [Main] //
    job = Job::PROCESS_SEARCH;
    randomSeed = -1;
    temperature = 300.0;
    checkpoint = false;
    quiet = false;
    iniFilename = "config.ini";
    conFilename = "reactant.con";
    finiteDifference = 0.01;

    // [Prefactor]
    prefactorDefaultValue = 0.0;
    prefactorMaxValue = 1e+21;
    prefactorMinValue = 1e+9;
    prefactorWithinRadius = 5.0;
    prefactorMinDisplacement = 0.25;

    // [Potential]
    potential = Potential::POT_LJ;
    MPIPollPeriod = 0.25;
    MPIPotentialRank = -1;
    LAMMPSLogging = true;

    // [Structure Comparison] //
    distanceDifference = 0.1;
    neighborCutoff = 3.3;
    checkRotation = false;

    // [Debug] //
    writeMovies = false;
    writeMoviesSteps = 1;

    // [Process Search] //
    processSearchMinimizeFirst = false;
    processSearchMinimizationOffset = 0.2;

    // [Saddle Search] //
    saddleDisplaceType = EpiCenters::DISP_LOAD;
    saddleMinmodeMethod = LowestEigenmode::MINMODE_DIMER;
    saddleMaxStepSize = 0.2;
    saddleMaxEnergy = 20.0;
    saddleMaxIterations = 1000;
    saddleDisplaceRadius = 4.0;
    saddleDisplaceMagnitude = 0.1;
    saddleMaxSingleDisplace = 10.;
    saddleConvergedForce = 0.005;
    saddleNonnegativeDisplacementAbort = false;    
    saddlePerpForceRatio = 0.0; // undocumented
    saddleConfinePositive = false; // undocumented
    saddleConfinePositiveMinForce = 0.5; // undocumented
    saddleConfinePositiveScaleRatio = 0.9; // undocumented
    saddleConfinePositiveBoost = 10.; // undocumented
    saddleConfinePositiveMinActive = 30; // undocumented

    // [Optimizers] //
    optMethod = "cg";
    optMaxIterations = 1000;
    optConvergedForce = 0.005;
    optMaxMove = 0.2;
    optTimeStep = 0.1;
    optVariableTimeStep = false;
    optLBFGSMemory = 50;
    optLBFGSInverseCurvature = 0.1;

    // [Dimer] //
    dimerRotationAngle = 0.005;
    dimerImproved = true;
    dimerConvergedAngle = 5.0; // degrees
    dimerOptMethod = ImprovedDimer::OPT_CG;
    dimerTorqueMin = 0.1; // old dimer
    dimerTorqueMax = 1.0; // old dimer
    dimerRotationsMin = 1; // old dimer
    dimerRotationsMax = 10; // old dimer and new dimer

    // [Lanczos] //
    lanczosTolerance = 0.01;
    lanczosMaxIterations = 20;

    // [Hessian] //
    hessianAtomList = string("All");
    hessianZeroFreqValue = 1e-6;

    // [Nudged Elastic Band] //
    nebImages = 5;
    nebSpring = 5.0;
    nebClimbingImageMethod = true;
    nebOldTangent = false;
    nebMaxIterations = 1000;

    // [Dynamics] //
    mdTimeStep = 1;
    mdSteps = 1000;

    // [Thermostat] //
    thermostat = Dynamics::NONE;
    thermoAndersenAlpha = 0.2; // collision strength
    thermoAndersenTcol = 10.0; // collision frequency in unit of fs
    thermoNoseMass = 1.0;
    thermoLangvinFriction = 0.005;

    // [Parallel Replica] //
    paraRepRefine = true;
    paraRepAutoStop = false;
    paraRepRefineAccuracy = 1;
    paraRepCheckPeriod = 500;
    paraRepRecordPeriod = 50;
    paraRepRelaxSteps = 500;
    paraRepDephaseSteps = 200;
    paraRepDephaseLoopStop = false;
    paraRepDephaseLoopMax = 5;

    // [Distributed Replica] //
    drBalanceSteps = 500;
    drSamplingSteps = 500;
    drTargetTemperature = 300.0;

    // [Hyperdynamics] //
    biasPotential = Hyperdynamics::NONE;
    bondBoostBALS = string("ALL"); // boosted atom list string 
    bondBoostDVMAX = 0.0;
    bondBoostQRR = 0.2; // can not be set to 0
    bondBoostPRR = 0.95;
    bondBoostQcut = 3.0;
    bondBoostRMDS = 100;

    // [Basin Hopping] //
    basinHoppingMaxDisplacement = 0.5;
    basinHoppingSteps = 10000;
    basinHoppingQuenchingSteps = 0;
    basinHoppingSingleAtomDisplace = false;
    basinHoppingSignificantStructure = false;
    basinHoppingMaxDisplacementAlgorithm = "standard";
    basinHoppingDisplacementDistribution = "uniform";
    basinHoppingSwapProbability = 0.0;
    basinHoppingJumpMax = 10;
    basinHoppingJumpSteps = 0;
    basinHoppingInitialMD = false;
    basinHoppingInitialMDTemperature = 300.0;

}

Parameters::~Parameters(){
    return;
}

string Parameters::toLowerCase(string s)
{
    for (string::size_type i = 0; i < s.length(); ++i) {
      s[i] = tolower(s[i]);
    }
    return s;
}

int Parameters::load(string filename)
{
    FILE *fh;

    fh = fopen(filename.c_str(), "rb");
    if (fh == NULL) {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return 1;
    }
    int error = load(fh);
    fclose(fh);
    return error;
}

int Parameters::load(FILE *file){

    CIniFile ini;
    ini.CaseInsensitive();
    int error=0;

    if(ini.ReadFile(file))
    {

        // [Main] //

        job = toLowerCase(ini.GetValue("Main", "job"));
        temperature = ini.GetValueF("Main", "temperature", temperature);
        randomSeed = ini.GetValueL("Main", "random_seed", randomSeed);
        checkpoint = ini.GetValueB("Main", "checkpoint", checkpoint);
        quiet = ini.GetValueB("Main", "quiet", quiet);
        finiteDifference = ini.GetValueF("Main", "finite_difference", finiteDifference);
        // Initialize random generator
        if(randomSeed < 0){
            unsigned i = time(NULL);
            randomSeed = i;
            helper_functions::random(i);
        }else{
            helper_functions::random(randomSeed);
        }

        // [Potential] //

        potential = toLowerCase(ini.GetValue("Potential", "potential"));
        MPIPollPeriod = ini.GetValueF("Potential", "mpi_poll_period", MPIPollPeriod);
        LAMMPSLogging = ini.GetValueB("Potential", "lammps_logging", LAMMPSLogging);

        // [Debug] //

        writeMovies = ini.GetValueB("Debug", "write_movies", writeMovies);
        writeMoviesSteps = ini.GetValueL("Debug","write_movies_steps",writeMoviesSteps);

        // [Structure Comparison] //

        distanceDifference = ini.GetValueF("Structure Comparison", "distance_difference", distanceDifference);
        neighborCutoff = ini.GetValueF("Structure Comparison", "neighbor_cutoff", neighborCutoff);
        checkRotation = ini.GetValueB("Structure Comparison", "check_rotation", checkRotation);

        // [Process Search] //

        processSearchMinimizeFirst = ini.GetValueB("Process Search", "minimize_first", processSearchMinimizeFirst);
        processSearchMinimizationOffset = ini.GetValueF("Process Search", "minimization_offset", processSearchMinimizationOffset);

        // [Optimizers] //

        optMethod = toLowerCase(ini.GetValue("Optimizer", "opt_method", optMethod));
        optConvergedForce = ini.GetValueF("Optimizer", "converged_force", optConvergedForce);
        optMaxIterations = ini.GetValueL("Optimizer", "max_iterations", optMaxIterations);
        optMaxMove = ini.GetValueF("Optimizer","max_move", optMaxMove);
        optTimeStep = ini.GetValueF("Optimizer","time_step", optTimeStep);
        optVariableTimeStep = ini.GetValueB("Optimizer","variable_time_step", optVariableTimeStep);
        optLBFGSMemory = ini.GetValueL("Optimizer", "lbfgs_memory", optLBFGSMemory);
        optLBFGSInverseCurvature = ini.GetValueF("Optimizer", "lbfgs_inverse_curvature", optLBFGSInverseCurvature);

        // [Saddle Search] //

        saddleMinmodeMethod = toLowerCase(ini.GetValue("Saddle Search", "min_mode_method", saddleMinmodeMethod));
        saddleDisplaceMagnitude = ini.GetValueF("Saddle Search", "displace_magnitude", saddleDisplaceMagnitude);
        saddleDisplaceRadius = ini.GetValueF("Saddle Search", "displace_radius", saddleDisplaceRadius);
        saddleMaxEnergy = ini.GetValueF("Saddle Search", "max_energy", saddleMaxEnergy);
        saddleMaxStepSize = ini.GetValueF("Saddle Search", "max_step_size", optMaxMove);
        saddleMaxIterations = ini.GetValueL("Saddle Search", "max_iterations", optMaxIterations);
        saddleNonnegativeDisplacementAbort = ini.GetValueB("Saddle Search", "nonnegative_displacement_abort", saddleNonnegativeDisplacementAbort); 
        saddleMaxSingleDisplace = ini.GetValueF("Saddle Search", "max_single_displace", saddleMaxSingleDisplace);
        saddleConvergedForce = ini.GetValueF("Saddle Search", "converged_force", optConvergedForce);
        saddlePerpForceRatio = ini.GetValueF("Saddle Search", "perp_force_ratio", saddlePerpForceRatio); // undocumented
        saddleDisplaceType = toLowerCase(ini.GetValue("Saddle Search", "client_displace_type", EpiCenters::DISP_LOAD));
        // XXX: This is a result of mixing our server/client config files.
        if(saddleDisplaceType != EpiCenters::DISP_NOT_FCC_OR_HCP &&
           saddleDisplaceType != EpiCenters::DISP_MIN_COORDINATED &&
           saddleDisplaceType != EpiCenters::DISP_LAST_ATOM &&
           saddleDisplaceType != EpiCenters::DISP_RANDOM){
              saddleDisplaceType = EpiCenters::DISP_LOAD;
           }
        saddleConfinePositive = ini.GetValueB("Saddle Search", "confine_positive", saddleConfinePositive);
        if(saddleConfinePositive) {
            saddleConfinePositiveMinForce = ini.GetValueF("Saddle Search", "confine_positive_min_move", saddleConfinePositiveMinForce);
            saddleConfinePositiveScaleRatio = ini.GetValueF("Saddle Search", "confine_positive_scale_ratio", saddleConfinePositiveScaleRatio);
            saddleConfinePositiveBoost = ini.GetValueF("Saddle Search", "confine_positive_scale_boost", saddleConfinePositiveBoost);
            saddleConfinePositiveMinActive = ini.GetValueL("Saddle Search", "confine_positive_scale_min_active", saddleConfinePositiveMinActive);            
        }

        // [Dimer] //

        dimerRotationAngle = ini.GetValueF("Dimer", "finite_angle", dimerRotationAngle);
        dimerImproved = ini.GetValueB("Dimer", "improved", dimerImproved);
        dimerConvergedAngle = ini.GetValueF("Dimer", "converged_angle", dimerConvergedAngle);
        dimerMaxIterations = ini.GetValueL("Dimer", "max_iterations", dimerMaxIterations);
        dimerOptMethod = toLowerCase(ini.GetValue("Dimer", "opt_method", dimerOptMethod));
        dimerRotationsMin = ini.GetValueL("Dimer", "rotations_min", dimerRotationsMin); // old
        dimerRotationsMax = ini.GetValueL("Dimer", "rotations_max", dimerRotationsMax); // old & new
        dimerTorqueMin = ini.GetValueF("Dimer", "torque_min", dimerTorqueMin); // old
        dimerTorqueMax = ini.GetValueF("Dimer", "torque_max", dimerTorqueMax); // old

        // [Lanczos] //

        lanczosTolerance = ini.GetValueF("Lanczos", "tolerance", lanczosTolerance);
        lanczosMaxIterations = ini.GetValueL("Lanczos", "max_iterations", lanczosMaxIterations);

        // [Prefactor] //

        prefactorDefaultValue = ini.GetValueF("Prefactor", "default_value", prefactorDefaultValue);
        prefactorMaxValue = ini.GetValueF("Prefactor", "max_value", prefactorMaxValue);
        prefactorMinValue = ini.GetValueF("Prefactor", "min_value", prefactorMinValue);
        prefactorWithinRadius = ini.GetValueF("Prefactor", "within_radius", prefactorWithinRadius);
        prefactorMinDisplacement = ini.GetValueF("Prefactor", "min_displacement", prefactorMinDisplacement);

        // [Hessian] //

        hessianAtomList = toLowerCase(ini.GetValue("Hessian", "atom_list", hessianAtomList));
        hessianZeroFreqValue = ini.GetValueF("Hessian", "zero_freq_value", hessianZeroFreqValue);

        // [Nudged Elastic Band] //

        nebImages = ini.GetValueL("Nudged Elastic Band", "images", nebImages);
        nebSpring = ini.GetValueF("Nudged Elastic Band", "spring", nebSpring);
        nebClimbingImageMethod = ini.GetValueB("Nudged Elastic Band", "climbing_image_method", nebClimbingImageMethod);
        nebOldTangent = ini.GetValueB("Nudged Elastic Band", "old_tangent", nebOldTangent);
        nebMaxIterations = ini.GetValueL("Nudged Elastic Band", "max_iterations", optMaxIterations);

        // [Dynamics] //

        mdTimeStep = ini.GetValueF("Dynamics", "time_step", mdTimeStep);
        mdTimeStep = mdTimeStep * 0.09823; //transfer the time unit from fs to 10.18 fs
        mdSteps = ini.GetValueL("Dynamics", "steps", mdSteps);
        thermostat = toLowerCase(ini.GetValue("Dynamics", "thermostat", "andersen"));
        thermoAndersenAlpha = ini.GetValueF("Dynamics","andersen_alpha",thermoAndersenAlpha);
        thermoAndersenTcol = ini.GetValueF("Dynamics","andersen_collision_period",thermoAndersenTcol);
        thermoNoseMass = ini.GetValueF("Dynamics","nose_mass",thermoNoseMass);
        thermoLangvinFriction = ini.GetValueF("Dynamics","langevin_friction",thermoLangvinFriction);

        // [Parallel Replica]

        paraRepDephaseSteps = ini.GetValueL("Parallel Replica", "dephase_steps", paraRepDephaseSteps);
        paraRepRefine = ini.GetValueB("Parallel Replica", "refine_transition_time", paraRepRefine);
        paraRepAutoStop = ini.GetValueB("Parallel Replica", "stop_after_transition", paraRepAutoStop);
        paraRepCheckPeriod = ini.GetValueL("Parallel Replica", "state_check_period", paraRepCheckPeriod);
        paraRepRecordPeriod = int(0.1*paraRepCheckPeriod);
        paraRepRecordPeriod = ini.GetValueL("Parallel Replica", "state_save_period", paraRepRecordPeriod);
        paraRepRefineAccuracy = ini.GetValueL("Parallel Replica", "bisection_accuracy", paraRepRefineAccuracy);
        paraRepRelaxSteps = ini.GetValueL("Parallel Replica", "post_transition_steps", paraRepRelaxSteps);
        paraRepDephaseLoopStop = ini.GetValueB("Parallel Replica", "dephase_loop_stop", paraRepDephaseLoopStop);
        paraRepDephaseLoopMax = ini.GetValueL("Parallel Replica", "dephase_loop_max", paraRepDephaseLoopMax);

        // [Distributed Replica] //

        drBalanceSteps = ini.GetValueL("Distributed Replica", "balance_steps", drBalanceSteps);
        drSamplingSteps = ini.GetValueL("Distributed Replica", "sampling_steps", drSamplingSteps);
        drTargetTemperature = ini.GetValueF("Distributed Replica", "target_temperature", drTargetTemperature);

        // [Hyperdynamics] //

        bondBoostRMDS = ini.GetValueL("Hyperdynamics","bb_rmd_steps",bondBoostRMDS);
        bondBoostBALS = toLowerCase(ini.GetValue("Hyperdynamics","bb_boost_atomlist",bondBoostBALS));
        bondBoostDVMAX = ini.GetValueF("Hyperdynamics","bb_dvmax",bondBoostDVMAX);
        bondBoostQRR = ini.GetValueF("Hyperdynamics","bb_stretch_threshold",bondBoostQRR );
        bondBoostPRR = ini.GetValueF("Hyperdynamics","bb_ds_curvature",bondBoostPRR );
        bondBoostQcut= ini.GetValueF("Hyperdynamics","bb_rcut",bondBoostQcut);
        biasPotential = toLowerCase(ini.GetValue("Hyperdynamics","bias_potential",biasPotential));

        // [Basin Hopping] //

        basinHoppingMaxDisplacement = ini.GetValueF("Basin Hopping", "max_displacement", basinHoppingMaxDisplacement);
        basinHoppingSteps = ini.GetValueL("Basin Hopping", "steps", basinHoppingSteps);
        basinHoppingQuenchingSteps = ini.GetValueL("Basin Hopping", "quenching_steps", basinHoppingQuenchingSteps);
        basinHoppingSingleAtomDisplace = ini.GetValueB("Basin Hopping", "single_atom_displace", basinHoppingSingleAtomDisplace);
        basinHoppingSignificantStructure = ini.GetValueB("Basin Hopping", "significant_structure", basinHoppingSignificantStructure);
        basinHoppingMaxDisplacementAlgorithm = toLowerCase(ini.GetValue("Basin Hopping", "max_displacement_algorithm", basinHoppingMaxDisplacementAlgorithm));
        basinHoppingDisplacementDistribution = toLowerCase(ini.GetValue("Basin Hopping", "displacement_distribution", basinHoppingDisplacementDistribution));
        basinHoppingSwapProbability = ini.GetValueF("Basin Hopping", "swap_probability", basinHoppingSwapProbability);
        basinHoppingJumpMax = ini.GetValueL("Basin Hopping", "jump_max", basinHoppingJumpMax);
        basinHoppingJumpSteps = ini.GetValueL("Basin Hopping", "jump_steps", basinHoppingJumpSteps);
        basinHoppingInitialMD = ini.GetValueB("Basin Hopping", "initial_md", basinHoppingInitialMD);
        basinHoppingInitialMDTemperature = ini.GetValueF("Basin Hopping", "initial_md_temperature", temperature);

    }
    else
    {
        fprintf(stderr, "Couldn't parse the ini file. Perhaps you are "
                        "using the old style config?\n");
        error = 1;
    }
    return error;
}
