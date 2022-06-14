#ifndef GPRHELPERS_H
#define GPRHELPERS_H

#include "Matter.h"
#include "Parameters.h"
#include "HelperFunctions.h"
#include "Log.h"

#include "subprojects/gprdimer/structures/Structures.h"
#include "subprojects/gprdimer/gpr/auxiliary/ProblemSetUp.h"
#include "subprojects/gprdimer/gpr/ml/GaussianProcessRegression.h"
#include "subprojects/gprdimer/data_types/Field.h"

namespace helper_functions {
    /**
     * \brief Create a parameters object for gpr_dimer
     *
     * TODO: Get the cell size from a Matter object
     *
     * @param *parameters An EON parameters object
     */
    gpr::InputParameters eon_parameters_to_gpr(Parameters *parameters);

    /**
     * \brief Create a parameters object for gpr_pot
     *
     * TODO: Get the cell size from a Matter object
     *
     * @param *parameters An EON parameters object
     */
    gpr::InputParameters eon_parameters_to_gprpot(Parameters *parameters);

    /**
     * \brief Create a configuration of atoms for gpr_dimer
     *
     * @param *Matter An EON Matter object
     */
    gpr::AtomsConfiguration eon_matter_to_atmconf(Matter* matter);

    gpr::Field<double> generateAtomsConfigField(const Matter& mat);
    /**
     * \brief Create an initial Observation object for gpr_dimer
     *
     * Note that this is essentially only for the setup, and it does not
     * actually append or handle the Observation structure other than for
     * initialization of atomic gp dimer
     *
     * @param matter An EON Matter object
     */
    gpr::Observation eon_matter_to_init_obs(Matter& matter);

    /**
     * \brief Setup initial path
     *
     * This is essentially what the NEB initialization does, however, we need
     * this here to prepare the initial observations for the GPR surface
     * */
    std::vector<Matter> prepInitialPath(
                   Parameters *params,
                   std::string fname_reactant="reactant.con"s,
                   std::string fname_product="product.con"s);
    /**
     * \brief Setup initial observations from path
     *
     * This starts with the construction of an observation object and then uses
     * the linearly interpolated NEB path from prepInitialPath to populate the
     * observation object
     * */
    gpr::Observation prepInitialObs(std::vector<Matter> &vecmat);
    /**
     * \brief Generate a conf_info object from a matter object
     *
     * This actually calls eon_matter_to_atmconf under the hood
     * The difference being that this function sets the frozen atoms
     *  \param matter The object to use as a constructor
     * */
    std::pair<gpr::AtomsConfiguration, gpr::Coord> eon_matter_to_frozen_conf_info(Matter* matter, double activeRadius);
    /**
     * \brief Initializes a GPR potential
     *
     * The logic here is to simply prepare the GPR for training, it is a
     * super-set of EON parameter conversions
     *
     * Note that it is **initialized** but NOT trained.
     * To train this run setHyperParameters and optimize on the returned object
     * */
    gpr::GaussianProcessRegression& initializeGPR(gpr::GaussianProcessRegression& gprfunc,
                                                                    gpr::AtomsConfiguration& atoms_config,
                                                                    gpr::Observation& obsPath,
                                                                    std::pair<Parameters, Matter>& eon_matter_params);
    } // namespace helper_functions
#endif /* GPRHELPERS_H */
