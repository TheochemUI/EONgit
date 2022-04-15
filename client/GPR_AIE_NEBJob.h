#ifndef GPR_AIE_NEBJOB_H_
#define GPR_AIE_NEBJOB_H_

#include "Matter.h"
#include "Parameters.h"
#include "NudgedElasticBand.h"
#include "SafeNudgedElasticBand.h"
#include "Job.h"
#include "GPRHelpers.h"

class GPR_AIE_NEBJob : public Job {

    public:

        GPR_AIE_NEBJob(Parameters *parametersPassed);
        ~GPR_AIE_NEBJob(void);
        std::vector<std::string> run(void);

    private:

        // functions
        void printEndState(int status);
        void saveData(int status, SafeNudgedElasticBand *neb);

        // variables
        std::vector<std::string> returnFiles;
        Parameters *parameters;
        std::unique_ptr<gpr::GaussianProcessRegression> gprfunc;
        size_t fCallsNEB;
        size_t fCallsGPR;

};

#endif // GPR_AIE_NEBJOB_H_
