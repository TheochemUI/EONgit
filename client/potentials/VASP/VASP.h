//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef VASP_POTENTIAL
#define VASP_POTENTIAL

#include "../../Potential.h"

class VASP : public Potential
{

    public:
    VASP(Parameters *params): Potential(params);
		~VASP();
        void initialize() {};
        void cleanMemory(void);    
        void force(long N, const double *R, const int *atomicNrs, 
                   double *F, double *U, const double *box);

        std::pair<double, AtomMatrix> get_ef(const AtomMatrix pos,
                                             const VectorXi atmnrs,
                                             const Matrix3d box) override {
          double energy{std::numeric_limits<double>::infinity()};
          long nAtoms{pos.rows()};
          AtomMatrix forces{Eigen::MatrixXd::Zero(nAtoms, 3)};
          this->force(nAtoms, pos.data(), atmnrs.data(), forces.data(), &energy,
                      box.data());
          return std::make_pair(energy, forces);
        };

    private:
        void writePOSCAR(long N, const double *R, const int *atomicNrs,
                         const double *box);
        void readFU(long N, double *F, double *U);
        void spawnVASP();
        bool vaspRunning();
        static bool firstRun;
        static long vaspRunCount;
        static pid_t vaspPID;
};

#endif

