//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include"GPRPotential.h"
#include"../../subprojects/gprdimer/structures/Structures.h"
#include"../../subprojects/gprdimer/gpr/auxiliary/AdditionalFunctionality.h"

namespace {

    const char *elementArray[] = {"Unknown", "H","He","Li","Be","B","C","N","O",
           "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc",
           "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
           "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag",
           "Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
           "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta",
           "W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
           "Fr","Ra","Ac","Th","Pa","U", NULL};

    // guess the atom type from the atomic mass,
    std::string mass2atom(double atomicmass) {
        return elementArray[int(atomicmass+.5)];
    }


    int symbol2atomicNumber(char const * symbol)
    {
        int i=0;

        while (elementArray[i] != NULL) {
            if (strcmp(symbol, elementArray[i]) == 0) {
                return i;
            }
            i++;
        }
        // invalid symbol
        return -1;
    }

    char const *atomicNumber2symbol(int n)
    {
        return elementArray[n];
    }
}

GPRPotential::GPRPotential(Parameters *p){
    gpr_model = nullptr;
}

void GPRPotential::registerGPRObject(gpr::GaussianProcessRegression *_gpr_model){
    gpr_model = _gpr_model;
}

void GPRPotential::initialize(void){
}

void GPRPotential::cleanMemory(void){
}

std::pair<double, AtomMatrix> GPRPotential::force(AtomMatrix pos, size_t nAtoms,  int nImages){
    // const int nFreeAtoms = posFree.rows();
    // std::cout<<"Hello from GPRPot with "<< nFreeAtoms << " free and "<<nAtoms << "\n";
    gpr::Observation obs;
    // TODO: Be better with the number of images
    obs.R.resize(1, 3 * nAtoms);
    obs.G.resize(1, 3 * nAtoms);
    obs.E.resize(1);
  // for(size_t idx = 0; idx < nAtoms; idx++) {
      obs.R.assignFromEigenMatrix(pos);
  // }

    // for (size_t idx{0}; idx < nFreeAtoms*3; ++idx){
        // std::cout<<obs.R.extractEigenMatrix()<<" ";
    // }
    // TODO: Benchmark this, see Potential.cpp

    // See GPRTrainTest.cpp for the functions to be called before this
    // std::cout<<"Energy before GPR "<<obs.E.extractEigenMatrix()(0)<<"\n";
    this->gpr_model->calculatePotential(obs);
    double calc_energy = obs.E.extractEigenMatrix()(0);
    AtomMatrix calc_forces = obs.G.extractEigenMatrix() * -1;

    //         double energy = 0;
    //  AtomMatrix forces = AtomMatrix::Constant(nFreeAtoms, 3, 0);
    // std::cout<<"Energy "<<calc_energy<<"\n";
    return std::make_pair(calc_energy, calc_forces);
    // return std::make_pair(energy, forces);
}

// pointer to number of atoms, pointer to array of positions	
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void GPRPotential::force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box, int nImages){
    gpr::Observation obs;
    // TODO: Be better with the number of images
    obs.R.resize(1, 3 * N);
    obs.G.resize(1, 3 * N);
    obs.E.resize(1);
    for (size_t idx{0}; idx < N*3; idx++){
        obs.R[idx] = R[idx];
        obs.G[idx] = 0.0;
    }
    this->gpr_model->calculatePotential(obs);
    double calc_energy = obs.E.extractEigenMatrix()(0);
    AtomMatrix calc_forces = obs.G.extractEigenMatrix() * -1;
    *U = calc_energy;
    std::cout<<"energy "<<*U<<"\n";
    for (size_t idx{0}; idx < N*3; idx++){
        F[idx] = calc_forces.data()[idx];
    }
}
