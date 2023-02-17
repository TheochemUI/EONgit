//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef SW_POTENTIAL
#define SW_POTENTIAL

#include "../../Potential.h"

    /** External function implemented in Fortran to calculate interactions between atoms using the Stillinger Weber forcefield
    @param[in]	N           number of atoms
    @param[in]	R           array to positions of the atoms in Angstrom
    @param[out]	F           array used to return the forces resulting from interactions between molecules. Forces are in eV/Angstrom
    @param[out]	U           pointer to energy in eV
    @param[in]  bx, by, bz  pointer to box dimensions in Angstrom
    */
extern "C" {
    void sw_(const long int *N, const double *R, double *F, double *U, const double* bx, const double* by, const double* bz);
}    

/** SW potential.*/
class SW : public Potential{    
public:
// Functions
	// constructor
    SW(Parameters* params) : Potential(params){};
	
    // To satisfy interface
    std::pair<double, AtomMatrix> get_ef(const AtomMatrix pos,
                                         const VectorXi atmnrs,
                                         const Matrix3d m_box) override;
};
#endif

