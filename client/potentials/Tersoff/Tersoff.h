//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef TERSOFF_POTENTIAL
#define TERSOFF_POTENTIAL

#include "../../Potential.h"

    /** External function implemented in Fortran to calculate interactions between atoms using the Tersoff forcefield
    @param[in]	N           Number of atoms
    @param[in]	R           Array to positions of the atoms in Angstrom
    @param[out]	F           Array used to return the forces between atoms, in eV/Angstrom
    @param[out]	U           Pointer to energy in eV
    @param[in]  bx, by, bz  Pointer to box dimensions in Angstrom
    */
extern "C" {
    void tersoff_(const long int *N, const double *R, double *F, double *U, const double* bx, const double* by, const double* bz);
}    

/** Tersoff potential */
class Tersoff : public Potential{    
private:
    Parameters *parameters;
public:
// Functions
	// constructor
    Tersoff(Parameters *p) : Potential(p), parameters{p} {}

    // To satisfy interface
    void initialize(void);
    void cleanMemory(void);
    void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
};
#endif

