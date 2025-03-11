#pragma once

#include "Definitions.h"
#include "Lattice.h"

void CalculateTotalEnergy(Lattice* lattice, float magnetic_field);

double SpinFlipEnergyChange(Lattice* lattice,int id, float magentic_field);

void CalculateMagnetizationPerSite(Lattice* lattice,UnitCell cell);

void CalulateAverageParameters(Lattice *lattice, UnitCell cell, double beta, float magnetic_field);