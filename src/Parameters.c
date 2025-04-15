#include <stdio.h>
#include <math.h>
#include "include/Parameters.h"
#include "include/Random.h"


void CalculateMagnetizationPerSite(Lattice* lattice,UnitCell cell)
{
    int no_of_atoms = NO_OF_ATOMS_SIDE;
    int no_of_spins = LATTICESIZE*LATTICESIZE*LATTICESIZE*(cell.no_of_atoms_1+cell.no_of_atoms_2);
    double magnetization = 0;
    for (int i = 0; i < no_of_atoms; i++)
    {
        for (int j = 0; j < no_of_atoms; j++)
        {
            for (int k = 0; k < no_of_atoms; k++)
            {
                Atom *current_atom = &lattice->atoms[i][j][k];

                // Skip empty sites (no neighbors, spin = 0)
                if (current_atom->spin == 0 || (current_atom->first_nearest_neighbours_size == 0 &&
                                                current_atom->second_nearest_neighbours_size == 0))
                {
                    continue;
                }
                magnetization += (double)current_atom->spin;
            }
        }
    }
    lattice->Magnetization = magnetization/no_of_spins;
}

void CalulateAverageParameters(Lattice *lattice, UnitCell cell, double beta, float magnetic_field)
{
    int no_of_atoms = NO_OF_ATOMS_SIDE;
    double delta_E;
    double average_energy = 0;
    double average_magnetization = 0;
    int no_of_spins = LATTICESIZE*LATTICESIZE*LATTICESIZE*(cell.no_of_atoms_1+cell.no_of_atoms_2);
    CalculateMagnetizationPerSite(lattice,cell);
    for(int i=0;i<(MONTECARLO_CYCLES_EQ);i++)
    {
        for(int j=0;j<(no_of_spins);)
        {
            int id = GenerateRandomAtomID();
            delta_E = SpinFlipEnergyChange(lattice,id,magnetic_field);
            if (delta_E==0) continue;
            else if(delta_E<0)
            {
                int i = id % no_of_atoms;
                int j = (id / no_of_atoms) % no_of_atoms;
                int k = (id / (no_of_atoms * no_of_atoms)) % no_of_atoms;
                lattice->atoms[i][j][k].spin *= -1;
                float spin = lattice->atoms[i][j][k].spin;
                lattice->Energy += delta_E;
                lattice->Magnetization += (2*spin)/(no_of_spins);
            }
            else
            {
                if(exp(-delta_E*beta)>Random())
                {
                    int i = id % no_of_atoms;
                    int j = (id / no_of_atoms) % no_of_atoms;
                    int k = (id / (no_of_atoms * no_of_atoms)) % no_of_atoms;
                    lattice->atoms[i][j][k].spin *= -1;
                    float spin = lattice->atoms[i][j][k].spin;
                    lattice->Energy += delta_E;
                    lattice->Magnetization += (2*spin)/(no_of_spins);
                }
            }
            if((j%SAMPLING_SIZE)==0)
            {
                average_energy+=lattice->Energy;
                average_magnetization+=lattice->Magnetization;
            }
            j++;
        }
    }
    int total_samples = (MONTECARLO_CYCLES_EQ*no_of_spins)/SAMPLING_SIZE;
    lattice->Energy = average_energy/total_samples;
    lattice->Magnetization = average_magnetization/total_samples;
}


void CalculateTotalEnergy(Lattice* lattice, float magnetic_field)
{
    float J1_Array[2][2] = {{J1MoMo, J1FeMo}, {J1FeMo, J1FeFe}};
    float J2_Array[2][2] = {{J2MoMo, J2FeMo}, {J2FeMo, J2FeFe}};
    
    lattice->Energy = 0;
    int no_of_atoms = NO_OF_ATOMS_SIDE;
    
    for (int i = 0; i < no_of_atoms; i++)
    {
        for (int j = 0; j < no_of_atoms; j++)
        {
            for (int k = 0; k < no_of_atoms; k++)
            {
                Atom *current_atom = &lattice->atoms[i][j][k];

                // Skip empty sites (no neighbors, spin = 0)
                if (current_atom->spin == 0 || (current_atom->first_nearest_neighbours_size == 0 &&
                                                current_atom->second_nearest_neighbours_size == 0))
                {
                    continue;
                }

                uint8_t current_type = current_atom->atom;
                float current_spin = current_atom->spin;

                // First Nearest Neighbors
                for (int a = 0; a < current_atom->first_nearest_neighbours_size; a++)
                {
                    int id = current_atom->first_nearest_neighbours[a];

                    int ni = id % (no_of_atoms);
                    int nj = (id / no_of_atoms) % no_of_atoms;
                    int nk = (id / (no_of_atoms * no_of_atoms)) % no_of_atoms;

                    Atom *neighbor = &lattice->atoms[ni][nj][nk];

                    if (neighbor->spin != 0) // Ensure neighbor isn't empty
                    {
                        lattice->Energy += J1_Array[current_type][neighbor->atom] * current_spin * neighbor->spin;
                    }
                }

                // Second Nearest Neighbors
                for (int a = 0; a < current_atom->second_nearest_neighbours_size; a++)
                {
                    int id = current_atom->second_nearest_neighbours[a];

                    int ni = id % (no_of_atoms);
                    int nj = (id / no_of_atoms) % no_of_atoms;
                    int nk = (id / (no_of_atoms * no_of_atoms)) % no_of_atoms;

                    Atom *neighbor = &lattice->atoms[ni][nj][nk];

                    if (neighbor->spin != 0) // Ensure neighbor isn't empty
                    {
                        lattice->Energy += J2_Array[current_type][neighbor->atom] * current_spin * neighbor->spin;
                    }
                }
                lattice->Energy += 2*magnetic_field*current_spin;
            }
        }
    }
    lattice->Energy /= -2; // Avoid double counting
}

double SpinFlipEnergyChange(Lattice* lattice, int id, float magnetic_field)
{
    // Removed negative sign from calculation of energy (flipped the sign)
    int no_of_atoms_side = NO_OF_ATOMS_SIDE;
    float J1_Array[2][2] = {{J1MoMo, J1FeMo}, {J1FeMo, J1FeFe}};
    float J2_Array[2][2] = {{J2MoMo, J2FeMo}, {J2FeMo, J2FeFe}};
    
    // Compute 3D indices once
    int i = id % no_of_atoms_side;
    int j = (id / no_of_atoms_side) % no_of_atoms_side;
    int k = (id / (no_of_atoms_side * no_of_atoms_side)) % no_of_atoms_side;

    Atom *atom1 = &lattice->atoms[i][j][k];

    // Skip empty sites
    if (atom1->spin == 0) return 0;

    float spin1 = atom1->spin;
    uint8_t type1 = atom1->atom;
    float Energy = 0;

    // First Nearest Neighbors
    for (int a = 0; a < atom1->first_nearest_neighbours_size; a++)
    {
        int id2 = atom1->first_nearest_neighbours[a];

        int ni = id2 % no_of_atoms_side;
        int nj = (id2 / no_of_atoms_side) % no_of_atoms_side;
        int nk = (id2 / (no_of_atoms_side * no_of_atoms_side)) % no_of_atoms_side;

        Atom *atom2 = &lattice->atoms[ni][nj][nk];

        if (atom2->spin != 0) // Ensure neighbor isn't empty
        {
            Energy += J1_Array[type1][atom2->atom] * spin1 * atom2->spin;
        }
    }

    // Second Nearest Neighbors
    for (int a = 0; a < atom1->second_nearest_neighbours_size; a++)
    {
        int id2 = atom1->second_nearest_neighbours[a];

        int ni = id2 % no_of_atoms_side;
        int nj = (id2 / no_of_atoms_side) % no_of_atoms_side;
        int nk = (id2 / (no_of_atoms_side * no_of_atoms_side)) % no_of_atoms_side;

        Atom *atom2 = &lattice->atoms[ni][nj][nk];

        if (atom2->spin != 0) // Ensure neighbor isn't empty
        {
            Energy += J2_Array[type1][atom2->atom] * spin1 * atom2->spin;
        }
    }
    Energy += magnetic_field*spin1;
    return (2*Energy);
}
