#include <stdio.h>
#include <math.h>
#include "include/Lattice.h"
#include "include/Parameters.h"
#include "include/Random.h"


void InitalizeLattice(Lattice* Latt,UnitCell cell)
{
    int no_of_atoms_side = NO_OF_ATOMS_SIDE; 
    for(int i=0;i<no_of_atoms_side;i++)
    {
        for(int j=0;j<no_of_atoms_side;j++)
        {
            for(int k=0;k<no_of_atoms_side;k++)
            {
                int l = (i+1)%2+1, m = (j+1)%2+1, n = (k+1)%2+1;
                if(i==0) {l=0;}
                if(j==0) {m=0;}
                if(k==0) {n=0;}
                Latt->atoms[i][j][k].atom = cell.atoms[n][m][l].atom;
                Latt->atoms[i][j][k].spin = cell.atoms[n][m][l].spin;
                Latt->atoms[i][j][k].first_nearest_neighbours_size = cell.atoms[n][m][l].first_nearest_neighbours_size;
                Latt->atoms[i][j][k].second_nearest_neighbours_size = cell.atoms[n][m][l].second_nearest_neighbours_size;
                Latt->atoms[i][j][k].id = k*no_of_atoms_side*no_of_atoms_side+j*no_of_atoms_side+i;
            }
        }
    }
}


void Calculate_Neighbours(Lattice* Latt,
    UnitCell cell)
{
    int no_of_atoms_side = NO_OF_ATOMS_SIDE; 
    for(int i=0;i<no_of_atoms_side;i++)
    {
        for(int j=0;j<no_of_atoms_side;j++)
        {
            for(int k=0;k<no_of_atoms_side;k++)
            {
                int l = (i+1)%2+1, m = (j+1)%2+1, n = (k+1)%2+1;
                if(i==0) {l=0;}
                if(j==0) {m=0;}
                if(k==0) {n=0;}
                for(int a=0;a<cell.atoms[n][m][l].first_nearest_neighbours_size;a++)
                {
                    int x = ((*((int*)(cell.atoms[n][m][l].first_nearest_neighbours)+3*a) + i) +
                    no_of_atoms_side);
                    int y = ((*((int*)(cell.atoms[n][m][l].first_nearest_neighbours)+3*a+1) + j) +
                    no_of_atoms_side);
                    int z = ((*((int*)(cell.atoms[n][m][l].first_nearest_neighbours)+3*a+2) + k) +
                    no_of_atoms_side);
                    x = x%no_of_atoms_side;
                    y = y%no_of_atoms_side;
                    z = z%no_of_atoms_side;
                    Latt->atoms[i][j][k].first_nearest_neighbours[a] = Latt->atoms[x][y][z].id;
                }
                for(int a=0;a<cell.atoms[n][m][l].second_nearest_neighbours_size;a++)
                {
                    int x = ((*((int*)(cell.atoms[n][m][l].second_nearest_neighbours)+3*a) + i) +
                    no_of_atoms_side);
                    int y = ((*((int*)(cell.atoms[n][m][l].second_nearest_neighbours)+3*a+1) + j) +
                    no_of_atoms_side);
                    int z = ((*((int*)(cell.atoms[n][m][l].second_nearest_neighbours)+3*a+2) + k) +
                    no_of_atoms_side);
                    x = x%no_of_atoms_side;
                    y = y%no_of_atoms_side;
                    z = z%no_of_atoms_side;
                    Latt->atoms[i][j][k].second_nearest_neighbours[a] = Latt->atoms[x][y][z].id;
                }
            }
        }
    }
}


void PutOnEquilibrium(Lattice *lattice, UnitCell cell, double beta, float magnetic_field)
{
    int no_of_atoms = NO_OF_ATOMS_SIDE;
    double delta_E;
    InitalizeLattice(lattice,cell);
    Calculate_Neighbours(lattice,cell);
    CalculateTotalEnergy(lattice,magnetic_field);
    for(int i=0;i<NO_OF_TRIALS;)
    {
        int id = GenerateRandomAtomID();
        delta_E = SpinFlipEnergyChange(lattice,id,magnetic_field);
        if(delta_E==0) continue;
        else if(delta_E<0)
        {
            lattice->atoms[id%no_of_atoms][(id/no_of_atoms)%no_of_atoms]
            [((id/no_of_atoms)/no_of_atoms)%no_of_atoms].spin*=-1;
            lattice->Energy += delta_E;
        }
        else
        {
            if(exp(-delta_E*beta)>Random())
            {
                lattice->atoms[id%no_of_atoms][(id/no_of_atoms)%no_of_atoms]
                [((id/no_of_atoms)/no_of_atoms)%no_of_atoms].spin*=-1;
                lattice->Energy += delta_E;
            }
        }
        i++;
    }
}