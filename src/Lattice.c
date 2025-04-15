#include <stdio.h>
#include <math.h>
#include "include/Lattice.h"
#include "include/Parameters.h"
#include "include/Random.h"
#include "include/HelperFunctions.h"



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
                if(Latt->atoms[i][j][k].atom==BYTEMAX)
                {
                    Latt->atoms[i][j][k].id = 0xffffffff;
                    continue;
                }
                int l = (i+1)%2+1, m = (j+1)%2+1, n = (k+1)%2+1;
                if(i==0) {l=0;}
                if(j==0) {m=0;}
                if(k==0) {n=0;}
                for(int a=0;a<cell.atoms[n][m][l].first_nearest_neighbours_size;a++)
                {
                    int x = cell.atoms[n][m][l].first_nearest_neighbours[3*a]+i+no_of_atoms_side;
                    int y = cell.atoms[n][m][l].first_nearest_neighbours[3*a+1]+j+no_of_atoms_side;
                    int z = cell.atoms[n][m][l].first_nearest_neighbours[3*a+2]+k+no_of_atoms_side;
                    x = x%no_of_atoms_side;
                    y = y%no_of_atoms_side;
                    z = z%no_of_atoms_side;
                    Latt->atoms[i][j][k].first_nearest_neighbours[a] = Latt->atoms[x][y][z].id;
                }
                for(int a=0;a<cell.atoms[n][m][l].second_nearest_neighbours_size;a++)
                {
                    int x = cell.atoms[n][m][l].second_nearest_neighbours[3*a] +i+  no_of_atoms_side;
                    int y = cell.atoms[n][m][l].second_nearest_neighbours[3*a+1] +j+  no_of_atoms_side;
                    int z = cell.atoms[n][m][l].second_nearest_neighbours[3*a+2] +k+ no_of_atoms_side;
                    x = x%no_of_atoms_side;
                    y = y%no_of_atoms_side;
                    z = z%no_of_atoms_side;
                    Latt->atoms[i][j][k].second_nearest_neighbours[a] = Latt->atoms[x][y][z].id;
                }
            }
        }
    }
}


void PutOnEquilibrium(Lattice *lattice, double beta, float magnetic_field)
{
    int no_of_atoms = NO_OF_ATOMS_SIDE;
    int no_of_spins = LATTICESIZE*LATTICESIZE*LATTICESIZE*4; //no of atoms in a cell
    double delta_E;
    for(int i=0;i<MONTECARLO_CYCLES_AVERAGE;i++)
    {
        for(int j=0;j<no_of_spins;)
        {
            int id = GenerateRandomAtomID();
            delta_E = SpinFlipEnergyChange(lattice,id,magnetic_field);
            if(delta_E==0) continue;
            else if(delta_E<0)
            {
                int i = id % no_of_atoms;
                int j = (id / no_of_atoms) % no_of_atoms;
                int k = (id / (no_of_atoms * no_of_atoms)) % no_of_atoms;
                lattice->atoms[i][j][k].spin *= -1;
                lattice->Energy += delta_E;
            }
            else
            {
                if(exp(-delta_E*beta)>Random())
                {
                    int i = id % no_of_atoms;
                    int j = (id / no_of_atoms) % no_of_atoms;
                    int k = (id / (no_of_atoms * no_of_atoms)) % no_of_atoms;
                    lattice->atoms[i][j][k].spin *= -1;
                    lattice->Energy += delta_E;
                }
            }
            j++;
        }
    }
}


void CreateDeficiency(Lattice * lattice, float DEFECIENCY_PERCENTAGE)
{
    int no_of_atoms = NO_OF_ATOMS_SIDE;
    int deficiency = (int)((DEFECIENCY_PERCENTAGE/100)*(LATTICESIZE*LATTICESIZE*LATTICESIZE*2)); //2 is no of atoms of a type in a unit cell
    for(;deficiency>0;)
    {
        int id = GenerateRandomAtomID();
        uint32_t indices [3];
        CalculateIndicesFromId(id,indices);
        if(lattice->atoms[indices[0]][indices[1]][indices[2]].atom==1)
        {
            lattice->atoms[indices[0]][indices[1]][indices[2]].atom=BYTEMAX;
            lattice->atoms[indices[0]][indices[1]][indices[2]].spin=0;
            deficiency--;
        }
    }
}