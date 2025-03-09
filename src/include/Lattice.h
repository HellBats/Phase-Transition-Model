#pragma once
#include <stdint.h>
#include "include/Definitions.h"

typedef struct 
{
    uint8_t atom;// Atom type array BYTEMAX means no atom there
    float spin; // This array accounts for spin of atoms
    uint32_t id;
    int first_nearest_neighbours[14];
    int second_nearest_neighbours[14];
    uint8_t first_nearest_neighbours_size;
    uint8_t second_nearest_neighbours_size;
}Atom;


typedef struct 
{
    uint8_t atom;// Atom type array BYTEMAX means no atom there
    float spin; // This array accounts for spin of atoms
    void* first_nearest_neighbours;
    void* second_nearest_neighbours;
    uint8_t first_nearest_neighbours_size;
    uint8_t second_nearest_neighbours_size;
}Atom_Template;


typedef struct UnitCell
{
    /* This structure only acounts for atoms at corners, face centers, edge centers and body 
    center. */
    Atom_Template atoms[3][3][3];
    uint8_t no_of_atoms_1,no_of_atoms_2,no_of_atoms_3; 
}UnitCell;


typedef struct 
{
    Atom atoms[NO_OF_ATOMS_SIDE][NO_OF_ATOMS_SIDE][NO_OF_ATOMS_SIDE];
    float Energy;
    double Magnetization;
}Lattice;

void InitalizeLattice(Lattice* Lattice,
     UnitCell cell);

void Calculate_Neighbours(Lattice* Lattice,
        UnitCell cell);

void PutOnEquilibrium(Lattice *lattice, UnitCell cell, double beta, float magnetic_field);