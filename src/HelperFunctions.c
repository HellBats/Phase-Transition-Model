
#include "include/HelperFunctions.h"

void CalculateIndicesFromId(uint32_t id, uint32_t indices[3])
{
    int no_of_atoms = NO_OF_ATOMS_SIDE;
    indices[0]= id % no_of_atoms;
    indices[1] = (id / no_of_atoms) % no_of_atoms;
    indices[2] = (id / (no_of_atoms * no_of_atoms)) % no_of_atoms;
}