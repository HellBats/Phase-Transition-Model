#include "include/Random.h"

static gsl_rng *rng = NULL;
uint32_t GenerateRandomAtomID()
{
    

    if(rng==NULL)
    {
        rng = gsl_rng_alloc(gsl_rng_mt19937);

        gsl_rng_set(rng,time(NULL));
    }
    uint32_t no_of_atoms_side = NO_OF_ATOMS_SIDE;
    uint32_t random_id = gsl_rng_uniform_int(rng,(no_of_atoms_side*no_of_atoms_side*no_of_atoms_side));
    return random_id;
}   


double Random()
{
    if(rng==NULL)
    {
        rng = gsl_rng_alloc(gsl_rng_mt19937);

        gsl_rng_set(rng,time(NULL));
    }
    double random = gsl_rng_uniform(rng);
    return random;
}