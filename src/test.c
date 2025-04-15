#include"include/tests.h"
#include"include/Lattice.h"
#include"include/UnitCellTemplate.h"
#include"include/HelperFunctions.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

// Insert your function and struct definitions here

void test_Calculate_Neighbours() {
    int no_of_atoms_side = NO_OF_ATOMS_SIDE;
    UnitCell cell = InitializeUnitCell();
    Lattice *Latt = (Lattice *)malloc(sizeof(Lattice));
    InitalizeLattice(Latt, cell);
    Calculate_Neighbours(Latt, cell);
    for(int id=0;id<no_of_atoms_side*no_of_atoms_side*no_of_atoms_side;id++)
    {
        int counter = 0;
        uint32_t indices[3];
        CalculateIndicesFromId(id,indices);
        if(Latt->atoms[indices[0]][indices[1]][indices[2]].atom==BYTEMAX) continue;
        for(int i=0;i<no_of_atoms_side;i++)
        {
            for(int j=0;j<no_of_atoms_side;j++)
            {
                for(int k=0;k<no_of_atoms_side;k++)
                {
                    // printf("%d %d %d\n",i,j,k);
                    if(Latt->atoms[i][j][k].atom==BYTEMAX) continue;
                    for(int a=0;a<cell.atoms[0][0][0].first_nearest_neighbours_size;a++)
                    {
                        // printf("%d ",Latt->atoms[i][j][k].first_nearest_neighbours[a]);
                        if(Latt->atoms[i][j][k].first_nearest_neighbours[a] == id) counter++;;
                    }
                    // printf("\n");
                    for(int a=0;a<cell.atoms[0][0][0].second_nearest_neighbours_size;a++)
                    {
                        // printf("%d ",Latt->atoms[i][j][k].second_nearest_neighbours[a]);
                        if(Latt->atoms[i][j][k].second_nearest_neighbours[a] == id) counter++;
                    }
                    // printf("\n");
                }
            }
        }
        assert(counter==18);
    }
    free(Latt);
    printf("All Neighbours are correct");   
}