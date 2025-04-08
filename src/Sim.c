#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include "include/Lattice.h"
#include "include/Parameters.h"
#include "include/Random.h"
#include "include/UnitCellTemplate.h"

void Print2DLattice(Lattice *Lattice)
{
    int positive = 0, negative = 0;
    for (int i = 0; i < NO_OF_ATOMS_SIDE; i++)
    {
        int j = 0;
        for (int k = 0; k < NO_OF_ATOMS_SIDE; k++)
        {
            printf("%f ", Lattice->atoms[0][i][k].spin);
            if (Lattice->atoms[0][i][k].spin > 0)
                positive++;
            else if (Lattice->atoms[0][i][k].spin < 0)
                negative++;
        }
        printf("\n");
    }
    printf("%d %d\n", positive, negative);
}

int main()
{
    #pragma omp parallel
    {
    #pragma omp single
        printf("OpenMP is using %d threads.\n", omp_get_num_threads());
    }
    uint8_t readings = 120;
    float magnetic_field = 0;
    float magnetization[120] = {0};
    float temperature[120] = {0};
    double Energies[120] = {0};
    FILE *fpt;
    fpt = fopen("graph.csv", "w+");
    fprintf(fpt, "Temprature,Manetization,Energies\n");
    #pragma omp parallel for
    for (int i = 0; i < readings; i++)
    {
        temperature[i] = ((i + 1) * 10);
        double beta = 1. / (temperature[i]);
        Lattice *lattice;
        UnitCell cell = InitializeUnitCell();
        lattice = (Lattice *)malloc(sizeof(Lattice));
        InitalizeLattice(lattice, cell);
        Calculate_Neighbours(lattice, cell);
        // CreateDeficiency(lattice);
        CalculateTotalEnergy(lattice, magnetic_field);
        PutOnEquilibrium(lattice, beta, magnetic_field);
        CalulateAverageParameters(lattice, cell, beta, magnetic_field);
        Energies[i] = lattice->Energy / (LATTICESIZE * LATTICESIZE * LATTICESIZE * (cell.no_of_atoms_1 + cell.no_of_atoms_2));
        magnetization[i] = (double)(lattice->Magnetization);
        free(lattice);
    }
    for (int i = 0; i < readings; i++)
    {
        fprintf(fpt, "%f,%f,%f\n", temperature[i], magnetization[i], Energies[i]);
    }
    fclose(fpt);
}
