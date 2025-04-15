#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include "include/tests.h"
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
    // test_Calculate_Neighbours();
    int readings = 130;  // Change this as needed
    int def_readings = 12;
    float magnetic_field = 0;
    float magnetizations[def_readings][MAX_READINGS];
    float temperature[MAX_READINGS] = {0};
    double Energies[MAX_READINGS] = {0};
    FILE *fpt;
    fpt = fopen("graph.csv", "w+");
    fprintf(fpt, "Temperature,");
    for(int j = 0; j < def_readings; j++)
    {
        fprintf(fpt, "%f,", (float)(j+16)/2.);
    }
    fprintf(fpt, "\n");
    UnitCell cell = InitializeUnitCell();

    #pragma omp parallel
    {
        #pragma omp single
        printf("OpenMP is using %d threads.\n", omp_get_num_threads());
    }
    for(int j=0;j<def_readings;j++)
    {
        #pragma omp parallel for
        for (int i = 0; i < readings; i++) {
            temperature[i] = (i+1) * 10;
            double beta = 1. / temperature[i];
            
            // Initialize UnitCell and Lattice
            Lattice *lattice = (Lattice *)malloc(sizeof(Lattice));
            
            if (lattice == NULL) {
                printf("Memory allocation failed for lattice!\n");
                continue;  // Proceed to next iteration if malloc fails
            }
            InitalizeLattice(lattice, cell);
            Calculate_Neighbours(lattice, cell);
            CreateDeficiency(lattice,(float)(j+12)/2.);
            CalculateTotalEnergy(lattice, magnetic_field);
            PutOnEquilibrium(lattice, beta, magnetic_field);
            CalulateAverageParameters(lattice, cell, beta, magnetic_field);
            
            Energies[i] = lattice->Energy / (LATTICESIZE * LATTICESIZE * LATTICESIZE * (cell.no_of_atoms_1 + cell.no_of_atoms_2));
            magnetizations[j][i] = lattice->Magnetization;
            free(lattice);  // Free the memory after using it
        }
    }
    // double Energy_mean = 0;
    // double Energy_variance = 0;
    // for (int i = 0; i < readings; i++) 
    // {
    //     Energy_mean+=Energies[i];
    //     Energy_variance+=pow(Energies[i],2);
    // }
    // Energy_mean = Energy_mean/readings;
    // Energy_variance = Energy_variance/readings;
    // printf("%f\n",Energy_variance - pow(Energy_mean,2));

    // Write results to file
    for (int i = 0; i < readings; i++) {
        fprintf(fpt, "%f,", temperature[i]);
        for(int j = 0; j < def_readings; j++)
        {
            fprintf(fpt, "%f,", magnetizations[j][i]);
        }
        fprintf(fpt, "\n");
    }
    fclose(fpt);
    return 0;
}