#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include "include/Plots.h"
#include "include/Lattice.h"
#include "include/Parameters.h"
#include "include/Random.h"
#include "include/UnitCellTemplate.h"

void Print2DLattice(Lattice Lattice)
{
    int positive = 0,negative = 0;
    for(int i=0;i<NO_OF_ATOMS_SIDE;i++)
    {
            int j = 0;
            for(int k=0;k<NO_OF_ATOMS_SIDE;k++)
            {
                printf("%f ",Lattice.atoms[0][i][k].spin);
                if(Lattice.atoms[0][i][k].spin>0) positive++;
                else if(Lattice.atoms[0][i][k].spin<0) negative++;
            }
        printf("\n");
    }
    printf("%d %d\n",positive,negative);
}


int main()
{
    #pragma omp parallel
    {
        #pragma omp single
        printf("OpenMP is using %d threads.\n", omp_get_num_threads());
    }
    int readings = 4;
    float magnetic_field = 0;
    int no_of_atoms = NO_OF_ATOMS_SIDE;
    double magnetization[100] = {0};
    double temperature[100] = {0};
    FILE *fpt;
    fpt = fopen("graph.csv","w+");
    fprintf(fpt,"Temprature,Manetization\n");
    #pragma omp parallel for
    for(int i=0;i<readings;i++)
    {
        temperature[i] = (double)((i+20)*10);
        double beta = 1./(temperature[i]);
        UnitCell cell = InitializeUnitCell();
        Lattice lattice = {.Energy=0, .Magnetization = 0};
        PutOnEquilibrium(&lattice,cell,beta,magnetic_field);
        CalulateAverageParameters(&lattice,cell,beta,magnetic_field);
        magnetization[i] = (double)(lattice.Magnetization);
            
    }
    for(int i=0;i<readings;i++)
    {
        fprintf(fpt,"%f,%f\n",temperature[i],magnetization[i]);
    }
    fclose(fpt);
    // int size = sizeof(magnetization) / sizeof(magnetization[0]);
    // printf("%d",size);
    // plotLineGraph(temperature, magnetization, size, "line_plot.png");
}
