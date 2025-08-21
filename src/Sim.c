#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "include/tests.h"
#include "include/Lattice.h"
#include "include/Parameters.h"
#include "include/Random.h"
#include "include/UnitCellTemplate.h"

// run simulation, estimate Tc, and write results to file
double estimate_Tc(
    int readings, double temp_start, double temp_step,
    UnitCell cell, float magnetic_field,
    const char *label, FILE *fpt)
{
    float magnetizations[MAX_READINGS];
    float temperature[MAX_READINGS];

    // --- Run simulation for each T ---
    #pragma omp parallel for
    for (int i = 0; i < readings; i++) {
        temperature[i] = temp_start + i * temp_step;
        double beta = 1.0 / temperature[i];

        Lattice *lattice = (Lattice *)malloc(sizeof(Lattice));
        if (lattice == NULL) continue;

        InitalizeLattice(lattice, cell);
        Calculate_Neighbours(lattice, cell);
        CalculateTotalEnergy(lattice, magnetic_field);
        PutOnEquilibrium(lattice, beta, magnetic_field);
        CalulateAverageParameters(lattice, cell, beta, magnetic_field);

        magnetizations[i] = lattice->Magnetization;

        free(lattice);
    }

    // --- Write to CSV ---
    for (int i = 0; i < readings; i++) {
        fprintf(fpt, "%f,%f,%s\n", temperature[i], magnetizations[i], label);
    }

    // --- Central difference slope computation ---
    double dMdT[MAX_READINGS] = {0};
    for (int i = 1; i < readings - 1; ++i) {
        double dT = (temperature[i+1] - temperature[i-1]) * 0.5;
        dMdT[i] = (magnetizations[i+1] - magnetizations[i-1]) / dT;
    }

    // --- Find max slope location ---
    int im = 1;
    double best = fabs(dMdT[1]);
    for (int i = 2; i < readings - 1; ++i) {
        double v = fabs(dMdT[i]);
        if (v > best) { best = v; im = i; }
    }

    return temperature[im];
}

int main() {
    float magnetic_field = 0;
    UnitCell cell = InitializeUnitCell();

    FILE *fpt = fopen("graph.csv", "w+");
    fprintf(fpt, "Temperature,Magnetization,Pass\n");

    // First coarse pass
    double coarse_step = 10.0;
    int coarse_points = MAX_READINGS; // 30 points up to 300 K
    double Tc_est1 = estimate_Tc(coarse_points, coarse_step, coarse_step, cell, magnetic_field, "coarse", fpt);
    printf("First coarse Tc estimate ≈ %f\n", Tc_est1);

    // Second finer pass around Tc_est1
    double fine_step = 1.0;
    double T_low  = Tc_est1 - 10; // search window
    int fine_points = 21;         // Tc_est1-10 to Tc_est1+10
    double Tc_est2 = estimate_Tc(fine_points, T_low, fine_step, cell, magnetic_field, "fine", fpt);
    printf("Refined Tc estimate ≈ %f\n", Tc_est2);

    // Third ultra-fine pass
    double ultra_step = 0.2;
    double T_low2  = Tc_est2 - 2; 
    int ultra_points = 21; // Tc_est2-2 to Tc_est2+2
    double Tc_est3 = estimate_Tc(ultra_points, T_low2, ultra_step, cell, magnetic_field, "ultra", fpt);
    printf("Ultra-fine Tc estimate ≈ %f\n", Tc_est3);

    fclose(fpt);
    return 0;
}
