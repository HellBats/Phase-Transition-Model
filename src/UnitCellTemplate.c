#include"include/UnitCellTemplate.h"



UnitCell InitializeUnitCell() {
    // Define atom templates
    Atom_Template fe = {.atom=1, .spin=5./2};
    int fe_corner_first_nearest_neighbours[24] = {2,0,0,-2,0,0,0,2,0,0,-2,0,
                                            1,1,0,1,-1,0,-1,1,0,-1,-1,0};
    int fe_corner_second_nearest_neighbours[30] = {1,1,1,1,1,-1,1,-1,1,1,-1,-1,
                                            -1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1,
                                            0,0,1,0,0,-1};
    fe.first_nearest_neighbours_size = 8;
    fe.second_nearest_neighbours_size = 10;
    for(int i=0;i<24;i++) fe.first_nearest_neighbours[i] = fe_corner_first_nearest_neighbours[i];
    for(int i=0;i<30;i++) fe.second_nearest_neighbours[i] = fe_corner_second_nearest_neighbours[i];

    Atom_Template mo = {.atom=0, .spin=1./2};
    int mo_corner_first_nearest_neighbours[24] = {2,0,0,-2,0,0,0,2,0,0,-2,0,
                                            1,1,0,1,-1,0,-1,1,0,-1,-1,0};
    int mo_corner_second_nearest_neighbours[30] = {1,1,1,1,1,-1,1,-1,1,1,-1,-1,
                                                -1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1,
                                                0,0,1,0,0,-1};
    mo.first_nearest_neighbours_size = 8;
    mo.second_nearest_neighbours_size = 10;
    for(int i=0;i<24;i++) mo.first_nearest_neighbours[i] = mo_corner_first_nearest_neighbours[i];
    for(int i=0;i<30;i++) mo.second_nearest_neighbours[i] = mo_corner_second_nearest_neighbours[i];

    Atom_Template null = {.atom=BYTEMAX, .spin=0, .first_nearest_neighbours_size=0,
                            .second_nearest_neighbours_size=0};

    // Initialize unit cell
    UnitCell cell = {
        .atoms = {
            {{fe, null, fe}, {null, mo, null}, {fe, null, fe}},
            {{mo, null, mo}, {null, fe, null}, {mo, null, mo}},
            {{fe, null, fe}, {null, mo, null}, {fe, null, fe}}
        }
    };
    cell.no_of_atoms_1 = 2;
    cell.no_of_atoms_2 = 2;
    cell.no_of_atoms_3 = 0;

    return cell;
}