#include"include/UnitCellTemplate.h"



UnitCell InitializeUnitCell() {
    // Define atom templates
    Atom_Template fe_corners = {.atom=1, .spin=5./2};
    int fe_corner_first_nearest_neighbours[][3] = {{2,0,0},{-2,0,0},{0,2,0},{0,-2,0},
                                            {1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0}};
    int fe_corner_second_nearest_neighbours[][3] = {{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},
                                            {-1,1,1},{-1,1,-1},{-1,-1,1},{-1,-1,-1},
                                            {0,0,1},{0,0,-1}};
    fe_corners.first_nearest_neighbours = *fe_corner_first_nearest_neighbours;
    fe_corners.second_nearest_neighbours = *fe_corner_second_nearest_neighbours;
    fe_corners.first_nearest_neighbours_size = 8;
    fe_corners.second_nearest_neighbours_size = 10;

    Atom_Template mo_edgecenters = {.atom=0, .spin=1./2};
    int mo_corner_first_nearest_neighbours[][3] = {{2,0,0},{-2,0,0},{0,2,0},{0,-2,0},
                                            {1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0}};
    int mo_corner_second_nearest_neighbours[][3] = {{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},
                                                {-1,1,1},{-1,1,-1},{-1,-1,1},{-1,-1,-1},
                                                {0,0,1},{0,0,-1}};
    mo_edgecenters.first_nearest_neighbours = *mo_corner_first_nearest_neighbours;
    mo_edgecenters.second_nearest_neighbours = *mo_corner_second_nearest_neighbours;
    mo_edgecenters.first_nearest_neighbours_size = 8;
    mo_edgecenters.second_nearest_neighbours_size = 10;

    Atom_Template fe = {.atom=1, .spin=5./2};
    int fe_second_nearest_neighbours[][3] = {{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},
                                            {-1,1,1},{-1,1,-1},{-1,-1,1},{-1,-1,-1}};
    fe.second_nearest_neighbours = *fe_second_nearest_neighbours;
    fe.first_nearest_neighbours_size = 0;
    fe.second_nearest_neighbours_size = 8;

    Atom_Template mo = {.atom=0, .spin=1./2};
    int mo_first_nearest_neighbours[][3] = {{1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0}};
    int mo_second_nearest_neighbours[][3] = {{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},
                                            {-1,1,1},{-1,1,-1},{-1,-1,1},{-1,-1,-1}};
    mo.first_nearest_neighbours = *mo_first_nearest_neighbours;
    mo.second_nearest_neighbours = *mo_second_nearest_neighbours;
    mo.first_nearest_neighbours_size = 4;
    mo.second_nearest_neighbours_size = 8;

    Atom_Template null = {.atom=BYTEMAX, .spin=0, .first_nearest_neighbours_size=0,
                            .second_nearest_neighbours_size=0};

    // Initialize unit cell
    UnitCell cell = {
        .atoms = {
            {{fe_corners, null, fe_corners}, {null, mo, null}, {fe_corners, null, fe_corners}},
            {{mo_edgecenters, null, mo_edgecenters}, {null, fe, null}, {mo_edgecenters, null, mo_edgecenters}},
            {{fe_corners, null, fe_corners}, {null, mo, null}, {fe_corners, null, fe_corners}}
        }
    };
    cell.no_of_atoms_1 = 2;
    cell.no_of_atoms_2 = 2;
    cell.no_of_atoms_3 = 0;

    return cell;
}