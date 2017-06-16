/*
This part deals with initializing the variables from a file
it can read in a file and read/construct/fill the
necessary variables with the appropriate values

alternatively if this is a restart from a saved state it can
also read a restart and place the system in the state at the restart
the difference between the restart and the initialization is only
that with a restart we read all particle position information out
from a restart file, while for an initialization we might calculate
new positions ourselves

*/

#pragma once

void read_parameters_and_variables(int argc, char *argv[]);
void read_parameter_line(char expected_type_specifier[], char expected_subtype_specifier[], char expected_variable_type[], char conversion_yes_no[],
	int *target_int, double *target_double, char target_char[]);
void initialize_diagnostics_file();
void make_copy_of_parameter_file(int argc, char *argv[]);

void initialize_pinning_sites();
void initialize_pinning_sites_spins();
void initialize_square_pinning();
void initialize_hexagonal_pinning();
void initialize_which_sites_will_be_double_defects(int empty_defects_exist,int double_defects_exist);
void initialize_defect_line_square();
void initialize_single_defect_line(int x0, int y0, int defect_line_length,int direction);

void initialize_particles();
void initialize_square_particle();
void initialize_hexagonal_particle();
void initialize_vertices();
void initialize_vertex_z();
void initialize_square_vertices();
void initialize_hexagonal_vertices();

void initialize_files();
void close_files();

void initialize_tabulated_forces();
void write_contour_file();
void write_triple_contour_file();

int check_site_if_decimated(int site_number);
void return_pinningsite_neighbors(int pinningsite, int *neighbors);
void return_pinningsite_vertices(int pinningsite, int *vert1, int *vert2);
void decimator();

void test_by_coloring();
void write_testfile();
void checkvertexfinder();
