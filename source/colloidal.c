/* 
Colloidal BD Simulation Code

Date: 2014.05.01
Author: Andras Libal
*/

#include "global_variables.h"
#include "initialization.h"
#include "running.h"
#include "timing.h"
#include "random_numrec.h"

int main(int argc, char *argv[])
{

read_parameters_and_variables(argc, argv);
make_copy_of_parameter_file(argc, argv);
setseed(global_parameters.seedmodifier);

initialize_pinning_sites();
write_contour_file();
//write_triple_contour_file();
initialize_particles();
initialize_pinning_sites_spins();

initialize_vertices();
initialize_vertex_z();
//checkvertexfinder();
initialize_tabulated_forces();
initialize_files();
build_Verlet_list();

//test_by_coloring();

get_starting_time();
run_system();
get_finishing_time();
echo_running_time();

close_files();
//write_testfile();
return 0;
}
