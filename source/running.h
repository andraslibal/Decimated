/*
This part deals with running the simulation

*/
#pragma once

void run_system();
void move_particles();
void write_cmovie_frame();
void write_xmgrace_files();

void calculate_temperature();
void calculate_external();
void calculate_magnetic_field();

void calculate_pinning_forces();
void calculate_temperature_forces();
void calculate_pairwise_forces();
void calculate_external_forces();
void calculate_spin_flipping(int pinning_i,int particle_i,double dxprime);
void calculate_external_forces();

void build_Verlet_list();

void recalculate_vertex_types();
void recalculate_square_vertex_type(int vertex_ID);
void recalculate_hexagonal_vertex_type(int vertex_ID);
void recalculate_square_gs(int vertex_ID);
void recalculate_square_topological_charge_on_vertex(int vertex_ID);

void initialize_statistics();
void calculate_statistics_per_step();
void write_statistics();

void writeout_hopping_file_during_run();
void writeout_hopping_heightmap();
void writeout_hophistogram();


