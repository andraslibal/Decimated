#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

struct parameters_struct 	
	{
	char 	parameter_file_name[255];
	FILE 	*parameter_file;

	int 	 initialization_type;
	int seedmodifier;

	char 	running_name[255];

	char 	diagnostics_file_name[255];
	FILE 	*diagnostics_file;
	int 	diagnostics_file_ready;

	char 	movie_file_name[255];
	FILE 	*movie_file;
	char 	statistics_file_name[255];
	FILE 	*statistics_file;
	char 	save_file_name[255];
	FILE 	*save_file;

	int 	total_runtime;
	int 	movie_write_time;
	int 	echo_write_time;
	int 	statistics_write_time;

	double 	dt;					//this will be in actual seconds

	int 	time;					//this is still in MD steps (time*dt is the ellapsed time)

	double  boltzmann_constant;			//for real units 0.0138064852 expressed in pN*nm/K
	double	kbT;					//at 16C this is about 4 pN*nm
	double  diffusion_constant;			//possible to set from parameter file	 about 36000 nm^2/s	
	double  mobility;				//calculated from the diffusion constant around 9000 nm/pN/s 
	double  magnetic_susceptibility;		//calculated from the diffusion constant 
	double	magnetic_field;				//in mT, about 10mT to 30mT
    double  initial_magnetic_field;
	double	magnetic_permeability;			//4pi*10^5 pN/A^2	

	double 	temperature;
	double	thermal_force_term;			//in pN, multiplied by gasdev() to get the force

	int 	N_temperature_set_points;
	double 	set_temperature[100];
	int 	set_temperature_time[100];

	double	external_fx, external_fy;
	int 	N_external_set_points;
	double 	set_external_fx[100];
	double 	set_external_fy[100];
	int 	set_external_time[100];

	double  SX,SY;

	int 	recalculate_SX, recalculate_SY;
    //1 means recalculate is needed from pins - spin ice
    //0 means do not recalculate teh SX SY, it is given

	int 	boundarycondition_x, boundarycondition_y;
    //0 means there is no PBC - spin ice does noo need it
    //1 means there is PBC  - not yet implemented

	char 	pinning_setup_type[255];				//spin ice
	char	pinning_lattice_type[255];				//square, hexagonal, or special

	int 	pinning_lattice_Nx, pinning_lattice_Ny;
	double 	pinning_lattice_ax, pinning_lattice_ay;
	double  pinning_xshift,	pinning_yshift;
	double  pinning_lx2, pinning_ly2;

	double	pinning_substrate_height;		//3000 nm
	double 	pinning_k_factor;				//a factor for the pinning spring constant
	double 	pinning_k;						//spring constant in pN/um

	double 	pinning_middle_k;				//middle bump spring constant / h in pN/um^2
	double	pinning_middle_height;			//middle bump height in um
	double 	pinning_middle_height_spread;	//spread arpund it

	double  particle_radius;				//about 10 microns
	double  particle_volume;				//about 4200 um^3
	double  particle_density;				//about 1.9
	double  particle_apparent_density;		//about 0.9
	double  particle_mass;					//about 8000 pg
	double  particle_apparent_weight;		//about 37 pN
	
	double  particle_magnetization;				//0.637 um^2*A
	double  particle_particle_force_pre_term;	// pN/A^2
	double  particle_particle_force_term;		// pN*um^4

	char 	particle_type[255];
	char 	particle_lattice_type[255];
	char 	particle_starting_position[255];	//random, biased or ground state
	char 	particle_defect_type[255];
	int 	particle_defect_number;
	int 	*defect_pinningsite_IDs;

	char 	interaction_type[255];
	double	multiplicator;
	double 	particle_q;
	double 	particle_a0;
	double 	r_cutoff_short;
	double 	r_cutoff_short2;	
	double 	r_cutoff_long;
	double 	r_cutoff_long2;	
	double 	r_cutoff_Verlet;
	double 	r_cutoff_Verlet2;

	char 	Verlet_rebuild_type[255];		//never; in the future: when moved a given dr2

	};

struct pinningsite_struct
	{
	double 		x,y;                //position of the center of the pinning site
	double 		cosfi,sinfi;        //fi is the angle with the horizontal
	double 		lx,ly;              //half length in x,y
	double 		middle_height;      //middle bump height
	int    		N_particles_in;     //number of particles inside the pinning site
	int    		particle_in_ID[2];	//so far max 2 particles per pinning site
	int    		groundstate_dir;	//+1 or -1 depending

    //special lattice related
    //the pinning site knows if it has been decimated
    int         decimated;          //0-filled 1-empty (Cristiano lattices)
    
	//spin related
	int spin_value;
	int spin_old_value;
	int spin_N_hops;
	double 	spin_distance_from_dd;
	int 	spin_distance_from_dd_bin;
	};

struct particle_struct
	{
	double 		x,y;
	double 		fx,fy;	
	double 		q;
	int 		pinningsite_ID;
    
    double      E; //energy of the particle due to the pp
    
	//each particle belongs to only one pinningsite max						
	//this might change in the future
	int 		color;
	//uint32_t	Z_index;
    int wasclosetosomeone;
	};

struct vertex_struct
	{
	double x,y;
	int type;
	int N_particles_in;
	int particle_in_ID[8];
	int chess_table_color; 	    // this is for square spin ice so that I can diff the 2 GS states
	int z;                      //this is for the decimated lattice, z can be 3 or 4
   
    	int N_neighbor_vertices;    //how many neighboring vertices it has 3 or 4
    	int N_neighbor_z3;          //how many of these neighbors are z3
  	int N_neighbor_z4;          //how many of these neighbors are z4
  	int neighbor_vertex_ID[4];  //which one are the neighboring vertices
    	int q;                      //this is the topological charge on the vertex
    	int q_neighbor_z3;          //total charge on the z3 neighbors
    	int q_neighbor_z4;          //total charge on the z4 neighbors
    	};

struct variables_struct  	
	{
	int N_pinningsites;
	int N_particles;
	int N_vertices;
	int N_pairs;
    
	//decimated project
	int N_vertices3;
    	int N_vertices4;
    	double eta_vertices;
    	double xi_vertices;
    
	struct pinningsite_struct	*pinningsites;
	struct particle_struct		*particles;
	struct vertex_struct		*vertices;
	int *Verlet_list_i;
	int *Verlet_list_j;

	double max_movement_x,max_movement_y;
	int N_tabulated;
	double *f_tabulated;
	double step_tabulated;
	};

struct statistics_struct
	{
	long int N_average;
	long int N_vertextype[7];							//maximum 7 types of vertices in the system
    	long int N_vertextype3[7];							//maximum 7 types of vertices in the system
	
	//this is for the decimated project
	long int N3,N4; 
	double xi, eta;   

	int q33,q34,q43,q44;	//topological charge on the neighbors
	

    	double pos1x,pos1y,pos3x,pos3y;
    	double total_energy;
	};

extern struct parameters_struct 	global_parameters;
extern struct variables_struct  	global_variables;
extern struct statistics_struct		global_statistics;
