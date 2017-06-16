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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "random_numrec.h"
#include "global_variables.h"
#include "initialization.h"

/* ==============================================================================================
void initialize_diagnostics_file()
Initializes a diagnostics file. This file will contain all the parameters read in and all the
echo that was sent out to console, to help later identify all parameters pertaining to a run
============================================================================================== */
void initialize_diagnostics_file()
{

strcpy(global_parameters.diagnostics_file_name,"results/diag_");
strcat(global_parameters.diagnostics_file_name,global_parameters.running_name);
strcat(global_parameters.diagnostics_file_name,".txt");

global_parameters.diagnostics_file = fopen(global_parameters.diagnostics_file_name,"wt");
if (global_parameters.diagnostics_file ==NULL)
	{
	printf("Could not create diagnostics file %s\n",global_parameters.diagnostics_file_name);
	exit(1);
	}
printf("Created diagnostics file %s\n",global_parameters.diagnostics_file_name);
global_parameters.diagnostics_file_ready = 1;
//if (global_parameters.diagnostics_file_ready) printf("%d\n",global_parameters.diagnostics_file_ready);

fprintf(global_parameters.diagnostics_file,"This is a diagnostics file for the colloidal code\n");
fprintf(global_parameters.diagnostics_file,"name running %s\n",global_parameters.running_name);
fprintf(global_parameters.diagnostics_file,"name diagnostics %s\n",global_parameters.diagnostics_file_name);

}


/* ==============================================================================================
void make_copy_of_parameter_file()
makes a copy of the parameter file with the name param_<runname>
============================================================================================== */
void make_copy_of_parameter_file(int argc, char *argv[])
{
FILE *copyfile;
char copyfile_name[255];
char copychar;
int number_of_variables_read;

if (argc<2)
	{
	strcpy(global_parameters.parameter_file_name,"parameters/parameters.txt");
	global_parameters.parameter_file = fopen(global_parameters.parameter_file_name,"rt");
	if (global_parameters.parameter_file==NULL)
		{
		printf("Could not find standard parameter file at parameters/parameter.txt");
		exit(1);
		}
	}
else
	{
	strcpy(global_parameters.parameter_file_name,argv[1]);
	global_parameters.parameter_file = fopen(global_parameters.parameter_file_name,"rt");
	if (global_parameters.parameter_file==NULL)
		{
		printf("Could not find specified parameter file at %s",global_parameters.parameter_file_name);
		exit(1);
		}
	}

printf("Re-opened parameter file %s for copying\n",global_parameters.parameter_file_name);
fprintf(global_parameters.diagnostics_file,"Re-opened parameter file %s for copying\n",global_parameters.parameter_file_name);

strcpy(copyfile_name,"results/param_");
strcat(copyfile_name,global_parameters.running_name);
strcat(copyfile_name,".txt");

copyfile = fopen(copyfile_name,"wt");
if (copyfile ==NULL)
	{
	printf("Could not create copy of the parameter file %s\n",copyfile_name);
	fprintf(global_parameters.diagnostics_file,"Could not create copy of the parameter file %s\n",copyfile_name);
	exit(1);
	}

while (!feof(global_parameters.parameter_file))
	{
	number_of_variables_read = fscanf(global_parameters.parameter_file,"%c",&copychar);
	//printf("%c",copychar);
	fprintf(copyfile,"%c",copychar);
	}

fclose(copyfile);
fclose(global_parameters.parameter_file);
printf("Copy complete %s\n",copyfile_name);
fprintf(global_parameters.diagnostics_file,"Copy complete %s\n",copyfile_name);
}


/* ==============================================================================================
read_parameters_and_variables(int argc, char *argv[])

reads in the parameters from the initialization or restart file
uses the line-by-line reading providede in the read_parameter_line()
============================================================================================== */
void read_parameters_and_variables(int argc, char *argv[])
{
int i;
int dummy_int;
double dummy_double;
char dummy_string[255];

if (strcmp(argv[1],"--help")==0)
	{
	printf("Colloidal spin ice simulation program\n");
	printf("Usage: colloidal [parameter file] [runname] [seed] [B] [Bias] [h] [h_spread]\n");
	exit(1);
	}

printf("==============================================================================================\n");
printf("Colloidal BD simulation\n");

printf("==============================================================================================\n");
printf("Parameters meaning\n");
if (argc<2) printf("[parameter file name] = by default parameters/parameters.txt\n");
else printf("[parameter file name] = %s\n",argv[1]);

if (argc>2) printf("[new_running_name] = %s\n", argv[2]);
if (argc>3) printf("[new_seed] = %s\n", argv[3]);
if (argc>4) printf("[new_magnetic_field] B = %s mT\n", argv[4]);
if (argc>5) printf("[new_biasing_force] Fx,y = %s\n", argv[5]);
if (argc>6) printf("[new_middle_height] = %s\n", argv[5]);
if (argc>7) printf("[new_middle_height_spread] = %s\n", argv[6]);


printf("==============================================================================================\n");
printf("Initialization\n");

if (argc<2)
	{
	strcpy(global_parameters.parameter_file_name,"parameters/parameters.txt");
	global_parameters.parameter_file = fopen(global_parameters.parameter_file_name,"rt");
	if (global_parameters.parameter_file==NULL)
		{
		printf("Could not find standard parameter file at parameters/parameter.txt");
		exit(1);
		}
	}
  	else
		{
		strcpy(global_parameters.parameter_file_name,argv[1]);
		global_parameters.parameter_file = fopen(global_parameters.parameter_file_name,"rt");
		if (global_parameters.parameter_file==NULL)
			{
			printf("Could not find specified parameter file at %s",global_parameters.parameter_file_name);
			exit(1);
			}
		}

printf("Opened parameter file %s\n",global_parameters.parameter_file_name);

global_parameters.diagnostics_file_ready = 0;

read_parameter_line("name","running","string","do not convert",&dummy_int,&dummy_double,global_parameters.running_name);
//ability to modify runname from command line
if (argc>2)
	{
	strcpy(global_parameters.running_name,argv[2]);
	printf("Command line overwrites runnig name to = %s\n",global_parameters.running_name);
	}

initialize_diagnostics_file();
read_parameter_line("name","movie","string","convert",&dummy_int,&dummy_double,global_parameters.movie_file_name);
read_parameter_line("name","statistics","string","convert",&dummy_int,&dummy_double,global_parameters.statistics_file_name);
read_parameter_line("name","savefile","string","convert",&dummy_int,&dummy_double,global_parameters.save_file_name);

read_parameter_line("init","type","int","convert",&global_parameters.initialization_type,&dummy_double,dummy_string);

read_parameter_line("seed","modifier","int","do not convert",&global_parameters.seedmodifier,&dummy_double,dummy_string);


//ability t modify seed from command line
if (argc>3)
	{

    if (strcmp(argv[3],"same")==0)
        {
            printf("Command line keeps the same random seed\n");
            printf("Random Seed = %d\n",global_parameters.seedmodifier);
        }
    else
        {
            global_parameters.seedmodifier = atoi(argv[3]);
            printf("Command line overwrites seed modifier to = %d (%d decimals)\n",global_parameters.seedmodifier,(int)floor(log10(global_parameters.seedmodifier)));
            printf("Seed size (int) is %ld bytes, max number is %.0lf (%d decimals recommended)\n",sizeof(int),pow(2,8*sizeof(int))/2,(int)floor(log10(pow(2,8*sizeof(int))/2))-1);
        }
    }


read_parameter_line("timing","total_runtime","int","do not convert",&global_parameters.total_runtime,&dummy_double,dummy_string);
read_parameter_line("timing","movie_write_time","int","do not convert",&global_parameters.movie_write_time,&dummy_double,dummy_string);
read_parameter_line("timing","echotime","int","do not convert",&global_parameters.echo_write_time,&dummy_double,dummy_string);
read_parameter_line("timing","statistics","int","do not convert",&global_parameters.statistics_write_time,&dummy_double,dummy_string);

//temperature related reading
read_parameter_line("temperature","initial","double","do not convert",&dummy_int,&global_parameters.temperature,dummy_string);
read_parameter_line("temperature","N_setpoints","int","do not convert",&global_parameters.N_temperature_set_points,&dummy_double,dummy_string);
global_parameters.set_temperature_time[0] = 0;
global_parameters.set_temperature[0] = global_parameters.temperature;

for(i=0;i<global_parameters.N_temperature_set_points;i++)
	{
	printf("Temperature setpoint # %d\n",i);
	fprintf(global_parameters.diagnostics_file,"Temperature setpoint # %d\n",i);
	read_parameter_line("temperature","setpoint_time","int","do not convert",&global_parameters.set_temperature_time[i+1],&dummy_double,dummy_string);
	read_parameter_line("temperature","setpoint_temp","double","do not convert",&dummy_int,&global_parameters.set_temperature[i+1],dummy_string);
	}

global_parameters.N_temperature_set_points++;

if (global_parameters.set_temperature_time[global_parameters.N_temperature_set_points-1]!= global_parameters.total_runtime)
	{
	global_parameters.set_temperature_time[global_parameters.N_temperature_set_points] = global_parameters.total_runtime;
	global_parameters.set_temperature[global_parameters.N_temperature_set_points] = global_parameters.set_temperature[global_parameters.N_temperature_set_points-1];
	global_parameters.N_temperature_set_points++;
	}

//external force related reading
read_parameter_line("external","initial_fx","double","do not convert",&dummy_int,&global_parameters.external_fx,dummy_string);
read_parameter_line("external","initial_fy","double","do not convert",&dummy_int,&global_parameters.external_fy,dummy_string);


if (argc>5)
    {
        //printf("%s\n",argv[4]);
        global_parameters.external_fx = atof(argv[5]);
        global_parameters.external_fy = global_parameters.external_fx;
        printf("Command line overwrites external force toto = %lf\n",global_parameters.external_fx);
    }


read_parameter_line("external","N_setpoints","int","do not convert",&global_parameters.N_external_set_points,&dummy_double,dummy_string);
//even if there are no setpoints the initial fx, fy are still stored as the first setpoint
//this is for runs with a constant bias
global_parameters.set_external_time[0] = 0;
global_parameters.set_external_fx[0] = global_parameters.external_fx;
global_parameters.set_external_fy[0] = global_parameters.external_fy;

for(i=0;i<global_parameters.N_external_set_points;i++)
	{
	printf("External force setpoint # %d\n",i);
	fprintf(global_parameters.diagnostics_file,"External force setpoint # %d\n",i);
	read_parameter_line("external","setpoint_time","int","do not convert",&global_parameters.set_external_time[i+1],&dummy_double,dummy_string);
	read_parameter_line("external","setpoint_fx","double","do not convert",&dummy_int,&global_parameters.set_external_fx[i+1],dummy_string);
	read_parameter_line("external","setpoint_fy","double","do not convert",&dummy_int,&global_parameters.set_external_fy[i+1],dummy_string);
	}
global_parameters.N_external_set_points++;
if (global_parameters.set_external_time[global_parameters.N_external_set_points-1]!= global_parameters.total_runtime)
	{
	global_parameters.set_external_time[global_parameters.N_external_set_points] = global_parameters.total_runtime;
	global_parameters.set_external_fx[global_parameters.N_external_set_points] = global_parameters.set_external_fx[global_parameters.N_external_set_points-1];
	global_parameters.set_external_fy[global_parameters.N_external_set_points] = global_parameters.set_external_fy[global_parameters.N_external_set_points-1];
	global_parameters.N_external_set_points++;
	}

//real units in the simulation
read_parameter_line("running","dt","double","do not convert",&dummy_int,&global_parameters.dt,dummy_string);
read_parameter_line("running","diffusion_const","double","do not convert",&dummy_int,&global_parameters.diffusion_constant,dummy_string);
//calculate derivative real units here from diffusion
global_parameters.boltzmann_constant = 0.0138064852;
printf("running kB (Boltzmann constant) = %lf pN*nm/K\n",global_parameters.boltzmann_constant);
fprintf(global_parameters.diagnostics_file,"kB (Boltzmann constant) = %lf pN*nm/K\n",global_parameters.boltzmann_constant);
global_parameters.kbT = global_parameters.boltzmann_constant * (273.15 + global_parameters.temperature); //temperature is in C
printf("running kBT (Boltzmann constant*T) = %lf pN*nm\n",global_parameters.kbT);
fprintf(global_parameters.diagnostics_file,"kBT (Boltzmann constant*T) = %lf pN*nm\n",global_parameters.kbT);
global_parameters.mobility = global_parameters.diffusion_constant / global_parameters.kbT / 1000.0;
printf("running mobility = %lf um/pN/s\n",global_parameters.mobility);
fprintf(global_parameters.diagnostics_file,"running mobility = %lf um/pN/s\n",global_parameters.mobility);
global_parameters.thermal_force_term = sqrt ( 2.0 / (global_parameters.diffusion_constant * global_parameters.dt) ) * global_parameters.kbT;
printf("running thermal force term = %lf pN\n",global_parameters.thermal_force_term);
fprintf(global_parameters.diagnostics_file,"running thermal force term = %lf pN\n",global_parameters.thermal_force_term);

//no longer calculating it from the diffusion constant
//global_parameters.magnetic_susceptibility = sqrt(79.5 / global_parameters.diffusion_constant );
//global_parameters.magnetic_susceptibility = 0.061;

global_parameters.magnetic_permeability = 4 * 3.14159265 * 1e5;
printf("running magnetic permeability = %lf pN/A^2\n",global_parameters.magnetic_permeability);
fprintf(global_parameters.diagnostics_file,"running magnetic permeability = %lf pN/A^2\n",global_parameters.magnetic_permeability);
read_parameter_line("running","magnetic_susceptibility","double","do not convert",&dummy_int,&global_parameters.magnetic_susceptibility,dummy_string);
//printf("running magnetic susceptibity = %lf [non-dimensional]\n",global_parameters.magnetic_susceptibility);
//fprintf(global_parameters.diagnostics_file,"running magnetic susceptibity = %lf [non-dimensional]\n",global_parameters.magnetic_susceptibility);
read_parameter_line("running","magnetic_field","double","do not convert",&dummy_int,&global_parameters.magnetic_field,dummy_string);
//ability to change the magnetic field from command line
if (argc>4)
	{
	//printf("%s\n",argv[4]);
	global_parameters.magnetic_field = atof(argv[4]);
	printf("Command line overwrites magnetic field to = %lf\n",global_parameters.magnetic_field);
	}

read_parameter_line("box","SX","double","convert",&dummy_int,&global_parameters.SX,dummy_string);
read_parameter_line("box","SY","double","convert",&dummy_int,&global_parameters.SY,dummy_string);
read_parameter_line("box","boundaryx","int","convert",&global_parameters.boundarycondition_x,&dummy_double,dummy_string);
read_parameter_line("box","boundaryy","int","convert",&global_parameters.boundarycondition_y,&dummy_double,dummy_string);

read_parameter_line("pinning","type","string","do not convert",&dummy_int,&dummy_double,global_parameters.pinning_setup_type);
read_parameter_line("pinning","lattice","string","do not convert",&dummy_int,&dummy_double,global_parameters.pinning_lattice_type);
read_parameter_line("pinning","lattice_Nx","int","do not convert",&global_parameters.pinning_lattice_Nx,&dummy_double,dummy_string);
read_parameter_line("pinning","lattice_Ny","int","do not convert",&global_parameters.pinning_lattice_Ny,&dummy_double,dummy_string);
read_parameter_line("pinning","lattice_ax","double","do not convert",&dummy_int,&global_parameters.pinning_lattice_ax,dummy_string);
read_parameter_line("pinning","lattice_ay","double","do not convert",&dummy_int,&global_parameters.pinning_lattice_ay,dummy_string);
read_parameter_line("pinning","xshift","double","do not convert",&dummy_int,&global_parameters.pinning_xshift,dummy_string);
read_parameter_line("pinning","yshift","double","do not convert",&dummy_int,&global_parameters.pinning_yshift,dummy_string);
read_parameter_line("pinning","lx2","double","do not convert",&dummy_int,&global_parameters.pinning_lx2,dummy_string);
read_parameter_line("pinning","ly2","double","do not convert",&dummy_int,&global_parameters.pinning_ly2,dummy_string);
read_parameter_line("pinning","substrate_height","double","do not convert",&dummy_int,&global_parameters.pinning_substrate_height,dummy_string);
//calculate the spring constant factor
global_parameters.pinning_k_factor = global_parameters.pinning_substrate_height / 5000.0 / 5000.0; // in 1/nm units
printf("pinning spring constant factor = %lf 1/nm\n",global_parameters.pinning_k_factor);
fprintf(global_parameters.diagnostics_file,"pinning spring constant factor = %lf 1/nm\n",global_parameters.pinning_k_factor);
//calculate the middle bump height
//global_parameters.pinning_middle_height = 31800 / global_parameters.diffusion_constant;
//global_parameters.pinning_middle_height = 0.870;
read_parameter_line("pinning","middle_height","double","do not convert",&dummy_int,&global_parameters.pinning_middle_height,dummy_string);
//printf("pinning middle height = %lf um\n",global_parameters.pinning_middle_height);
//fprintf(global_parameters.diagnostics_file,"pinning middle height = %lf um\n",global_parameters.pinning_middle_height);

if (argc>6)
	{
	global_parameters.pinning_middle_height = atof(argv[6]);
	printf("Command line overwrites pinning middle height to = %lf\n",global_parameters.pinning_middle_height);
	}

read_parameter_line("pinning","middle_h_spread","double","do not convert",&dummy_int,&global_parameters.pinning_middle_height_spread,dummy_string);
if (argc>7)
	{
	global_parameters.pinning_middle_height_spread = atof(argv[7]);
	printf("Command line overwrites pinning middle height spread to = %lf\n",global_parameters.pinning_middle_height_spread);
	}

//printf("pinning middle height spread = %lf um\n",global_parameters.pinning_middle_height_spread);
//fprintf(global_parameters.diagnostics_file,"pinning middle height spread = %lf um\n",global_parameters.pinning_middle_height_spread);
read_parameter_line("particle","radius","double","do not convert",&dummy_int,&global_parameters.particle_radius,dummy_string);
//calculate global particle parameters (same for all particles)
global_parameters.particle_volume = 4*3.141592/3.0*global_parameters.particle_radius*global_parameters.particle_radius*global_parameters.particle_radius;
printf("particle volume = %lf um^3\n",global_parameters.particle_volume);
fprintf(global_parameters.diagnostics_file,"particle volume = %lf [um^3]\n",global_parameters.particle_volume);
read_parameter_line("particle","density","double","do not convert",&dummy_int,&global_parameters.particle_density,dummy_string);
global_parameters.particle_apparent_density = global_parameters.particle_density - 1.0;
printf("particle apparent density = %lf g/cm^3\n",global_parameters.particle_apparent_density);
fprintf(global_parameters.diagnostics_file,"particle apparent density = %lf g/cm^3\n",global_parameters.particle_apparent_density);
global_parameters.particle_mass = global_parameters.particle_volume * global_parameters.particle_density;
printf("particle mass = %lf pg\n",global_parameters.particle_mass);
fprintf(global_parameters.diagnostics_file,"particle mass = %lf pg\n",global_parameters.particle_mass);
global_parameters.particle_apparent_weight = global_parameters.particle_volume * global_parameters.particle_apparent_density * 9.81 / 1000.0;
//divided by 1000 to make it into pN units of force it was um^3*g/cm^3*m/s^2 which is 0.001 pN
printf("particle apparent weight = %lf pN\n",global_parameters.particle_apparent_weight);
fprintf(global_parameters.diagnostics_file,"particle apparent weight = %lf pN\n",global_parameters.particle_apparent_weight);

global_parameters.particle_magnetization = 1000.0 * global_parameters.magnetic_susceptibility * global_parameters.particle_volume * global_parameters.magnetic_field / global_parameters.magnetic_permeability;
printf("particle magnetization = %lf um^2*A\n",global_parameters.particle_magnetization);
fprintf(global_parameters.diagnostics_file,"particle magnetization = %lf um^2*A\n",global_parameters.particle_magnetization);

global_parameters.particle_particle_force_pre_term = 3.0 * global_parameters.magnetic_permeability / 2.0 / 3.14159265;
printf("particle particle force pre term = %lf pN/A^2\n",global_parameters.particle_particle_force_pre_term);
fprintf(global_parameters.diagnostics_file,"particle particle force pre term = %lf pN/A^2\n",global_parameters.particle_particle_force_pre_term);


global_parameters.particle_particle_force_term = global_parameters.particle_particle_force_pre_term * global_parameters.particle_magnetization * global_parameters.particle_magnetization;


printf("particle particle force term = %lf pN*um^4\n",global_parameters.particle_particle_force_term);
fprintf(global_parameters.diagnostics_file,"particle particle force term = %lf pN*um^4\n",global_parameters.particle_particle_force_term);

printf("estimated max p-p force at 10um = %lf pN\n",global_parameters.particle_particle_force_term/10.0/10/10/10);
fprintf(global_parameters.diagnostics_file,"estimated max p-p force at 10um = %lf pN\n",global_parameters.particle_particle_force_term/10.0/10/10/10);
printf("estimated max p-p force at 20um = %lf pN\n",global_parameters.particle_particle_force_term/20.0/20/20/20);
fprintf(global_parameters.diagnostics_file,"estimated max p-p force at 20um = %lf pN\n",global_parameters.particle_particle_force_term/20.0/20/20/20);
printf("estimated max p-p force at 30um = %lf pN\n",global_parameters.particle_particle_force_term/30.0/30/30/30);
fprintf(global_parameters.diagnostics_file,"estimated max p-p force at 30um = %lf pN\n",global_parameters.particle_particle_force_term/30.0/30/30/30);

printf("pinning related calculations based on particle weight\n");
fprintf(global_parameters.diagnostics_file,"pinning related calculations based on particle weight\n");
global_parameters.pinning_k = 2000.0 * global_parameters.particle_apparent_weight * global_parameters.pinning_k_factor;
printf("pinning spring constant = %lf pN/um\n",global_parameters.pinning_k);
fprintf(global_parameters.diagnostics_file,"pinning spring constant = %lf pN/um\n",global_parameters.pinning_k);
global_parameters.pinning_middle_k = global_parameters.particle_apparent_weight * 2.0 / global_parameters.pinning_lx2 / global_parameters.pinning_lx2;
printf("pinning middle spring constant / middle height = %lf pN/um^2\n",global_parameters.pinning_middle_k);
fprintf(global_parameters.diagnostics_file,"pinning middle spring constant / middle height = %lf pN/um^2\n",global_parameters.pinning_middle_k);
printf("expected pinning middle force max = %lf pN\n", global_parameters.pinning_middle_k * global_parameters.pinning_lx2 * global_parameters.pinning_middle_height);
fprintf(global_parameters.diagnostics_file,"expected pinning middle force max = %lf pN\n", global_parameters.pinning_middle_k * global_parameters.pinning_lx2 * global_parameters.pinning_middle_height);

read_parameter_line("particle","type","string","do not convert",&dummy_int,&dummy_double,global_parameters.particle_type);
read_parameter_line("particle","lattice","string","do not convert",&dummy_int,&dummy_double,global_parameters.particle_lattice_type);
read_parameter_line("particle","starting","string","do not convert",&dummy_int,&dummy_double,global_parameters.particle_starting_position);
read_parameter_line("particle","defect_type","string","do not convert",&dummy_int,&dummy_double,global_parameters.particle_defect_type);
read_parameter_line("particle","defect_number","int","do not convert",&global_parameters.particle_defect_number,&dummy_double,dummy_string);
read_parameter_line("particle","interaction","string","do not convert",&dummy_int,&dummy_double,global_parameters.interaction_type);
read_parameter_line("particle","multiplicator","double","do not convert",&dummy_int,&global_parameters.multiplicator,dummy_string);
read_parameter_line("particle","charge","double","do not convert",&dummy_int,&global_parameters.particle_q,dummy_string);
read_parameter_line("particle","screening_a","double","do not convert",&dummy_int,&global_parameters.particle_a0,dummy_string);
read_parameter_line("particle","short_cutoff","double","do not convert",&dummy_int,&global_parameters.r_cutoff_short,dummy_string);
global_parameters.r_cutoff_short2 = global_parameters.r_cutoff_short * global_parameters.r_cutoff_short;
printf("particle short_cutoff^2 %lf\n",global_parameters.r_cutoff_short2);
fprintf(global_parameters.diagnostics_file,"particle short_cutoff^2 %lf\n",global_parameters.r_cutoff_short2);
read_parameter_line("particle","long_cutoff","double","do not convert",&dummy_int,&global_parameters.r_cutoff_long,dummy_string);
global_parameters.r_cutoff_long2 = global_parameters.r_cutoff_long * global_parameters.r_cutoff_long;
printf("particle long_cutoff^2 %lf\n",global_parameters.r_cutoff_long2);
fprintf(global_parameters.diagnostics_file,"particle long_cutoff^2 %lf\n",global_parameters.r_cutoff_long2);
read_parameter_line("particle","Verlet_cutoff","double","do not convert",&dummy_int,&global_parameters.r_cutoff_Verlet,dummy_string);
global_parameters.r_cutoff_Verlet2 = global_parameters.r_cutoff_Verlet * global_parameters.r_cutoff_Verlet;
printf("particle verlet_cutoff^2 %lf\n",global_parameters.r_cutoff_Verlet2);
fprintf(global_parameters.diagnostics_file,"particle verlet_cutoff^2 %lf\n",global_parameters.r_cutoff_Verlet2);
read_parameter_line("particle","Verlet_rebuild","string","do not convert",&dummy_int,&dummy_double,global_parameters.Verlet_rebuild_type);

global_parameters.time = 0;
//zero magnetic field
global_parameters.initial_magnetic_field = global_parameters.magnetic_field;
global_parameters.magnetic_field = 0.0;

printf("Initializing random variable generator\n");
fprintf(global_parameters.diagnostics_file,"Initializing random variable generator\n");
printf("Initialization complete\n");
fprintf(global_parameters.diagnostics_file,"Initialization complete\n");
fflush(global_parameters.diagnostics_file);

fclose(global_parameters.parameter_file);
}

/* ==============================================================================================
void read_parameter_line(char expected_type_specifier[], char expected_subtype_specifier[],
char expected_variable_type[], char conversion_yes_no[], int *target_int, double *target_double,
char target_string[255])

reads in one line from the parameter file
never checks for the expected type specifier!!!!
============================================================================================== */
void read_parameter_line(char expected_type_specifier[], char expected_subtype_specifier[],
	char expected_variable_type[], char conversion_yes_no[], int *target_int, double *target_double,
	char target_string[255])
{

char type_specifier[255];
char subtype_specifier[255];
char variable_value[255];

char dummy_string[255];
char *dummy_string_pointer;

int number_of_variables_read;
int is_a_comment_line;

int  int_variable;
double double_variable;
char string_variable[255];

char *charpointer;
int charposition;
int i;
char remainder_string[255];

is_a_comment_line = 1;

do
	{
	number_of_variables_read = fscanf(global_parameters.parameter_file,"%s",type_specifier);
	if (number_of_variables_read==1)
		{
		if (strcmp(type_specifier,"#")==0)
			{
			dummy_string_pointer = fgets(dummy_string,255,global_parameters.parameter_file);
			}
		else
			{
			is_a_comment_line = 0;
			}
		}
	}
while (is_a_comment_line==1);

if (strcmp(type_specifier,expected_type_specifier)!=0)
	{
	printf("I was expecting a different type specifier\n");
	printf("I was expecting: %s\n",expected_type_specifier);
	printf("I found: %s\n",type_specifier);
	exit(1);
	}


number_of_variables_read = fscanf(global_parameters.parameter_file,"%s",subtype_specifier);
if (number_of_variables_read!=1)
	{
	printf("I cannot find the subtype specifier\n");
	exit(1);
	}

if (strcmp(subtype_specifier,expected_subtype_specifier)!=0)
	{
	printf("I was expecting a different subtype specifier\n");
	printf("I was expecting: %s\n",expected_subtype_specifier);
	printf("I found: %s\n",subtype_specifier);
	exit(1);
	}

number_of_variables_read = fscanf(global_parameters.parameter_file,"%s",variable_value);
if (number_of_variables_read!=1)
	{
	printf("I cannot find the value\n");
	exit(1);
	}

if (strcmp(conversion_yes_no,"convert")==0)
	{

	//init clean or restart conversion to 0 or 1
	if (( strcmp(type_specifier,"init")==0) && (strcmp(subtype_specifier,"type")==0))
		{
		if (strcmp(variable_value,"clean")==0) int_variable = 0;
		if (strcmp(variable_value,"restart")==0) int_variable = 1;
		}
	//insert the base name instead of <running> to make it easy to change
	if (( strcmp(type_specifier,"name")==0))
		{
		/*<running> 9 chars is replaced by the running name*/
		charpointer = strstr(variable_value,"<running>");
		if (charpointer!=NULL) charposition = charpointer - variable_value;

		for(i=0;i<charposition;i++)
			string_variable[i] = variable_value[i];
		string_variable[charposition] = '\0';
		for(i=charposition+strlen("<running>");i<strlen(variable_value);i++)
			remainder_string[i-charposition-strlen("<running>")] = variable_value[i];
		remainder_string[strlen(variable_value)-strlen("<running>")-charposition] = '\0';
		strcat(string_variable,global_parameters.running_name);
		strcat(string_variable,remainder_string);
		}
	//SX SY if it needs to be recalculated
	if (( strcmp(type_specifier,"box")==0) && (strcmp(subtype_specifier,"SX")==0))
		{

		if (strcmp(variable_value,"calculate_spin_ice")==0)
			{
			global_parameters.recalculate_SX = 1;
			double_variable = 0.0;
			}
		else
			{
			global_parameters.recalculate_SX = 0;
			double_variable = atof(variable_value);
			}
		}
	if (( strcmp(type_specifier,"box")==0) && (strcmp(subtype_specifier,"SY")==0))
		{
		if (strcmp(variable_value,"calculate_spin_ice")==0)
			{
			global_parameters.recalculate_SY = 1;
			double_variable = 0.0;
			}
		else
			{
			global_parameters.recalculate_SY = 0;
			double_variable = atof(variable_value);
			}
		}
    if (( strcmp(type_specifier,"box")==0) && (strcmp(subtype_specifier,"boundaryy")==0))
		{
		if (strcmp(variable_value,"NO_PBC")==0) int_variable = 0;
		if (strcmp(variable_value,"PBC")==0) 	int_variable = 1;
		}
    if (( strcmp(type_specifier,"box")==0) && (strcmp(subtype_specifier,"boundaryx")==0))
		{
		if (strcmp(variable_value,"NO_PBC")==0) int_variable = 0;
		if (strcmp(variable_value,"PBC")==0) 	int_variable = 1;
		}

	}
else
	{
	//convert an int value
	if (strcmp(expected_variable_type,"int")==0)
		int_variable = atoi(variable_value);
	if (strcmp(expected_variable_type,"double")==0)
		double_variable = atof(variable_value);
	if (strcmp(expected_variable_type,"string")==0)
		strcpy(string_variable,variable_value);
	}


if (strcmp(expected_variable_type,"int")==0)
	{
	*target_int = int_variable;
	printf("%s %s %d\n",type_specifier, subtype_specifier, int_variable);
	if (global_parameters.diagnostics_file_ready) { fprintf(global_parameters.diagnostics_file,"%s %s %d\n",type_specifier, subtype_specifier, int_variable);fflush(global_parameters.diagnostics_file);}

	}
	if (strcmp(expected_variable_type,"double")==0)
	{
	*target_double = double_variable;
	printf("%s %s %lf\n",type_specifier, subtype_specifier, double_variable);
	if (global_parameters.diagnostics_file_ready) { fprintf(global_parameters.diagnostics_file,"%s %s %lf\n",type_specifier, subtype_specifier, double_variable);fflush(global_parameters.diagnostics_file);}
	}
	if (strcmp(expected_variable_type,"string")==0)
	{
	strcpy(target_string,string_variable);
	printf("%s %s %s\n",type_specifier, subtype_specifier, string_variable);
	if (global_parameters.diagnostics_file_ready) { fprintf(global_parameters.diagnostics_file,"%s %s %s\n",type_specifier, subtype_specifier, string_variable);fflush(global_parameters.diagnostics_file);}
	}

}


/* ==============================================================================================
void initialize_pinning_sites()

initializes the pinning sites based on the pinning type specified
============================================================================================== */
void initialize_pinning_sites()
{

//printf("%s %s\n",global_parameters.pinning_setup_type,global_parameters.pinning_lattice_type);

if ( (strcmp(global_parameters.pinning_setup_type,"spin_ice_pin")==0)
		&&(strcmp(global_parameters.pinning_lattice_type,"square")==0) )
            {
            initialize_square_pinning();
            decimator();
            }

if ( (strcmp(global_parameters.pinning_setup_type,"spin_ice_pin")==0)
		&&(strcmp(global_parameters.pinning_lattice_type,"hexagonal")==0) ) initialize_hexagonal_pinning();

}

/* ==============================================================================================
void initialize_pinning_sites_spins()
initializes the pinning sites spin related information
this needs to be called separately from the original pin initialization
as the particles need to be already placed in the pins
============================================================================================== */
void initialize_pinning_sites_spins()
{
int i,j;
int i_particle;
double dx,dy;
double dr,dr2;
double dr_lattice;
double dxprime,dyprime;
double sinfi,cosfi;
int N_at_a_dd_distance[3];

for(i=0;i<global_variables.N_pinningsites;i++)
	{
	//this is an empty pinning site
	if (global_variables.pinningsites[i].N_particles_in==0)
		{
		global_variables.pinningsites[i].spin_value = 0;
		global_variables.pinningsites[i].spin_old_value = 0;
	 	global_variables.pinningsites[i].spin_N_hops = 0;
		}
	else if (global_variables.pinningsites[i].N_particles_in==2)
		{
		global_variables.pinningsites[i].spin_value = 2;
		global_variables.pinningsites[i].spin_old_value = 2;
	 	global_variables.pinningsites[i].spin_N_hops = 0;
		}
	else //should be 1 particle in the site no more no less
		{
		i_particle = global_variables.pinningsites[i].particle_in_ID[0];
		dx = global_variables.particles[i_particle].x - global_variables.pinningsites[i].x;
		dy = global_variables.particles[i_particle].y - global_variables.pinningsites[i].y;
		//PBC fold back
		if (dx>(global_parameters.SX/2.0)) 	dx = dx - global_parameters.SX;
		if (dx<(-global_parameters.SX/2.0)) dx = dx + global_parameters.SX;
		if (dy>(global_parameters.SY/2.0)) 	dy = dy - global_parameters.SY;
		if (dy<(-global_parameters.SY/2.0)) dy = dy + global_parameters.SY;

		sinfi = global_variables.pinningsites[i].sinfi;
		cosfi = global_variables.pinningsites[i].cosfi;

		//rotate all coordinates so the pin lies horizontally
		dxprime =  cosfi*dx + sinfi*dy;
		dyprime = -sinfi*dx + cosfi*dy;

		if (dxprime>=0)
			{
			global_variables.pinningsites[i].spin_value = 1;
			global_variables.pinningsites[i].spin_old_value = 1;
	 		global_variables.pinningsites[i].spin_N_hops = 0;

	 		//global_variables.particles[i_particle].color = 3;
	 		//printf("3");
			}
		else
			{
			global_variables.pinningsites[i].spin_value = -1;
			global_variables.pinningsites[i].spin_old_value = -1;
	 		global_variables.pinningsites[i].spin_N_hops = 0;

	 		//global_variables.particles[i_particle].color = 2;
	 		//printf("3");
			}
		}

	}

//find the distance of the pinning site to the nearest double defect
//also decide on which bin it's in (currently 2 bins: close=0 and far=1)

for(i=0;i<global_variables.N_pinningsites;i++)
	{
	//the largest possible distance is insterted here (x2) because we will be searching for a distance smaller than this
	global_variables.pinningsites[i].spin_distance_from_dd = sqrt(global_parameters.SX*global_parameters.SX + global_parameters.SY*global_parameters.SY);
	//bin -1 means it is not in a bin
	global_variables.pinningsites[i].spin_distance_from_dd_bin = -1;
	}


for(i=0;i<global_variables.N_pinningsites;i++)
	{
	//if this pinningsite is a double defect the distance is 0
	//but there is going to be no hopping counted
	//(by definition hopping is only counted for singly occupied)
	//there is no need to search all other pinningsites for a neighboring double defect
	if (global_variables.pinningsites[i].N_particles_in==2)
		{
		global_variables.pinningsites[i].spin_distance_from_dd = 0.0;
		global_variables.pinningsites[i].spin_distance_from_dd_bin = 0;
		}
	else
		for(j=0;j<global_variables.N_pinningsites;j++)
			if ((i!=j)&&(global_variables.pinningsites[j].N_particles_in==2))
				{
				dx = global_variables.pinningsites[j].x - global_variables.pinningsites[i].x;
				dy = global_variables.pinningsites[j].y - global_variables.pinningsites[i].y;

				//PBC fold back
				if (dx>(global_parameters.SX/2.0)) 	dx = dx - global_parameters.SX;
				if (dx<(-global_parameters.SX/2.0)) dx = dx + global_parameters.SX;
				if (dy>(global_parameters.SY/2.0)) 	dy = dy - global_parameters.SY;
				if (dy<(-global_parameters.SY/2.0)) dy = dy + global_parameters.SY;

				dr2 = dx*dx + dy*dy;
				dr = sqrt(dr2);
				//if we find one that is closer, this is the distance
				if (dr<global_variables.pinningsites[i].spin_distance_from_dd) global_variables.pinningsites[i].spin_distance_from_dd = dr;
				}
	}


//this is going to determine whether the particle is close or it's far from the double defect
//I take only the x direction lattice constant, if the lattice is elongated in y for some reason this needs to be rewritten
if ( (strcmp(global_parameters.pinning_setup_type,"spin_ice_pin")==0)
		&&(strcmp(global_parameters.pinning_lattice_type,"square")==0) )
	dr_lattice = 1.1*sqrt(global_parameters.pinning_lattice_ax*global_parameters.pinning_lattice_ax);
else if ( (strcmp(global_parameters.pinning_setup_type,"spin_ice_pin")==0)
		&&(strcmp(global_parameters.pinning_lattice_type,"hexagonal")==0) )
	dr_lattice = 0.9*sqrt(global_parameters.pinning_lattice_ax*global_parameters.pinning_lattice_ax);

printf("Distance to consider a spin close to a double defect %lf\n",dr_lattice);
fprintf(global_parameters.diagnostics_file,"Distance to consider a spin close to a double defect %lf\n",dr_lattice);

for(i=0;i<3;i++)
 N_at_a_dd_distance[i]=0;

for(i=0;i<global_variables.N_pinningsites;i++)
	{
	if (global_variables.pinningsites[i].spin_distance_from_dd<dr_lattice) global_variables.pinningsites[i].spin_distance_from_dd_bin = 1;
	else global_variables.pinningsites[i].spin_distance_from_dd_bin = 2;

	if (global_variables.pinningsites[i].N_particles_in==2) global_variables.pinningsites[i].spin_distance_from_dd_bin = 0;

	N_at_a_dd_distance[global_variables.pinningsites[i].spin_distance_from_dd_bin]++;
	}


for(i=0;i<3;i++)
	{
	printf("Pinningsites at dist = %d N  = %d",i,N_at_a_dd_distance[i]);if (i<2) printf("\n");
	fprintf(global_parameters.diagnostics_file,"Pinningsites at dist = %d N  = %d\n",i,N_at_a_dd_distance[i]);
	if (i<2) fprintf(global_parameters.diagnostics_file,"\n");
	}
printf("\n");
fprintf(global_parameters.diagnostics_file,"\n");

/*
FILE *testfile;

testfile = fopen("distance_0.txt","wt");
for(i=0;i<global_variables.N_pinningsites;i++)
	if (global_variables.pinningsites[i].spin_distance_from_dd_bin==0)
			fprintf(testfile,"%lf %lf\n",global_variables.pinningsites[i].x,global_variables.pinningsites[i].y);
fclose(testfile);
testfile = fopen("distance_1.txt","wt");
for(i=0;i<global_variables.N_pinningsites;i++)
	if (global_variables.pinningsites[i].spin_distance_from_dd_bin==1)
			fprintf(testfile,"%lf %lf\n",global_variables.pinningsites[i].x,global_variables.pinningsites[i].y);
fclose(testfile);
testfile = fopen("distance_2.txt","wt");
for(i=0;i<global_variables.N_pinningsites;i++)
	if (global_variables.pinningsites[i].spin_distance_from_dd_bin==2)
			fprintf(testfile,"%lf %lf\n",global_variables.pinningsites[i].x,global_variables.pinningsites[i].y);
fclose(testfile);
*/
fflush(global_parameters.diagnostics_file);
}


/* ==============================================================================================
void initialize_square_pinning()
============================================================================================== */
void initialize_square_pinning()
{
int i,j,k;

printf("==============================================================================================\n");
fprintf(global_parameters.diagnostics_file,"==============================================================================================\n");
printf("Initializing square ice pinning sites\n");
fprintf(global_parameters.diagnostics_file,"Initializing square ice pinning sites\n");

if (global_parameters.recalculate_SX==1) global_parameters.SX = global_parameters.pinning_lattice_Nx * global_parameters.pinning_lattice_ax;
if (global_parameters.recalculate_SY==1) global_parameters.SY = global_parameters.pinning_lattice_Ny * global_parameters.pinning_lattice_ay;

if ((global_parameters.recalculate_SX==1)||(global_parameters.recalculate_SY==1))
	{
	printf("Recalculate in x direction: %d x %lf = %lf\n",global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_ax,global_parameters.SX);
	fprintf(global_parameters.diagnostics_file,"Recalculate in x direction: %d x %lf = %lf\n",global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_ax,global_parameters.SX);
	printf("Recalculate in y direction: %d x %lf = %lf\n",global_parameters.pinning_lattice_Ny,global_parameters.pinning_lattice_ay,global_parameters.SY);
	fprintf(global_parameters.diagnostics_file,"Recalculate in y direction: %d x %lf = %lf\n",global_parameters.pinning_lattice_Ny,global_parameters.pinning_lattice_ay,global_parameters.SY);
	printf("Recalculated box size: SX = %lf SY = %lf\n",global_parameters.SX, global_parameters.SY);
	fprintf(global_parameters.diagnostics_file,"Recalculated box size: SX = %lf SY = %lf\n",global_parameters.SX, global_parameters.SY);
	}
global_variables.N_pinningsites = global_parameters.pinning_lattice_Nx * global_parameters.pinning_lattice_Ny*2;

printf("Number of pinning sites on a square %d %d lattice is %d %d x 2 = %d\n",
	global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_Ny,
	global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_Ny,
	global_variables.N_pinningsites);
fprintf(global_parameters.diagnostics_file,"Number of pinning sites on a square %d %d lattice is %d %d x 2 = %d\n",
	global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_Ny,
	global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_Ny,
	global_variables.N_pinningsites);

printf("Allocating space for %d pinning sites [%ld bytes]\n",global_variables.N_pinningsites,global_variables.N_pinningsites*sizeof(struct pinningsite_struct));
fprintf(global_parameters.diagnostics_file,"Allocating space for %d pinning sites [%ld bytes]\n",global_variables.N_pinningsites,global_variables.N_pinningsites*sizeof(struct pinningsite_struct));

global_variables.pinningsites = (struct pinningsite_struct *)malloc(global_variables.N_pinningsites*sizeof(struct pinningsite_struct));

k=0;
for(i=0;i<global_parameters.pinning_lattice_Nx;i++)
	for(j=0;j<global_parameters.pinning_lattice_Ny;j++)
		{
		//horizontal pinning site
		global_variables.pinningsites[k].x = (i + 0.5 + global_parameters.pinning_xshift) * global_parameters.pinning_lattice_ax;
		global_variables.pinningsites[k].y = (j + global_parameters.pinning_yshift) * global_parameters.pinning_lattice_ay;
		global_variables.pinningsites[k].lx = global_parameters.pinning_lx2;
		global_variables.pinningsites[k].ly = global_parameters.pinning_ly2;
		global_variables.pinningsites[k].cosfi = 1.0;
		global_variables.pinningsites[k].sinfi = 0.0;
        	global_variables.pinningsites[k].middle_height = global_parameters.pinning_middle_height + gasdev() * global_parameters.pinning_middle_height_spread;
        	if (global_variables.pinningsites[k].middle_height < 0) {global_variables.pinningsites[k].middle_height = 0.0;printf("Warning:middle_height clipped\n");}
		global_variables.pinningsites[k].N_particles_in = 0;
		global_variables.pinningsites[k].particle_in_ID[0] = 0;
		global_variables.pinningsites[k].particle_in_ID[1] = 0;
		if ( (i+j)%2==0 ) global_variables.pinningsites[k].groundstate_dir = -1;
		if ( (i+j)%2==1 ) global_variables.pinningsites[k].groundstate_dir = 1;
        global_variables.pinningsites[k].decimated = 0;
		//if (k==0)
		//printf("test: %lf %lf\n",global_variables.pinningsites[k].x,global_variables.pinningsites[k].y);
		k++;

		//vertical pinning site
		global_variables.pinningsites[k].x = (i+global_parameters.pinning_xshift) * global_parameters.pinning_lattice_ax;
		global_variables.pinningsites[k].y = (j+0.5+global_parameters.pinning_yshift) * global_parameters.pinning_lattice_ay;
		global_variables.pinningsites[k].lx = global_parameters.pinning_lx2;
		global_variables.pinningsites[k].ly = global_parameters.pinning_ly2;
		global_variables.pinningsites[k].cosfi = 0.0;
		global_variables.pinningsites[k].sinfi = 1.0;
		global_variables.pinningsites[k].middle_height = global_parameters.pinning_middle_height + gasdev() * global_parameters.pinning_middle_height_spread;
        	if (global_variables.pinningsites[k].middle_height < 0) {global_variables.pinningsites[k].middle_height = 0.0;printf("Warning:middle_height clipped\n");}
		global_variables.pinningsites[k].N_particles_in = 0;
		global_variables.pinningsites[k].particle_in_ID[0] = 0;
		global_variables.pinningsites[k].particle_in_ID[1] = 0;
		if ( (i+j)%2==0 ) global_variables.pinningsites[k].groundstate_dir = 1;
		if ( (i+j)%2==1 ) global_variables.pinningsites[k].groundstate_dir = -1;
        global_variables.pinningsites[k].decimated = 0;
		k++;
		}

fflush(global_parameters.diagnostics_file);
}

/* ==============================================================================================
void initialize_hexagonal_pinning()
============================================================================================== */
void initialize_hexagonal_pinning()
{
int i,j,k;

printf("==============================================================================================\n");
fprintf(global_parameters.diagnostics_file,"==============================================================================================\n");
printf("Initializing hexagonal ice pinning sites\n");
fprintf(global_parameters.diagnostics_file,"==============================================================================================\n");

if (global_parameters.pinning_lattice_ay!=sqrt(3)/2.0*global_parameters.pinning_lattice_ax)
	{
	global_parameters.pinning_lattice_ay=sqrt(3)/2.0*global_parameters.pinning_lattice_ax;
	printf("Fixed ay to be a regular hexagonal ax = %lf ay = %lf\n",
		global_parameters.pinning_lattice_ax,global_parameters.pinning_lattice_ay);
	fprintf(global_parameters.diagnostics_file,"Fixed ay to be a regular hexagonal ax = %lf ay = %lf\n",
		global_parameters.pinning_lattice_ax,global_parameters.pinning_lattice_ay);
	}

if (global_parameters.recalculate_SX==1) global_parameters.SX = global_parameters.pinning_lattice_Nx * global_parameters.pinning_lattice_ax * 3.0;
if (global_parameters.recalculate_SY==1) global_parameters.SY = global_parameters.pinning_lattice_Ny * global_parameters.pinning_lattice_ay * 2.0;

if ((global_parameters.recalculate_SX==1)||(global_parameters.recalculate_SY==1))
	{
	printf("Recalculate in x direction: %d x %lf = %lf\n",global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_ax,global_parameters.SX);
	fprintf(global_parameters.diagnostics_file,"Recalculate in x direction: %d x %lf = %lf\n",global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_ax,global_parameters.SX);
	printf("Recalculate in y direction: %d x %lf = %lf\n",global_parameters.pinning_lattice_Ny,global_parameters.pinning_lattice_ay,global_parameters.SY);
	fprintf(global_parameters.diagnostics_file,	"Recalculate in y direction: %d x %lf = %lf\n",global_parameters.pinning_lattice_Ny,global_parameters.pinning_lattice_ay,global_parameters.SY);
	printf("Recalculated box size: SX = %lf SY = %lf\n",global_parameters.SX, global_parameters.SY);
	fprintf(global_parameters.diagnostics_file,"Recalculated box size: SX = %lf SY = %lf\n",global_parameters.SX, global_parameters.SY);
	}

global_variables.N_pinningsites = global_parameters.pinning_lattice_Nx * global_parameters.pinning_lattice_Ny*6;

printf("Allocating space for %d pinning sites [%ld bytes]\n",global_variables. N_pinningsites,global_variables. N_pinningsites*sizeof(struct pinningsite_struct));
fprintf(global_parameters.diagnostics_file,"Allocating space for %d pinning sites [%ld bytes]\n",global_variables. N_pinningsites,global_variables. N_pinningsites*sizeof(struct pinningsite_struct));
global_variables.pinningsites = (struct pinningsite_struct *)malloc(global_variables. N_pinningsites*sizeof(struct pinningsite_struct));

printf("Number of pinning sites on a hexagonal %d %d lattice is %d %d x 6 = %d\n",
	global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_Ny,
	global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_Ny,
	global_variables.N_pinningsites);
fprintf(global_parameters.diagnostics_file,"Number of pinning sites on a hexagonal %d %d lattice is %d %d x 6 = %d\n",
	global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_Ny,
	global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_Ny,
	global_variables.N_pinningsites);


k=0;
for(i=0;i<global_parameters.pinning_lattice_Nx;i++)
	for(j=0;j<global_parameters.pinning_lattice_Ny;j++)
		{
		//first horizontal pinning site
		global_variables.pinningsites[k].x = ( (i*3) + 0.5 + global_parameters.pinning_xshift ) * global_parameters.pinning_lattice_ax;
		//printf(">>%lf %lf\n",(i*3) + 0.5 + global_parameters.pinning_xshift,global_variables.pinningsites[k].x);
		global_variables.pinningsites[k].y = ((j*2)+1+global_parameters.pinning_yshift) * global_parameters.pinning_lattice_ay;

		global_variables.pinningsites[k].lx = global_parameters.pinning_lx2;
		global_variables.pinningsites[k].ly = global_parameters.pinning_ly2;
		global_variables.pinningsites[k].cosfi = 1.0;
		global_variables.pinningsites[k].sinfi = 0.0;
		global_variables.pinningsites[k].middle_height = global_parameters.pinning_middle_height + gasdev() * global_parameters.pinning_middle_height_spread;
        	if (global_variables.pinningsites[k].middle_height < 0) {global_variables.pinningsites[k].middle_height = 0.0;printf("Warning:middle_height clipped\n");}
		global_variables.pinningsites[k].N_particles_in = 0;
		global_variables.pinningsites[k].particle_in_ID[0] = 0;
		global_variables.pinningsites[k].particle_in_ID[1] = 0;
		global_variables.pinningsites[k].groundstate_dir = 1;
		k++;
		//slanted down pinning site
		global_variables.pinningsites[k].x = ((i*3)+1.25+global_parameters.pinning_xshift) * global_parameters.pinning_lattice_ax;
		global_variables.pinningsites[k].y = ((j*2)+0.5+global_parameters.pinning_yshift) * global_parameters.pinning_lattice_ay;
		global_variables.pinningsites[k].lx = global_parameters.pinning_lx2;
		global_variables.pinningsites[k].ly = global_parameters.pinning_ly2;
		global_variables.pinningsites[k].sinfi = -sqrt(3)/2.0;
		global_variables.pinningsites[k].cosfi = 0.5;
		global_variables.pinningsites[k].middle_height = global_parameters.pinning_middle_height + gasdev() * global_parameters.pinning_middle_height_spread;
        	if (global_variables.pinningsites[k].middle_height < 0) {global_variables.pinningsites[k].middle_height = 0.0;printf("Warning:middle_height clipped\n");}
		global_variables.pinningsites[k].N_particles_in = 0;
		global_variables.pinningsites[k].particle_in_ID[0] = 0;
		global_variables.pinningsites[k].particle_in_ID[1] = 0;
		global_variables.pinningsites[k].groundstate_dir = 1;
		k++;
		//slanted up pinning site
		global_variables.pinningsites[k].x = ((i*3)+1.25+global_parameters.pinning_xshift) * global_parameters.pinning_lattice_ax;
		global_variables.pinningsites[k].y = ((j*2)+1.5+global_parameters.pinning_yshift) * global_parameters.pinning_lattice_ay;
		global_variables.pinningsites[k].lx = global_parameters.pinning_lx2;
		global_variables.pinningsites[k].ly = global_parameters.pinning_ly2;
		global_variables.pinningsites[k].sinfi = sqrt(3)/2.0;
		global_variables.pinningsites[k].cosfi = 0.5;
		global_variables.pinningsites[k].middle_height = global_parameters.pinning_middle_height + gasdev() * global_parameters.pinning_middle_height_spread;
        	if (global_variables.pinningsites[k].middle_height < 0) {global_variables.pinningsites[k].middle_height = 0.0;printf("Warning:middle_height clipped\n");}
		global_variables.pinningsites[k].N_particles_in = 0;
		global_variables.pinningsites[k].particle_in_ID[0] = 0;
		global_variables.pinningsites[k].particle_in_ID[1] = 0;
		global_variables.pinningsites[k].groundstate_dir = 1;
		k++;
		//second horizontal pinning site
		global_variables.pinningsites[k].x = ((i*3)+2+global_parameters.pinning_xshift) * global_parameters.pinning_lattice_ax;
		global_variables.pinningsites[k].y = ((j*2)+global_parameters.pinning_yshift) * global_parameters.pinning_lattice_ay;
		global_variables.pinningsites[k].lx = global_parameters.pinning_lx2;
		global_variables.pinningsites[k].ly = global_parameters.pinning_ly2;
		global_variables.pinningsites[k].cosfi = 1.0;
		global_variables.pinningsites[k].sinfi = 0.0;
		global_variables.pinningsites[k].middle_height = global_parameters.pinning_middle_height + gasdev() * global_parameters.pinning_middle_height_spread;
        	if (global_variables.pinningsites[k].middle_height < 0) {global_variables.pinningsites[k].middle_height = 0.0;printf("Warning:middle_height clipped\n");}
		global_variables.pinningsites[k].N_particles_in = 0;
		global_variables.pinningsites[k].particle_in_ID[0] = 0;
		global_variables.pinningsites[k].particle_in_ID[1] = 0;
		global_variables.pinningsites[k].groundstate_dir = 1;
		k++;
		//second slanted down pinning site
		global_variables.pinningsites[k].x = ((i*3)+2.75+global_parameters.pinning_xshift) * global_parameters.pinning_lattice_ax;
		global_variables.pinningsites[k].y = ((j*2)+1.5+global_parameters.pinning_yshift) * global_parameters.pinning_lattice_ay;
		global_variables.pinningsites[k].lx = global_parameters.pinning_lx2;
		global_variables.pinningsites[k].ly = global_parameters.pinning_ly2;
		global_variables.pinningsites[k].sinfi = -sqrt(3)/2.0;
		global_variables.pinningsites[k].cosfi = 0.5;
		global_variables.pinningsites[k].middle_height = global_parameters.pinning_middle_height + gasdev() * global_parameters.pinning_middle_height_spread;
        	if (global_variables.pinningsites[k].middle_height < 0) {global_variables.pinningsites[k].middle_height = 0.0;printf("Warning:middle_height clipped\n");}
		global_variables.pinningsites[k].N_particles_in = 0;
		global_variables.pinningsites[k].particle_in_ID[0] = 0;
		global_variables.pinningsites[k].particle_in_ID[1] = 0;
		global_variables.pinningsites[k].groundstate_dir = 1;
		k++;
		//slanted up pinning site
		global_variables.pinningsites[k].x = ((i*3)+2.75+global_parameters.pinning_xshift) * global_parameters.pinning_lattice_ax;
		global_variables.pinningsites[k].y = ((j*2)+0.5+global_parameters.pinning_yshift) * global_parameters.pinning_lattice_ay;
		global_variables.pinningsites[k].lx = global_parameters.pinning_lx2;
		global_variables.pinningsites[k].ly = global_parameters.pinning_ly2;
		global_variables.pinningsites[k].sinfi = sqrt(3)/2.0;
		global_variables.pinningsites[k].cosfi = 0.5;
		global_variables.pinningsites[k].middle_height = global_parameters.pinning_middle_height + gasdev() * global_parameters.pinning_middle_height_spread;
        	if (global_variables.pinningsites[k].middle_height < 0) {global_variables.pinningsites[k].middle_height = 0.0;printf("Warning:middle_height clipped\n");}
		global_variables.pinningsites[k].N_particles_in = 0;
		global_variables.pinningsites[k].particle_in_ID[0] = 0;
		global_variables.pinningsites[k].particle_in_ID[1] = 0;
		global_variables.pinningsites[k].groundstate_dir = 1;
		k++;
		}

fflush(global_parameters.diagnostics_file);
}


/* ==============================================================================================
void initialize_particles()
initializaes particles
============================================================================================== */
void initialize_particles()
{

printf("%s %s\n",global_parameters.particle_type,global_parameters.particle_lattice_type);
fprintf(global_parameters.diagnostics_file,"%s %s\n",global_parameters.particle_type,global_parameters.particle_lattice_type);
if ( (strcmp(global_parameters.particle_type,"spin_ice_particle")==0)
		&&(strcmp(global_parameters.particle_lattice_type,"in_spin_ice")==0)
		&&(strcmp(global_parameters.pinning_lattice_type,"square")==0) )
			{

			initialize_square_particle();

                //end speed
                //initialize_single_defect_line(20,20,20,1);

                //four lines figure 2a
                //initialize_single_defect_line(22,11,30,1);
                //initialize_single_defect_line(17,15,30,1);
                //initialize_single_defect_line(8,33,20,2);
                //initialize_single_defect_line(13,37,20,2);

                /*
                this was for the zigzag idea
                initialize_single_defect_line(10,10,6,1);
                initialize_single_defect_line(16,17,6,3);
                initialize_single_defect_line(17,17,6,1);
                initialize_single_defect_line(23,24,6,3);
                initialize_single_defect_line(24,24,6,1);
                initialize_single_defect_line(30,31,6,3);
                initialize_single_defect_line(31,31,6,1);
                initialize_single_defect_line(37,38,6,3);
                initialize_single_defect_line(38,38,6,1);
                initialize_single_defect_line(44,45,6,3);
                */


            //initialize_single_defect_line(10,36,40,2);

			//initialize_defect_line_square();
			}

if ( (strcmp(global_parameters.particle_type,"spin_ice_particle")==0)
		&&(strcmp(global_parameters.particle_lattice_type,"in_spin_ice")==0)
		&&(strcmp(global_parameters.pinning_lattice_type,"hexagonal")==0) ) 	initialize_hexagonal_particle();

}

/* ==============================================================================================
void initialize_which_sites_will_be_double_defects(int empty_defects_exist,int double_defects_exist,int defect_number)
picks the IDs of the spots where the defects will be
============================================================================================== */
void initialize_which_sites_will_be_double_defects(int empty_defects_exist,int double_defects_exist)
{
int i,j;
int N;
int proposed_ID;
int already_ID;

if (global_parameters.particle_defect_number<=global_variables.N_pinningsites)
	global_parameters.defect_pinningsite_IDs = (int *) malloc(global_parameters.particle_defect_number*sizeof(int));
else
	{
	printf("Too many defects are specified %d defects on %d pinningsites\n",global_parameters.particle_defect_number,global_variables.N_pinningsites);
	fprintf(global_parameters.diagnostics_file,"Too many defects are specified %d defects on %d pinningsites\n",global_parameters.particle_defect_number,global_variables.N_pinningsites);
	fflush(global_parameters.diagnostics_file);
	exit(1);
	}

N=0;

//fill the IDs with the pinning site IDs
//that will be the double defects
for(i=0;i<global_parameters.particle_defect_number;i++)
	{
	do	{
		already_ID = 0;
		proposed_ID  = (int)floor(global_variables.N_pinningsites*Rand());

		for(j=0;j<i;j++)
			if (global_parameters.defect_pinningsite_IDs[j]==proposed_ID)
				{
				already_ID=1;
				break;
				}
		} while (already_ID==1);
	global_parameters.defect_pinningsite_IDs[N]=proposed_ID;
	N++;
	}

//printf("Number of defects selected: %d",N);
//for(i=0;i<N;i++)
//	printf("Defect Id %d\n",global_parameters.defect_pinningsite_IDs[i]);

}

/* ==============================================================================================
int check_site_if_double_defect(int pinningsite_ID)
checks if the site is a double defect site 1 if yes 0 is no
============================================================================================== */
int check_site_if_double_defect(int pinningsite_ID)
{
int i;
int isitdoubledefect;

isitdoubledefect = 0;
for(i=0;i<global_parameters.particle_defect_number;i++)
	if (global_parameters.defect_pinningsite_IDs[i]==pinningsite_ID)
		{
		isitdoubledefect = 1;
		break;
		}
return isitdoubledefect;
}


/* ==============================================================================================
void initialize_square_particle()
initializaes particles in a square lattice corresponding to pinning sites
will
============================================================================================== */
void initialize_square_particle()
{
int initial_state;

int double_defects_exist;
int empty_defects_exist;
double defect_probability;

double shift_dir;

int i,j,k;

int count_doubles_placed;

printf("==============================================================================================\n");
fprintf(global_parameters.diagnostics_file,"==============================================================================================\n");
printf("Initializing square ice particles in pinning sites\n");
fprintf(global_parameters.diagnostics_file,"Initializing square ice particles in pinning sites\n");

if (strcmp(global_parameters.particle_starting_position,"random")==0) initial_state = 0;
else if (strcmp(global_parameters.particle_starting_position,"ground_state")==0) initial_state = 1;
else if (strcmp(global_parameters.particle_starting_position,"biased")==0) initial_state = 2;
else
	{
	printf("Unrecognized initial state\n");
	fprintf(global_parameters.diagnostics_file,"Unrecognized initial state\n");
	printf("Allowed keywords: random, ground_state, biased\n");
	fprintf(global_parameters.diagnostics_file,"Allowed keywords: random, ground_state, biased\n");
	exit(1);
	}

printf("Initializing particles in the %s(%d) configuration\n",global_parameters.particle_starting_position,initial_state);
fprintf(global_parameters.diagnostics_file,"Initializing particles in the %s(%d) configuration\n",global_parameters.particle_starting_position,initial_state);

double_defects_exist = 0;
empty_defects_exist = 0;

if ((strcmp(global_parameters.particle_defect_type,"double_defect")==0)&&(global_parameters.particle_defect_number>0)) double_defects_exist=1;
if ((strcmp(global_parameters.particle_defect_type,"empty_defect")==0)&&(global_parameters.particle_defect_number>0)) empty_defects_exist=1;

//this is for decimated lattices
if ( (strcmp(global_parameters.particle_defect_type,"decimated_cw")==0)||
     (strcmp(global_parameters.particle_defect_type,"decimated_ccw")==0)||
     (strcmp(global_parameters.particle_defect_type,"decimated_rand")==0) )
        {
        printf("Special decimated lattice placing\n");
        empty_defects_exist = 2;
        //decimator should have been called already when pinning sites were initialized
        }

//printf("%d %d\n",double_defects_exist,empty_defects_exist);

//check if pinning sites are initialized already
if (global_variables.N_pinningsites==0)
	{
	printf("Pinning sites need to be initialized before particles\n");
	fprintf(global_parameters.diagnostics_file,"Pinning sites need to be initialized before particles\n");
	exit(1);
	}


//double defects or empty defects
//this is the old way of doing the double defects
//the problem with this approach is that is will not produce the correct number of double or empty defects in the system
/*
if ((double_defects_exist==1)||(empty_defects_exist==1))
	{
	defect_probability = global_parameters.particle_defect_number/(double)global_variables.N_pinningsites;
	printf("N defects = %d N sites = %d Probability of defect = %lf\n", global_parameters.particle_defect_number,global_variables.N_pinningsites,defect_probability);
	fprintf(global_parameters.diagnostics_file,	"N defects = %d N sites = %d Probability of defect = %lf\n", global_parameters.particle_defect_number,global_variables.N_pinningsites,defect_probability);
	}
else
	{
	printf("No defects in the spin ice\n");
	fprintf(global_parameters.diagnostics_file, "No defects in the spin ice\n");
	}

global_variables.N_particles = global_variables.N_pinningsites;
if (double_defects_exist==1) global_variables.N_particles += 2*global_parameters.particle_defect_number;
else if (empty_defects_exist==1) global_variables.N_particles -= global_parameters.particle_defect_number/2;
//to make sure there is enough room even if there are more randomly produced defects
printf("Allocating space for possibly max %d particles [%ld bytes]\n",global_variables.N_particles,global_variables.N_particles*sizeof(struct particle_struct));
fprintf(global_parameters.diagnostics_file,"Allocating space for possibly max %d particles [%ld bytes]\n",global_variables.N_particles,global_variables.N_particles*sizeof(struct particle_struct));
*/

global_variables.N_particles = global_variables.N_pinningsites;
if (double_defects_exist==1) global_variables.N_particles += global_parameters.particle_defect_number;
else if (empty_defects_exist==1) global_variables.N_particles -= global_parameters.particle_defect_number;
//to make sure there is enough room even if there are more randomly produced defects
printf("Allocating space for %d particles [%ld bytes]\n",global_variables.N_particles,global_variables.N_particles*sizeof(struct particle_struct));
fprintf(global_parameters.diagnostics_file,"Allocating space for %d particles [%ld bytes]\n",global_variables.N_particles,global_variables.N_particles*sizeof(struct particle_struct));
initialize_which_sites_will_be_double_defects(empty_defects_exist,double_defects_exist);
global_variables.particles = (struct particle_struct *) malloc(global_variables.N_particles*sizeof(struct particle_struct));

count_doubles_placed  = 0;

k=0;
for(i=0;i<global_variables.N_pinningsites;i++)
		{

		if ((empty_defects_exist==1)&&(check_site_if_double_defect(i)==1))
			{
			//empty pinning site
			global_variables.pinningsites[i].N_particles_in = 0;
			//printf("pin %d remains empty\n",i);

			}
        else if ((empty_defects_exist==2)&&(check_site_if_decimated(i)==1))
            {
            global_variables.pinningsites[i].N_particles_in = 0;
            //this is where it gets with the decimated lattices


            }
		else if ((double_defects_exist==1)&&(check_site_if_double_defect(i)==1))
			{
			//doubly occupied pinning site
			global_variables.particles[k].x = global_variables.pinningsites[i].x;
			global_variables.particles[k].y = global_variables.pinningsites[i].y;
			global_variables.particles[k].x += global_variables.pinningsites[i].lx*global_variables.pinningsites[i].cosfi;
			global_variables.particles[k].y += global_variables.pinningsites[i].lx*global_variables.pinningsites[i].sinfi;
			//the pinning site knows the particles that are in it
			global_variables.pinningsites[i].particle_in_ID[global_variables.pinningsites[i].N_particles_in] = k;
			global_variables.pinningsites[i].N_particles_in += 1;
			//the particles know the pinning site they are in
			global_variables.particles[k].pinningsite_ID = i;
			global_variables.particles[k].fx = 0.0;
			global_variables.particles[k].fy = 0.0;
			global_variables.particles[k].q = global_parameters.particle_q;
			global_variables.particles[k].color = 3;
            global_variables.particles[k].E = 0.0;
            global_variables.particles[k].wasclosetosomeone = 0;
			k++;

			global_variables.particles[k].x = global_variables.pinningsites[i].x;
			global_variables.particles[k].y = global_variables.pinningsites[i].y;
			global_variables.particles[k].x -= global_variables.pinningsites[i].lx*global_variables.pinningsites[i].cosfi;
			global_variables.particles[k].y -= global_variables.pinningsites[i].lx*global_variables.pinningsites[i].sinfi;
			//the pinning site knows the particles that are in it
			global_variables.pinningsites[i].particle_in_ID[global_variables.pinningsites[i].N_particles_in] = k;
			global_variables.pinningsites[i].N_particles_in += 1;
			//the particles know the pinning site they are in
			global_variables.particles[k].pinningsite_ID = i;
			global_variables.particles[k].fx = 0.0;
			global_variables.particles[k].fy = 0.0;
			global_variables.particles[k].q = global_parameters.particle_q;
			global_variables.particles[k].color = 3;
            global_variables.particles[k].E = 0.0;
            global_variables.particles[k].wasclosetosomeone = 0;
			k++;

			//count how many double defects were placed
			count_doubles_placed++;

			}
		else
			{
			//single occupancy pinning site
			//initial state can be random, ground state or biased

			if (initial_state==0) //random
				if (Rand()<0.5) shift_dir = -1;
				else			shift_dir = 1;
			else if (initial_state==1) //ground state
				shift_dir = global_variables.pinningsites[i].groundstate_dir;
			else if (initial_state==2) //biased
				shift_dir = 1;

			global_variables.particles[k].x = global_variables.pinningsites[i].x;
			global_variables.particles[k].y = global_variables.pinningsites[i].y;
			global_variables.particles[k].x += shift_dir*global_variables.pinningsites[i].lx*global_variables.pinningsites[i].cosfi;
			global_variables.particles[k].y += shift_dir*global_variables.pinningsites[i].lx*global_variables.pinningsites[i].sinfi;
			//the pinning site knows the particles that are in it
			global_variables.pinningsites[i].particle_in_ID[global_variables.pinningsites[i].N_particles_in] = k;
			global_variables.pinningsites[i].N_particles_in += 1;
			//the particles know the pinning site they are in
			global_variables.particles[k].pinningsite_ID = i;
			global_variables.particles[k].fx = 0.0;
			global_variables.particles[k].fy = 0.0;
			global_variables.particles[k].q = global_parameters.particle_q;
			global_variables.particles[k].color = 2;
			global_variables.particles[k].E = 0.0;
            global_variables.particles[k].wasclosetosomeone = 0;
            k++;
			}
		}

global_variables.N_particles = k;
printf("Put down %d particles on spin ice square lattice\n",global_variables.N_particles);
fprintf(global_parameters.diagnostics_file,"Put down %d particles on spin ice square lattice\n",global_variables.N_particles);
printf("Double defects placed = %d\n",count_doubles_placed);
fprintf(global_parameters.diagnostics_file,"Double defects placed = %d\n",count_doubles_placed);
fflush(global_parameters.diagnostics_file);
}


/* ==============================================================================================
void initialize_defect_line_square()
initializes a defect line in square spin ice
there is no point to initialize it in hexagonal (so in hexagonal it will not work, also in biased)
also it needs some parameters as to where the defect line will go
it will skip doubly occupied sites
============================================================================================== */
void initialize_defect_line_square()
{
int pinningsite_ID;
int k[500], shift_dir;
int i;
int start_pos,curr_pos;
int defect_length;


printf("Initializing defect line\n");
fprintf(global_parameters.diagnostics_file,"Initializing defect line\n");

defect_length = 40;

//start_pos = 2116 + 46 - 46*2*20 - 20*2; //in the middle currently with a 80 long line

//center of sample horizontal pin
start_pos = global_parameters.pinning_lattice_Nx * global_parameters.pinning_lattice_Ny + global_parameters.pinning_lattice_Ny;
start_pos -= global_parameters.pinning_lattice_Nx * defect_length/2; //move back in the x direction half of defect length
start_pos -= defect_length/2; //move back in the y direction half of defect length

curr_pos = start_pos;

//curr_pos  = start_pos - 8;

for(i=0;i<defect_length;i++)
	{
 	if (i%2==0)
		{
		k[i] = curr_pos;
		curr_pos += global_parameters.pinning_lattice_Ny*2+1;
		}
 	else	{
		k[i] = curr_pos;
		curr_pos += 1;
		}
	}

/*
curr_pos = start_pos + 6;

for(i=defect_length;i<2*defect_length;i++)
    {
        if (i%2==0)
        {
            k[i] = curr_pos;
            curr_pos += global_parameters.pinning_lattice_Ny*2+1;
        }
        else	{
            k[i] = curr_pos;
            curr_pos += 1;
        }
    }
*/


for(i=0;i<defect_length;i++)
	{
	//printf("old coordinates: %lf %lf\n",global_variables.particles[k[i]].x,global_variables.particles[k[i]].y);

	pinningsite_ID = global_variables.particles[k[i]].pinningsite_ID;

	//shift_dir = global_variables.pinningsites[pinningsite_ID].groundstate_dir;
	//for the case when we have GS and we put in a biased
	shift_dir = -global_variables.pinningsites[pinningsite_ID].groundstate_dir;

	global_variables.particles[k[i]].x = global_variables.pinningsites[pinningsite_ID].x;
	global_variables.particles[k[i]].y = global_variables.pinningsites[pinningsite_ID].y;
	global_variables.particles[k[i]].x += shift_dir*global_variables.pinningsites[pinningsite_ID].lx*global_variables.pinningsites[pinningsite_ID].cosfi;
	global_variables.particles[k[i]].y += shift_dir*global_variables.pinningsites[pinningsite_ID].lx*global_variables.pinningsites[pinningsite_ID].sinfi;

	//printf("new coordinates: %lf %lf\n",global_variables.particles[k[i]].x,global_variables.particles[k[i]].y);
	}
}

/* ==============================================================================================
initialize_single_defect_line(int defect_line_lenght, int x0, int y0, int direction)
initializes a single defect line
starts at vertex position x0,y0
direction can be 4 possible values
0 going up and right
1 going down and right
2 going down and left
3 going up and left
============================================================================================== */
void initialize_single_defect_line(int x0, int y0, int defect_line_length, int direction)
{
    int start_pos,curr_pos;
    int k[500], shift_dir;
    int i;
    int pinningsite_ID;

    //figure out the starting position of the defect line, which pinningsite it is
    if (!((x0>0)&&(x0<global_parameters.pinning_lattice_Nx)
        &&(y0>0)&&(y0<global_parameters.pinning_lattice_Ny)))
        {
            printf("x0 = %d y0 = %d is not a valid staring position for the defect line\n",x0,y0);
            fprintf(global_parameters.diagnostics_file,
                    "x0 = %d y0 = %d is not a valid staring position for the defect line\n",x0,y0);
        }
    else
    {
        start_pos = x0 * global_parameters.pinning_lattice_Ny * 2 + y0 * 2;

        if ((direction==3)||(direction==4)) start_pos -= 2*global_parameters.pinning_lattice_Ny;

        curr_pos = start_pos;

        //get the positions of the pinning sites that will be modified
        //into the array k[i]
        for(i=0;i<defect_line_length;i++)
        {
            //i am on a horizontal pinning site
            if (i%2==0)
            {
                k[i] = curr_pos;
                if (direction==1) curr_pos += global_parameters.pinning_lattice_Ny*2+1;
                if (direction==2) curr_pos += global_parameters.pinning_lattice_Ny*2-1;
                if (direction==3) curr_pos += -1;
                if (direction==4) curr_pos += 1;
            }
            else
            {
                k[i] = curr_pos;
                if (direction==1) curr_pos += 1;
                if (direction==2) curr_pos += -1;
                if (direction==3) curr_pos += -global_parameters.pinning_lattice_Ny*2-1;
                if (direction==4) curr_pos += -global_parameters.pinning_lattice_Ny*2+1;
            }
        }

        //flip these pinning sites
        for(i=0;i<defect_line_length;i++)
        {
            if ((k[i]>=0)&&(k[i]<global_variables.N_pinningsites))
                {
                    pinningsite_ID = global_variables.particles[k[i]].pinningsite_ID;

                    shift_dir = -global_variables.pinningsites[pinningsite_ID].groundstate_dir;

                    global_variables.particles[k[i]].x = global_variables.pinningsites[pinningsite_ID].x;
                    global_variables.particles[k[i]].y = global_variables.pinningsites[pinningsite_ID].y;
                    global_variables.particles[k[i]].x +=
                    shift_dir*global_variables.pinningsites[pinningsite_ID].lx*global_variables.pinningsites[pinningsite_ID].cosfi;
                    global_variables.particles[k[i]].y +=
                    shift_dir*global_variables.pinningsites[pinningsite_ID].lx*global_variables.pinningsites[pinningsite_ID].sinfi;
                }
            else
                {
                    printf("Defect line cannot continue at flip %d\n",i);
                    fprintf(global_parameters.diagnostics_file,"Defect line cannot continue at flip %d\n",i);
                }
        }



    }



}





/* ==============================================================================================
void initialize_hexagonal_particle()
initializaes particles
============================================================================================== */
void initialize_hexagonal_particle()
{

int initial_state;

int double_defects_exist;
int empty_defects_exist;
double defect_probability;

double shift_dir;

int i,j,k;

int count_doubles_placed;

printf("==============================================================================================\n");
fprintf(global_parameters.diagnostics_file,"==============================================================================================\n");
printf("Initializing hexagonal ice particles in pinning sites\n");
fprintf(global_parameters.diagnostics_file,"Initializing hexagonal ice particles in pinning sites\n");

if (strcmp(global_parameters.particle_starting_position,"random")==0) initial_state = 0;
else if (strcmp(global_parameters.particle_starting_position,"ground_state")==0) initial_state = 1;
else if (strcmp(global_parameters.particle_starting_position,"biased")==0) initial_state = 2;
else
	{
	printf("Unrecognized initial state\n");
	fprintf(global_parameters.diagnostics_file,"Unrecognized initial state\n");
	printf("Allowed keywords: random, ground_state, biased\n");
	fprintf(global_parameters.diagnostics_file,"Allowed keywords: random, ground_state, biased\n");
	exit(1);
	}

printf("Initializing particles in the %s(%d) configuration\n",global_parameters.particle_starting_position,initial_state);
fprintf(global_parameters.diagnostics_file,"Initializing particles in the %s(%d) configuration\n",global_parameters.particle_starting_position,initial_state);

double_defects_exist = 0;
empty_defects_exist = 0;
if ((strcmp(global_parameters.particle_defect_type,"double_defect")==0)&&(global_parameters.particle_defect_number>0)) double_defects_exist=1;
if ((strcmp(global_parameters.particle_defect_type,"empty_defect")==0)&&(global_parameters.particle_defect_number>0)) empty_defects_exist=1;

//printf("%d %d\n",double_defects_exist,empty_defects_exist);

//check if pinning sites are initialized already
if (global_variables.N_pinningsites==0)
	{
	printf("Pinning sites need to be initialized before particles\n");
	fprintf(global_parameters.diagnostics_file,"Pinning sites need to be initialized before particles\n");
	exit(1);
	}

//double defects or empty defects
/*
if ((double_defects_exist==1)||(empty_defects_exist==1))
	{
	defect_probability = global_parameters.particle_defect_number/(double)global_variables.N_pinningsites;
	if (defect_probability>1.0)
			{
			defect_probability = 1.0;
			global_parameters.particle_defect_number = global_variables.N_pinningsites;
			printf("Defect number will be full double defects at %d number of defects\n",global_variables.N_pinningsites);
			fprintf(global_parameters.diagnostics_file,"Defect number will be full double defects at %d number of defects\n",global_variables.N_pinningsites);
			}
	printf("N defects = %d N sites = %d Probability of defect = %lf\n", global_parameters.particle_defect_number,global_variables.N_pinningsites,defect_probability);
	fprintf(global_parameters.diagnostics_file,"N defects = %d N sites = %d Probability of defect = %lf\n", global_parameters.particle_defect_number,global_variables.N_pinningsites,defect_probability);
	}
else
	{
	printf("No defects in the spin ice\n");
	fprintf(global_parameters.diagnostics_file,"No defects in the spin ice\n");
	}

global_variables.N_particles = global_variables.N_pinningsites;
if (double_defects_exist==1) global_variables.N_particles += 2*global_parameters.particle_defect_number;
else if (empty_defects_exist==1) global_variables.N_particles -= global_parameters.particle_defect_number/2;
//to make sure there is enough room even if there are more randomly produced defects
printf("Allocating space for possibly max %d particles [%ld bytes]\n",global_variables.N_particles,global_variables.N_particles*sizeof(struct particle_struct));
fprintf(global_parameters.diagnostics_file,"Allocating space for possibly max %d particles [%ld bytes]\n",global_variables.N_particles,global_variables.N_particles*sizeof(struct particle_struct));
*/

global_variables.N_particles = global_variables.N_pinningsites;
if (double_defects_exist==1) global_variables.N_particles += global_parameters.particle_defect_number;
else if (empty_defects_exist==1) global_variables.N_particles -= global_parameters.particle_defect_number;
//to make sure there is enough room even if there are more randomly produced defects
printf("Allocating space for %d particles [%ld bytes]\n",global_variables.N_particles,global_variables.N_particles*sizeof(struct particle_struct));
fprintf(global_parameters.diagnostics_file,"Allocating space for %d particles [%ld bytes]\n",global_variables.N_particles,global_variables.N_particles*sizeof(struct particle_struct));
initialize_which_sites_will_be_double_defects(empty_defects_exist,double_defects_exist);
global_variables.particles = (struct particle_struct *) malloc(global_variables.N_particles*sizeof(struct particle_struct));

count_doubles_placed = 0;

k=0;
for(i=0;i<global_variables.N_pinningsites;i++)
		{

	if ((empty_defects_exist==1)&&(check_site_if_double_defect(i)==1))
			{
			//empty pinning site
			global_variables.pinningsites[i].N_particles_in = 0;
			//printf("pin %d remains empty\n",i);
			}
		else if ((double_defects_exist==1)&&(check_site_if_double_defect(i)==1))
			{
			//doubly occupied pinning site
			global_variables.particles[k].x = global_variables.pinningsites[i].x;
			global_variables.particles[k].y = global_variables.pinningsites[i].y;

			//printf("%d %d %lf %lf %lf %lf\n",k,i,global_variables.particles[k].x,global_variables.particles[k].y, global_variables.pinningsites[i].x,global_variables.pinningsites[i].y);
			global_variables.particles[k].x += global_variables.pinningsites[i].lx*global_variables.pinningsites[i].cosfi;
			global_variables.particles[k].y += global_variables.pinningsites[i].lx*global_variables.pinningsites[i].sinfi;


			//the pinning site knows the particles that are in it
			global_variables.pinningsites[i].particle_in_ID[global_variables.pinningsites[i].N_particles_in] = k;
			global_variables.pinningsites[i].N_particles_in += 1;
			//the particles know the pinning site they are in
			global_variables.particles[k].pinningsite_ID = i;
			global_variables.particles[k].fx = 0.0;
			global_variables.particles[k].fy = 0.0;
			global_variables.particles[k].q = global_parameters.particle_q;
			global_variables.particles[k].color = 3;
			k++;


			global_variables.particles[k].x = global_variables.pinningsites[i].x;
			global_variables.particles[k].y = global_variables.pinningsites[i].y;
			global_variables.particles[k].x -= global_variables.pinningsites[i].lx*global_variables.pinningsites[i].cosfi;
			global_variables.particles[k].y -= global_variables.pinningsites[i].lx*global_variables.pinningsites[i].sinfi;
			//the pinning site knows the particles that are in it
			global_variables.pinningsites[i].particle_in_ID[global_variables.pinningsites[i].N_particles_in] = k;
			global_variables.pinningsites[i].N_particles_in += 1;
			//the particles know the pinning site they are in
			global_variables.particles[k].pinningsite_ID = i;
			global_variables.particles[k].fx = 0.0;
			global_variables.particles[k].fy = 0.0;
			global_variables.particles[k].q = global_parameters.particle_q;
			global_variables.particles[k].color = 3;
			k++;

			count_doubles_placed++;

			}
		else
			{
			//single occupancy pinning site
			if (initial_state==0) //random
				if (Rand()<0.5) shift_dir = -1;
				else			shift_dir = 1;
			else if (initial_state==1) //ground state
				shift_dir = global_variables.pinningsites[i].groundstate_dir;
			else if (initial_state==2) //biased
				shift_dir = 1;

			global_variables.particles[k].x = global_variables.pinningsites[i].x;
			global_variables.particles[k].y = global_variables.pinningsites[i].y;
			global_variables.particles[k].x += shift_dir*global_variables.pinningsites[i].lx*global_variables.pinningsites[i].cosfi;
			global_variables.particles[k].y += shift_dir*global_variables.pinningsites[i].lx*global_variables.pinningsites[i].sinfi;
			//the pinning site knows the particles that are in it
			global_variables.pinningsites[i].particle_in_ID[global_variables.pinningsites[i].N_particles_in] = k;
			global_variables.pinningsites[i].N_particles_in += 1;
			//the particles know the pinning site they are in
			global_variables.particles[k].pinningsite_ID = i;
			global_variables.particles[k].fx = 0.0;
			global_variables.particles[k].fy = 0.0;
			global_variables.particles[k].q = global_parameters.particle_q;
			global_variables.particles[k].color = 2;
			k++;
			}
		}

global_variables.N_particles = k;
printf("Placed %d particles on spin ice hexagonal lattice\n",global_variables.N_particles);
fprintf(global_parameters.diagnostics_file,"Placed %d particles on spin ice hexagonal lattice\n",global_variables.N_particles);
printf("Double defects placed = %d\n",count_doubles_placed);
fprintf(global_parameters.diagnostics_file,"Double defects placed = %d\n",count_doubles_placed);
//printf("SX SY %lf %lf\n",global_parameters.SX,global_parameters.SY);
fflush(global_parameters.diagnostics_file);
}

/* ==============================================================================================
void initialize_files()
============================================================================================== */
void initialize_files()
{

global_parameters.movie_file = fopen(global_parameters.movie_file_name,"wb");
if (global_parameters.movie_file ==NULL)
	{
	printf("Could not create file %s\n",global_parameters.movie_file_name);
	fprintf(global_parameters.diagnostics_file,"Could not create file %s\n",global_parameters.movie_file_name);
	exit(1);
	}
printf("Created file %s\n",global_parameters.movie_file_name);
fprintf(global_parameters.diagnostics_file,"Created file %s\n",global_parameters.movie_file_name);

global_parameters.statistics_file = fopen(global_parameters.statistics_file_name,"wb");
if (global_parameters.statistics_file ==NULL)
	{
	printf("Could not create file %s\n",global_parameters.statistics_file_name);
	fprintf(global_parameters.diagnostics_file,"Could not create file %s\n",global_parameters.statistics_file_name);
	exit(1);
	}
printf("Created file %s\n",global_parameters.statistics_file_name);
fprintf(global_parameters.diagnostics_file,"Created file %s\n",global_parameters.statistics_file_name);
fflush(global_parameters.diagnostics_file);
}

/* ==============================================================================================
void close_files()
============================================================================================== */
void close_files()
{

fclose(global_parameters.movie_file);
fclose(global_parameters.statistics_file);
fclose(global_parameters.diagnostics_file);
}

/* ==============================================================================================
void initialize_tabulated_forces()
============================================================================================== */
void initialize_tabulated_forces()
{
int i;
double low_limit;
double high_limit;
double step_size;
double dr;
double dr2;
double f;
FILE *tesztf;

low_limit = global_parameters.r_cutoff_short2;
high_limit = global_parameters.r_cutoff_long2;
global_variables.N_tabulated = 10000;
step_size = (high_limit - low_limit)/(double)global_variables.N_tabulated;
global_variables.step_tabulated = step_size;

printf("Tabulating pairwise interaction forces\n");
fprintf(global_parameters.diagnostics_file,"Tabulating pairwise interaction forces\n");
printf("Distances  : %lf to %lf\n",sqrt(low_limit),sqrt(high_limit));
printf("Distances^2: %lf to %lf\n",low_limit,high_limit);
printf("N_tabulated = %d\n",global_variables.N_tabulated);

printf("Tabulating step size in dr2 = %lf\n",step_size);
fprintf(global_parameters.diagnostics_file,"Tabulating step size in dr2 = %lf\n",step_size);

printf("Tabulating force type: %s\n",global_parameters.interaction_type);
fprintf(global_parameters.diagnostics_file,"Tabulating force type: %s\n",global_parameters.interaction_type);

global_variables.f_tabulated = (double *)malloc(global_variables.N_tabulated*sizeof(double));

//tesztf = fopen("tabulalt.txt","wt");
//Coulomb so far
for(i=0;i<10000;i++)
	{
	dr2 = low_limit + i*step_size;
	dr = sqrt(dr2);

	//this is for the Yukawa screened force
	if (strcmp(global_parameters.interaction_type,"Screened_Coulomb")==0)
		f = exp(-global_parameters.particle_a0 * dr) * (1.0 + global_parameters.particle_a0*dr)/dr2;
	//this is for the Magnetic colloids
	else if (strcmp(global_parameters.interaction_type,"Magnetic")==0)
		f = global_parameters.particle_particle_force_term /dr2/dr2 ;

	global_variables.f_tabulated[i] = f/dr;
	//fprintf(tesztf,"%lf %lf %lf\n",dr,global_variables.f_tabulated[i],global_variables.f_tabulated[i]*dr);
	}

printf("Largest force  at dr = %5.2lf is %lf\n",
	sqrt(low_limit),global_variables.f_tabulated[0]*sqrt(low_limit));

printf("Smallest force at dr = %5.2lf is %lf\n",
	sqrt(high_limit),global_variables.f_tabulated[global_variables.N_tabulated-1]*sqrt(high_limit));


//fclose(tesztf);
fflush(global_parameters.diagnostics_file);
}


/* ==============================================================================================
void initialize_vertices()
============================================================================================== */
void initialize_vertices()
{
if ( (strcmp(global_parameters.particle_type,"spin_ice_particle")==0)
		&&(strcmp(global_parameters.particle_lattice_type,"in_spin_ice")==0)
		&&(strcmp(global_parameters.pinning_lattice_type,"square")==0) )	initialize_square_vertices();

if ( (strcmp(global_parameters.particle_type,"spin_ice_particle")==0)
		&&(strcmp(global_parameters.particle_lattice_type,"in_spin_ice")==0)
		&&(strcmp(global_parameters.pinning_lattice_type,"hexagonal")==0) ) 	initialize_hexagonal_vertices();
}



/* ==============================================================================================
void initialize_square_vertices()
============================================================================================== */
void initialize_square_vertices()
{
int i,j,k;
int ii;
double dx,dy;

printf("==============================================================================================\n");
fprintf(global_parameters.diagnostics_file,"==============================================================================================\n");
printf("Initializing square ice vertices\n");
fprintf(global_parameters.diagnostics_file,"Initializing square ice vertices\n");

global_variables.N_vertices = global_parameters.pinning_lattice_Nx * global_parameters.pinning_lattice_Ny;

printf("Number of vertices on a square %d %d lattice is = %d\n",
	global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_Ny,
	global_variables.N_vertices);
fprintf(global_parameters.diagnostics_file,"Number of vertices on a square %d %d lattice is = %d\n",
	global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_Ny,
	global_variables.N_vertices);

printf("Allocating space for %d vertices [%ld bytes]\n",global_variables.N_vertices,global_variables.N_vertices*sizeof(struct vertex_struct));
fprintf(global_parameters.diagnostics_file,"Allocating space for %d vertices [%ld bytes]\n",global_variables.N_vertices,global_variables.N_vertices*sizeof(struct vertex_struct));

global_variables.vertices = (struct vertex_struct *)malloc(global_variables.N_vertices*sizeof(struct vertex_struct));

k=0;
for(i=0;i<global_parameters.pinning_lattice_Nx;i++)
	for(j=0;j<global_parameters.pinning_lattice_Ny;j++)
		{
		//vertex position
		global_variables.vertices[k].x = (i+global_parameters.pinning_xshift) * global_parameters.pinning_lattice_ax;
		global_variables.vertices[k].y = (j+global_parameters.pinning_yshift) * global_parameters.pinning_lattice_ay;
		global_variables.vertices[k].chess_table_color = (i+j)%2;
		global_variables.vertices[k].N_particles_in = 0;

		//find nearby particles
		if (global_parameters.pinning_lattice_ax!=global_parameters.pinning_lattice_ay)
			{
			printf("Vertex initialize particle finding routine can't handle different lattice ax,ay yet\n");
			fprintf(global_parameters.diagnostics_file,"Vertex initialize particle finding routine can't handle different lattice ax,ay yet\n");
			exit(1);
			}

		//find nearby particles
		for(ii=0;ii<global_variables.N_particles;ii++)
			{
			dx = global_variables.particles[ii].x - global_variables.vertices[k].x;
			dy = global_variables.particles[ii].y - global_variables.vertices[k].y;
			//PBC fold back
			if (dx>(global_parameters.SX/2.0)) 	dx = dx - global_parameters.SX;
			if (dx<(-global_parameters.SX/2.0)) dx = dx + global_parameters.SX;
			if (dy>(global_parameters.SY/2.0)) 	dy = dy - global_parameters.SY;
			if (dy<(-global_parameters.SY/2.0)) dy = dy + global_parameters.SY;

			if ( (dx*dx+dy*dy) < (global_parameters.pinning_lattice_ax*global_parameters.pinning_lattice_ax) )
				{
				global_variables.vertices[k].N_particles_in++;
				global_variables.vertices[k].particle_in_ID[global_variables.vertices[k].N_particles_in-1] = ii;
				}
			}
		global_variables.vertices[k].type = 0;
		//recalculate_square_vertex_type(k); //this will be done at running time
		k++;
		}

fflush(global_parameters.diagnostics_file);
}

/* ==============================================================================================
void initialize_hexagonal_vertices()
============================================================================================== */
void initialize_hexagonal_vertices()
{

int i,j,k;
int ii;
double dx,dy;

printf("==============================================================================================\n");
fprintf(global_parameters.diagnostics_file,"==============================================================================================\n");
printf("Initializing hexagonal ice vertices\n");
fprintf(global_parameters.diagnostics_file,"Initializing hexagonal ice vertices\n");

global_variables.N_vertices = global_parameters.pinning_lattice_Nx * global_parameters.pinning_lattice_Ny*4;

printf("Number of vertices on a hexagonal %d %d lattice is = %d\n",
	global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_Ny,
	global_variables.N_vertices);
fprintf(global_parameters.diagnostics_file,"Number of vertices on a hexagonal %d %d lattice is = %d\n",
	global_parameters.pinning_lattice_Nx,global_parameters.pinning_lattice_Ny,
	global_variables.N_vertices);

printf("Allocating space for %d vertices [%ld bytes]\n",global_variables.N_vertices,global_variables.N_vertices*sizeof(struct vertex_struct));
fprintf(global_parameters.diagnostics_file,"Allocating space for %d vertices [%ld bytes]\n",global_variables.N_vertices,global_variables.N_vertices*sizeof(struct vertex_struct));

global_variables.vertices = (struct vertex_struct *)malloc(global_variables.N_vertices*sizeof(struct vertex_struct));

k=0;
for(i=0;i<global_parameters.pinning_lattice_Nx;i++)
	for(j=0;j<global_parameters.pinning_lattice_Ny;j++)
		{
		//vertex 1 position
		global_variables.vertices[k].x = ( (i*3) ) * global_parameters.pinning_lattice_ax;
		global_variables.vertices[k].y = ((j*2)+1) * global_parameters.pinning_lattice_ay;

		global_variables.vertices[k].chess_table_color = 0;
		global_variables.vertices[k].N_particles_in = 0;

		//find nearby particles
		for(ii=0;ii<global_variables.N_particles;ii++)
			{
			dx = global_variables.particles[ii].x - global_variables.vertices[k].x;
			dy = global_variables.particles[ii].y - global_variables.vertices[k].y;
			//PBC fold back
			if (dx>(global_parameters.SX/2.0)) 	dx = dx - global_parameters.SX;
			if (dx<(-global_parameters.SX/2.0)) dx = dx + global_parameters.SX;
			if (dy>(global_parameters.SY/2.0)) 	dy = dy - global_parameters.SY;
			if (dy<(-global_parameters.SY/2.0)) dy = dy + global_parameters.SY;

			if ( (dx*dx+dy*dy) < (global_parameters.pinning_lattice_ax*global_parameters.pinning_lattice_ax) )
				{
				global_variables.vertices[k].N_particles_in++;
				global_variables.vertices[k].particle_in_ID[global_variables.vertices[k].N_particles_in-1] = ii;
				}
			}
		global_variables.vertices[k].type = 0;
		k++;

		//vertex 2 position
		global_variables.vertices[k].x = ((i*3)+1) * global_parameters.pinning_lattice_ax;
		global_variables.vertices[k].y = ((j*2)+1) * global_parameters.pinning_lattice_ay;
		global_variables.vertices[k].chess_table_color = 1;
		global_variables.vertices[k].N_particles_in = 0;

		//find nearby particles
		for(ii=0;ii<global_variables.N_particles;ii++)
			{
			dx = global_variables.particles[ii].x - global_variables.vertices[k].x;
			dy = global_variables.particles[ii].y - global_variables.vertices[k].y;
			//PBC fold back
			if (dx>(global_parameters.SX/2.0)) 	dx = dx - global_parameters.SX;
			if (dx<(-global_parameters.SX/2.0)) dx = dx + global_parameters.SX;
			if (dy>(global_parameters.SY/2.0)) 	dy = dy - global_parameters.SY;
			if (dy<(-global_parameters.SY/2.0)) dy = dy + global_parameters.SY;

			if ( (dx*dx+dy*dy) < (global_parameters.pinning_lattice_ax*global_parameters.pinning_lattice_ax) )
				{
				global_variables.vertices[k].N_particles_in++;
				global_variables.vertices[k].particle_in_ID[global_variables.vertices[k].N_particles_in-1] = ii;
				}
			}
		global_variables.vertices[k].type = 0;
		k++;

		//vertex 3 position
		global_variables.vertices[k].x = ( (i*3)+1.5 ) * global_parameters.pinning_lattice_ax;
		global_variables.vertices[k].y = ((j*2)) * global_parameters.pinning_lattice_ay;

		global_variables.vertices[k].chess_table_color = 0;
		global_variables.vertices[k].N_particles_in = 0;

		//find nearby particles
		for(ii=0;ii<global_variables.N_particles;ii++)
			{
			dx = global_variables.particles[ii].x - global_variables.vertices[k].x;
			dy = global_variables.particles[ii].y - global_variables.vertices[k].y;
			//PBC fold back
			if (dx>(global_parameters.SX/2.0)) 	dx = dx - global_parameters.SX;
			if (dx<(-global_parameters.SX/2.0)) dx = dx + global_parameters.SX;
			if (dy>(global_parameters.SY/2.0)) 	dy = dy - global_parameters.SY;
			if (dy<(-global_parameters.SY/2.0)) dy = dy + global_parameters.SY;

			if ( (dx*dx+dy*dy) < (global_parameters.pinning_lattice_ax*global_parameters.pinning_lattice_ax) )
				{
				global_variables.vertices[k].N_particles_in++;
				global_variables.vertices[k].particle_in_ID[global_variables.vertices[k].N_particles_in-1] = ii;
				}
			}
		global_variables.vertices[k].type = 0;
		k++;


		//vertex 4 position
		global_variables.vertices[k].x = ((i*3)+2.5) * global_parameters.pinning_lattice_ax;
		global_variables.vertices[k].y = ((j*2)) * global_parameters.pinning_lattice_ay;
		global_variables.vertices[k].chess_table_color = 0;
		global_variables.vertices[k].N_particles_in = 0;

		//find nearby particles
		for(ii=0;ii<global_variables.N_particles;ii++)
			{
			dx = global_variables.particles[ii].x - global_variables.vertices[k].x;
			dy = global_variables.particles[ii].y - global_variables.vertices[k].y;
			//PBC fold back
			if (dx>(global_parameters.SX/2.0)) 	dx = dx - global_parameters.SX;
			if (dx<(-global_parameters.SX/2.0)) dx = dx + global_parameters.SX;
			if (dy>(global_parameters.SY/2.0)) 	dy = dy - global_parameters.SY;
			if (dy<(-global_parameters.SY/2.0)) dy = dy + global_parameters.SY;

			if ( (dx*dx+dy*dy) < (global_parameters.pinning_lattice_ax*global_parameters.pinning_lattice_ax) )
				{
				global_variables.vertices[k].N_particles_in++;
				global_variables.vertices[k].particle_in_ID[global_variables.vertices[k].N_particles_in-1] = ii;
				}
			}
		global_variables.vertices[k].type = 0;
		k++;
		}

fflush(global_parameters.diagnostics_file);
}


/* ==============================================================================================
void write_contour_file()
============================================================================================== */
void write_contour_file()
{
int i;
FILE *f;
int N_notdecimated;
char filename[100];

strcpy(filename,"plotting/contour_");
strcat(filename,global_parameters.running_name);
strcat(filename,".txt");

f = fopen(filename,"wt");
if (f==NULL)
	{
	printf("Could not create contour file %s\n",filename);
	exit(1);
	}

N_notdecimated = 0;
for(i=0;i<global_variables.N_pinningsites;i++)
    if (global_variables.pinningsites[i].decimated==0) N_notdecimated++;

printf("Not decimated pinning sites for countour file = %d\n",N_notdecimated);

fprintf(f,"%d\n",N_notdecimated);


for(i=0;i<global_variables.N_pinningsites;i++)
    if (global_variables.pinningsites[i].decimated==0)
	{
	fprintf(f,"%e\n",global_variables.pinningsites[i].x);
	fprintf(f,"%e\n",global_variables.pinningsites[i].y);
	fprintf(f,"%e\n",global_variables.pinningsites[i].middle_height);
    fprintf(f,"%e\n",global_variables.pinningsites[i].lx);
    fprintf(f,"%e\n",global_variables.pinningsites[i].cosfi);
    
    }
fclose(f);
printf("Written contour to file %s\n",filename);

}

/* ==============================================================================================
void write_triple_contour_file()
============================================================================================== */
void write_triple_contour_file()
{
int i;
FILE *f;
double newx,newy;
char filename[100];

strcpy(filename,"contour3_");
strcat(filename,global_parameters.running_name);
strcat(filename,".txt");

f=fopen(filename,"wt");

fprintf(f,"%d\n",global_variables.N_pinningsites*3);


for(i=0;i<global_variables.N_pinningsites;i++)
	{

	//first pin regular location in the center
	fprintf(f,"%e\n",global_variables.pinningsites[i].x);
	fprintf(f,"%e\n",global_variables.pinningsites[i].y);
	fprintf(f,"%e\n",global_parameters.pinning_k * global_variables.pinningsites[i].ly);
	fprintf(f,"%e\n",global_variables.pinningsites[i].middle_height);
	fprintf(f,"%e\n",global_variables.pinningsites[i].ly);
	//the other two are offset at the ends of the pinning site

	newx = global_variables.pinningsites[i].x + global_variables.pinningsites[i].cosfi * 0.9 *(global_variables.pinningsites[i].lx-global_variables.pinningsites[i].ly);
	newy = global_variables.pinningsites[i].y + global_variables.pinningsites[i].sinfi * 0.9 *(global_variables.pinningsites[i].lx-global_variables.pinningsites[i].ly);

	fprintf(f,"%e\n",newx);
	fprintf(f,"%e\n",newy);
	fprintf(f,"%e\n",global_parameters.pinning_k * global_variables.pinningsites[i].ly);
	fprintf(f,"%e\n",global_variables.pinningsites[i].middle_height);
	fprintf(f,"%e\n",global_variables.pinningsites[i].ly);

	newx = global_variables.pinningsites[i].x - global_variables.pinningsites[i].cosfi * 0.9*(global_variables.pinningsites[i].lx-global_variables.pinningsites[i].ly);
	newy = global_variables.pinningsites[i].y - global_variables.pinningsites[i].sinfi * 0.9*(global_variables.pinningsites[i].lx-global_variables.pinningsites[i].ly);

	fprintf(f,"%e\n",newx);
	fprintf(f,"%e\n",newy);
	fprintf(f,"%e\n",global_parameters.pinning_k * global_variables.pinningsites[i].ly);
	fprintf(f,"%e\n",global_variables.pinningsites[i].middle_height);
	fprintf(f,"%e\n",global_variables.pinningsites[i].ly);
	}
fclose(f);
}



/* ==============================================================================================
void test_by_coloring()
generic tester by recoloring particles
============================================================================================== */
void test_by_coloring()
{
int i;

for(i=0;i<global_variables.vertices[30].N_particles_in;i++)
	global_variables.particles[global_variables.vertices[30].particle_in_ID[i]].color = 3;

}

/* ==============================================================================================
void write_testfile()
generic tester
============================================================================================== */
void write_testfile()
{
int i;
FILE *f;

f=fopen("testpins.txt","wt");

for(i=0;i<global_variables.N_pinningsites;i++)
	{
	fprintf(f,"%lf %lf\n",global_variables.pinningsites[i].x,global_variables.pinningsites[i].y);
	fprintf(f,"%lf %lf\n",
		global_variables.pinningsites[i].x+global_variables.pinningsites[i].lx*global_variables.pinningsites[i].cosfi,
		global_variables.pinningsites[i].y+global_variables.pinningsites[i].lx*global_variables.pinningsites[i].sinfi);
	fprintf(f,"%lf %lf\n",
		global_variables.pinningsites[i].x-global_variables.pinningsites[i].lx*global_variables.pinningsites[i].cosfi,
		global_variables.pinningsites[i].y-global_variables.pinningsites[i].lx*global_variables.pinningsites[i].sinfi);

	}
fclose(f);

f=fopen("testparticles.txt","wt");

for(i=0;i<global_variables.N_particles;i++)
	{
	fprintf(f,"%f %f\n",global_variables.particles[i].x,global_variables.particles[i].y);
	}
fclose(f);


f=fopen("testvertices.txt","wt");

for(i=0;i<global_variables.N_vertices;i++)
	{
	fprintf(f,"%f %f\n",global_variables.vertices[i].x,global_variables.vertices[i].y);
	}
fclose(f);
}

/* ==============================================================================================
void decimator()
checks if a site should be empty or not based on the site  decimation algorithm
site_number = pinning site number
============================================================================================== */
void decimator()
{
int i;
int neighbors[6];
int nr_trials,nr_max_trials;
int nr_defects,nr_defects_placed;
int picked_site,good_site;

//random decimation up to the defects_number
//it will also stop when there were too many failed attempts
if (strcmp(global_parameters.particle_defect_type,"decimated_rand")==0)
    {
        nr_max_trials = global_parameters.pinning_lattice_Nx * global_parameters.pinning_lattice_Ny * 100;
        nr_defects = global_parameters.particle_defect_number;
        printf("Will attempt to decimate %d times or %d failures\n",nr_defects,nr_max_trials);

        nr_trials = 0;
        nr_defects_placed = 0;

				if (nr_defects>0)
        while ((nr_trials<nr_max_trials)&&(nr_defects_placed<nr_defects))
            {
            //printf("placed=%d tri=%d\n",nr_defects_placed,nr_trials);fflush(stdout);

            //pick a pinning site
            picked_site = (int) floor(Rand() * global_variables.N_pinningsites);


            //debugging purposes
            /*
            if (nr_defects_placed>=210)
             {
             printf("%d : ",picked_site);
             fflush(stdout);
             }
            */

            return_pinningsite_neighbors(picked_site,neighbors);

            /*
            if (nr_defects_placed>=210)
                {
                for(i=0;i<6;i++)
                    printf("%d ",neighbors[i]);

                printf("\n");
                fflush(stdout);
                }
            */

            //check if its neighbors are already decimated or not
            good_site = 1;

            for(i=0;i<6;i++)
                if (global_variables.pinningsites[neighbors[i]].decimated==1) good_site = 0;

            //check if this site itself was picked already
            if (global_variables.pinningsites[picked_site].decimated==1) good_site = 0;


            if (good_site==1)
                {
                global_variables.pinningsites[picked_site].decimated = 1;
                nr_defects_placed++;
                nr_trials = 0;
                }
            else
                {
                nr_trials++;
                }
            }
        printf("Placed a total of %d defects out of %d\n",nr_defects_placed,nr_defects);

        //debug is this really how many defects it placed
        int nr_really_placed = 0;
        for(i=0;i<global_variables.N_pinningsites;i++)
            nr_really_placed += global_variables.pinningsites[i].decimated;
        printf("Really placed %d\n",nr_really_placed);
				printf("macs\n");

        //global_variables.pinningsites[picked_site].decimated = 1;
				fflush(stdout);

        printf("Trials %d out of max %d\n",nr_trials,nr_max_trials);
				fflush(stdout);
    }
/*
i = 2450;
return_pinningsite_neighbors(i,neighbors);
printf("Neighbor testing\n");
for(i=0;i<6;i++)
     {
        printf("neighbor %d = %d\n",i,neighbors[i]);
        global_variables.pinningsites[neighbors[i]].decimated = 1;
     }
*/

/*
if ( (strcmp(global_parameters.particle_defect_type,"decimated_cw")==0)||
     (strcmp(global_parameters.particle_defect_type,"decimated_ccw")==0)||
     (strcmp(global_parameters.particle_defect_type,"decimated_rand")==0) )

for(i=0;i<global_parameters.pinning_lattice_Nx;i++)
	for(j=0;j<global_parameters.pinning_lattice_Ny;j++)

global_variables.pinningsites[k].decimated = 0;
*/


}

/* ==============================================================================================
void return_pinningsite_neighbors(int pinningsite, int neighbors[6])
returns the 6 neighbors of the pinning site
assumes all are in the initialized order before any Morton or Z ordering
============================================================================================== */
void return_pinningsite_neighbors(int pinningsite, int *neighbors)
{

//this is a vertical pinning site
if (pinningsite%2==1)
    {
    //up right
    neighbors[0] = pinningsite + 1;
    //up
    neighbors[1] = pinningsite + 2;
    //up left
    neighbors[2] = pinningsite + 1 - 2 * global_parameters.pinning_lattice_Ny;
    //down left
    neighbors[3] = pinningsite - 1 - 2 * global_parameters.pinning_lattice_Ny;
    //down
    neighbors[4] = pinningsite - 2;
    //down right
    neighbors[5] = pinningsite - 1;

    //PBC checks

    //the pinning site is at the top of the simulated region
    if ((pinningsite + 1) % (2 * global_parameters.pinning_lattice_Ny) == 0)
        {
        neighbors[0] -= 2 * global_parameters.pinning_lattice_Ny;
        neighbors[1] -= 2 * global_parameters.pinning_lattice_Ny;
        neighbors[2] -= 2 * global_parameters.pinning_lattice_Ny;
        }

    //the pinning site is at the bottom of the simulated region
    if ((pinningsite - 1) % (2 * global_parameters.pinning_lattice_Ny) == 0)
        {
        neighbors[4] += 2 * global_parameters.pinning_lattice_Ny;
        }

    //the pinning site is on the right side //no adjustments needed

    //the pinning site is on the left side of the simulated region
    if (pinningsite < (2 * global_parameters.pinning_lattice_Ny))
        {
        neighbors[2] += 2 * global_parameters.pinning_lattice_Ny * global_parameters.pinning_lattice_Nx;
        neighbors[3] += 2 * global_parameters.pinning_lattice_Ny * global_parameters.pinning_lattice_Nx;
        }
    }

//horizontal pinning site
if (pinningsite%2==0)
    {
    //right
    neighbors[0] = pinningsite + 2 * global_parameters.pinning_lattice_Ny;
    //up right
    neighbors[1] = pinningsite + 2 * global_parameters.pinning_lattice_Ny + 1;
    //up left
    neighbors[2] = pinningsite + 1;
    //left
    neighbors[3] = pinningsite - 2 * global_parameters.pinning_lattice_Ny;
    //down left
    neighbors[4] = pinningsite - 1;
    //down right
    neighbors[5] = pinningsite + 2 * global_parameters.pinning_lattice_Ny - 1;

    //PBC checks

    //the pinning site is at the top of the simulated region//no adjustments needed

    //the pinning site is at the bottom of the simulated region
    if ((pinningsite) % (2 * global_parameters.pinning_lattice_Ny) == 0)
        {
        neighbors[4] += 2 * global_parameters.pinning_lattice_Ny;
        neighbors[5] += 2 * global_parameters.pinning_lattice_Ny;
        }

    //the pinning site is on the right side
    if ( 2 * global_parameters.pinning_lattice_Ny * global_parameters.pinning_lattice_Nx - pinningsite
            <= 2 * global_parameters.pinning_lattice_Ny)
        {
        neighbors[0] -= 2 * global_parameters.pinning_lattice_Ny * global_parameters.pinning_lattice_Nx;
        neighbors[1] -= 2 * global_parameters.pinning_lattice_Ny * global_parameters.pinning_lattice_Nx;
        neighbors[5] -= 2 * global_parameters.pinning_lattice_Ny * global_parameters.pinning_lattice_Nx;
        }

    //the pinning site is on the left side of the simulated region
    if (pinningsite < (2 * global_parameters.pinning_lattice_Ny))
        {
        neighbors[3] += 2 * global_parameters.pinning_lattice_Ny * global_parameters.pinning_lattice_Nx;
        }

    }

}

/* ==============================================================================================
return_pinningsite_vertices(int pinningsite, int *vert1, int *vert2)
returns the 2 neighboring vertices of the pinning site
assumes all are in the initialized order before any Morton or Z ordering
assumes return_pinningsite_neighbors works well
============================================================================================== */
void return_pinningsite_vertices(int pinningsite, int *vert1, int *vert2)
{
int neighbors[6];

return_pinningsite_neighbors(pinningsite,neighbors);

//horizontal pinning site
if (pinningsite%2==0)
    {
    *vert1 = pinningsite/2;
    *vert2 = neighbors[0]/2;
    }
//vertical pinning site
else if (pinningsite%1==0)
    {
    *vert1 = pinningsite/2;
    *vert2 = neighbors[1]/2;
    }

}


/* ==============================================================================================
int check_site_if_decimated(int site_number)
checks if a site should be empty or not based on the site  decimation algorithm
site_number = pinning site number
============================================================================================== */
int check_site_if_decimated(int site_number)
{
int decim;
decim = global_variables.pinningsites[site_number].decimated;

return decim;
}


/* ==============================================================================================
void initialize_vertex_z()
initializes the vertex z number
initializes N_neighbor_vertex which stores how many vertices are connected to a given vertex
initializes neighbor_vertex_ID which stores the ID of the neighboring vertices
============================================================================================== */
void initialize_vertex_z()
{
int i,j;
int vert1,vert2;
int vert_ID;

for(i=0;i<global_variables.N_vertices;i++)
    {
    global_variables.vertices[i].z = 4;
    global_variables.vertices[i].N_neighbor_vertices = 0;
    }
    
for(i=0;i<global_variables.N_pinningsites;i++)
        {
            //because this pin remained empty the vertices belonging to this pin are now z=3
            //and not z=4 vertices, so we need to set these vertices as z=3 vertices
            if ((check_site_if_decimated(i)==1))
                {
                return_pinningsite_vertices(i,&vert1,&vert2);
                global_variables.vertices[vert1].z = 3;
                global_variables.vertices[vert2].z = 3;
                }
            else
                {
                //only the non-decimated pinning sites make a connection between vertices
                //and make them neighbors
                return_pinningsite_vertices(i,&vert1,&vert2);
                
                global_variables.vertices[vert1].neighbor_vertex_ID[global_variables.vertices[vert1].N_neighbor_vertices] = vert2;
                global_variables.vertices[vert1].N_neighbor_vertices++;
                global_variables.vertices[vert2].neighbor_vertex_ID[global_variables.vertices[vert2].N_neighbor_vertices] = vert1;
                global_variables.vertices[vert2].N_neighbor_vertices++;
                }
         }
//check what I've done
int totalz3 = 0;
int totalz4 = 0;
double eta, xi;

for(i=0;i<global_variables.N_vertices;i++)
    {
    if (global_variables.vertices[i].z==3) totalz3++;
    else if (global_variables.vertices[i].z==4) totalz4++;
    }

eta = (double)totalz3 / (double)totalz4;
xi = eta/(eta+1.0)/4.0;

global_variables.N_vertices3 = totalz3;
global_variables.N_vertices4 = totalz4;
global_variables.eta_vertices = eta;
global_variables.xi_vertices = xi;

printf("Decimated lattice initialized\nz3 = %d z4 = %d eta = %lf xi=%lf\n",totalz3,totalz4,global_variables.eta_vertices,global_variables.xi_vertices);
printf("Total number of vertices is %d with z4 = %d and z3 = %d\n",global_variables.N_vertices,global_variables.N_vertices4,global_variables.N_vertices3);

//check the vertex neighbors if they are ok
for(i=0;i<global_variables.N_vertices;i++)
    {
         if (global_variables.vertices[i].z != global_variables.vertices[i].N_neighbor_vertices)
            printf("Something is wrong at vertex %d z=%d neighbors=%d\n",
                i,global_variables.vertices[i].z,global_variables.vertices[i].N_neighbor_vertices);
    }

//count how many z3 and z4 neighbors each vertex has
for(i=0;i<global_variables.N_vertices;i++)
    {
    global_variables.vertices[i].N_neighbor_z3 = 0;
    global_variables.vertices[i].N_neighbor_z4 = 0;
    for(j=0;j<global_variables.vertices[i].N_neighbor_vertices;j++)
        {
        vert_ID = global_variables.vertices[i].neighbor_vertex_ID[j];
        if (global_variables.vertices[vert_ID].z==3) global_variables.vertices[i].N_neighbor_z3++;
        else if (global_variables.vertices[vert_ID].z==4) global_variables.vertices[i].N_neighbor_z4++;
        }
    
    }



//test case
//tested and it works
/*
int k=240;

printf("%d %d\n",global_variables.vertices[k].N_neighbor_z3,global_variables.vertices[k].N_neighbor_z4);
for(i=0;i<global_variables.vertices[k].N_neighbor_vertices;i++)
    printf("neighbor vertex is %d [%d]\n",global_variables.vertices[k].neighbor_vertex_ID[i],
        global_variables.vertices[global_variables.vertices[k].neighbor_vertex_ID[i]].z);

*/
}



/* ==============================================================================================
void checkvertexfinder()
debug purposes check if the pinning site finds correct vertices
============================================================================================== */
void checkvertexfinder()
{
int i;

printf("CHECKVERTEXFINDER\n");

for(i=0;i<global_variables.N_pinningsites;i++)
        {
            //because this pin remained empty the vertices belonging to this pin are now z=3
            //and not z=4 vertices, so we need to set these vertices as z=3 vertices

            if ((check_site_if_decimated(i)==1))
            {
            int vert1,vert2;
            return_pinningsite_vertices(i,&vert1,&vert2);
            printf(">>>>>>>>pinnigsite %d vertices: %d %d\n",i,vert1,vert2);
            printf("%lf %lf\n",global_variables.pinningsites[i].x,global_variables.pinningsites[i].y);
            printf("%lf %lf\n",global_variables.vertices[vert1].x,global_variables.vertices[vert1].y);
            printf("%lf %lf\n",global_variables.vertices[vert2].x,global_variables.vertices[vert2].y);
            }
         }
}
