

#include <xmmintrin.h>
#include "global_variables.h"
#include "running.h"
#include "random_numrec.h"
#include <math.h>
#include <string.h>
#include "timing.h"

/* ==============================================================================================
run_system()
============================================================================================== */
void run_system()
{

printf("==============================================================================================\n");
fprintf(global_parameters.diagnostics_file,"==============================================================================================\n");
printf("Running %s\n",global_parameters.running_name);
fprintf(global_parameters.diagnostics_file,"Running %s\n",global_parameters.running_name);

global_variables.max_movement_x = 0.0;
global_variables.max_movement_y = 0.0;

//start time
global_parameters.time = 0;
while(global_parameters.time<=global_parameters.total_runtime)
	{
	//estimate timing
	if (global_parameters.time==global_parameters.echo_write_time)
		estimate_running_time();

	//echo time
	if (global_parameters.time%global_parameters.echo_write_time==0)
		{
		printf("%s %d/%d ",global_parameters.running_name,global_parameters.time,global_parameters.total_runtime);
		//printf("T = %.2lf ",global_parameters.temperature);
		printf("B = %.2lf ",global_parameters.magnetic_field);
        printf("Fx,y = (%.2lf %.2lf) ",global_parameters.external_fx,global_parameters.external_fy);
		printf("Frame = %d/%d\n",global_parameters.time/global_parameters.movie_write_time,global_parameters.total_runtime/global_parameters.movie_write_time);

		fprintf(global_parameters.diagnostics_file,"%s %d/%d ",
                global_parameters.running_name,global_parameters.time,global_parameters.total_runtime);
        //fprintf(global_parameters.diagnostics_file,"T = %.2lf ",global_parameters.temperature);
		fprintf(global_parameters.diagnostics_file,"B = %.2lf ",global_parameters.magnetic_field);
        fprintf(global_parameters.diagnostics_file,"Fx,y = (%.2lf %.2lf) ",
                global_parameters.external_fx,global_parameters.external_fy);
		fprintf(global_parameters.diagnostics_file,"Frame = %d/%d\n",global_parameters.time/global_parameters.movie_write_time,global_parameters.total_runtime/global_parameters.movie_write_time);
		fflush(global_parameters.diagnostics_file);
		}

	//movie
	if (global_parameters.time%global_parameters.movie_write_time==0)
		{
		write_cmovie_frame();
		}

	//testing
	//if (global_parameters.time==2000)
	//	write_xmgrace_files();

	//statistics
	//if (global_parameters.time==0)
	//	{
	//	writeout_hophistogram();
	//	}

	if ((global_parameters.time%global_parameters.statistics_write_time==0)&&(global_parameters.time!=0))
		{
		//writeout_hophistogram();
		write_statistics();
		calculate_temperature();
		//calculate_external();
		}

	calculate_temperature_forces();
	calculate_external_forces();
    calculate_pinning_forces();
    calculate_magnetic_field();
	calculate_pairwise_forces();

	calculate_statistics_per_step();

	//don't move on last frame
	if (global_parameters.time<global_parameters.total_runtime)
		move_particles();

	//next timestep
	global_parameters.time++;
	}

printf("maximum forces fx = %lf fy = %lf\n",global_variables.max_movement_x,global_variables.max_movement_y);
fprintf(global_parameters.diagnostics_file,"maximum forces fx = %lf fy = %lf\n",global_variables.max_movement_x,global_variables.max_movement_y);
fflush(global_parameters.diagnostics_file);
}

/* ==============================================================================================
write_cmovie_frame()
write legacy cmovie format for delpot
============================================================================================== */
void write_cmovie_frame()
{
int i;
float floatholder;
int intholder;

intholder = global_variables.N_particles + global_variables.N_vertices;
fwrite(&intholder,sizeof(int),1,global_parameters.movie_file);

intholder = global_parameters.time;
fwrite(&intholder,sizeof(int),1,global_parameters.movie_file);

for (i=0;i<global_variables.N_particles;i++)
	{
	//color decode
	/*
	switch (global_variables.particles[i].color)
			{
			case 1:
				intholder = 2; // color
				break;
			case 2:
				intholder = 11;
				break;
			default:
				intholder = 3;
			}
	*/
	intholder = global_variables.particles[i].color;
	fwrite(&intholder,sizeof(int),1,global_parameters.movie_file);
	intholder = i;//ID
	fwrite(&intholder,sizeof(int),1,global_parameters.movie_file);
	floatholder = (float)global_variables.particles[i].x;
	fwrite(&floatholder,sizeof(float),1, global_parameters.movie_file);
	floatholder = (float)global_variables.particles[i].y;
	fwrite(&floatholder,sizeof(float),1, global_parameters.movie_file);
	floatholder = 1.0;//cum_disp, cmovie format
	fwrite(&floatholder,sizeof(float),1,global_parameters.movie_file);
	}

for (i=0;i<global_variables.N_vertices;i++)
	{

	//intholder = global_variables.vertices[i].type + 4;//color

	//check if the vertex is z3 or z4 and then color it by the charge
	//coloring by the charge, differently for z3 and z4
	intholder = 1;
    
    
    //coloring scheme
    /*
    4   4 out z4      charge -4
    5   3 out z3      charge -3
    6   3 out 1 in z4 charge -2
    7   2 out 1 in z3 charge -1
    8   2 in 2 out z4 charge 0
    9   2 in 1 out z3 charge 1
    10  3 in 1 out z4 charge 2
    11  3 in z3       charge 3
    12  4 in z4       charge 4
    
    I still want to differentiate GS1,GS2 and biased
    13 GS1
    14 GS2
    15 biased
    so 8 will not be present in this case
    */
    
	if (global_variables.vertices[i].z==3)
		{
        //to never color a z3 vertex light or dark grey because that is both not GS but a 2in 1out state
		if ((global_variables.vertices[i].type ==2)||(global_variables.vertices[i].type ==5)||(global_variables.vertices[i].type ==6)) intholder = 9;
		if (global_variables.vertices[i].type ==1) intholder = 7;
		if (global_variables.vertices[i].type ==3) intholder = 11;
		if (global_variables.vertices[i].type ==0) intholder = 5;
		}
	else if (global_variables.vertices[i].z==4)
		{
        //old scheme no longer used
		//if ((global_variables.vertices[i].type ==2)||(global_variables.vertices[i].type ==5)||(global_variables.vertices[i].type ==6)) intholder = 8;
        if (global_variables.vertices[i].type ==2) intholder = 15;
        if (global_variables.vertices[i].type ==5) intholder = 13;
        if (global_variables.vertices[i].type ==6) intholder = 14;
    
		if (global_variables.vertices[i].type ==1) intholder = 6;
		if (global_variables.vertices[i].type ==3) intholder = 10;
		if (global_variables.vertices[i].type ==0) intholder = 4;
		if (global_variables.vertices[i].type ==4) intholder = 12;
		}


	fwrite(&intholder,sizeof(int),1,global_parameters.movie_file);
	intholder = i;//ID
	fwrite(&intholder,sizeof(int),1,global_parameters.movie_file);
	floatholder = (float)global_variables.vertices[i].x;
	fwrite(&floatholder,sizeof(float),1, global_parameters.movie_file);
	floatholder = (float)global_variables.vertices[i].y;
	fwrite(&floatholder,sizeof(float),1, global_parameters.movie_file);
	floatholder = 1.0;//cum_disp, cmovie format
	fwrite(&floatholder,sizeof(float),1,global_parameters.movie_file);
	}



}


/* ==============================================================================================
write_xmgrace_files()
this is for testing purposes, writes a separate file for the particles, for the pinning sites etc.
============================================================================================== */
void write_xmgrace_files()
{
int i;
FILE *xmgracefile;

xmgracefile = fopen("results/test_particles.dat","wt");
for (i=0;i<global_variables.N_particles;i++)
	fprintf(xmgracefile,"%lf %lf\n",global_variables.particles[i].x,global_variables.particles[i].y);
fclose(xmgracefile);

xmgracefile = fopen("results/test_vertices.dat","wt");
for (i=0;i<global_variables.N_vertices;i++)
	fprintf(xmgracefile,"%lf %lf\n",global_variables.vertices[i].x,global_variables.vertices[i].y);
fclose(xmgracefile);

xmgracefile = fopen("results/test_pinningsites.dat","wt");
for (i=0;i<global_variables.N_pinningsites;i++)
	fprintf(xmgracefile,"%lf %lf\n",global_variables.pinningsites[i].x,global_variables.pinningsites[i].y);
fclose(xmgracefile);
}



/* ==============================================================================================
move_particles()
============================================================================================== */
void move_particles()
{
int i;
for(i=0;i<global_variables.N_particles;i++)
	{
	//move particles

	//printf(">>%lf %lf\n",global_variables.particles[i].fx*global_parameters.dt , global_variables.particles[i].fy*global_parameters.dt);

	//if (fabs(global_variables.particles[i].fx * global_parameters.dt)>max_movement_x) max_movement_x = fabs(global_variables.particles[i].fx * global_parameters.dt);
	//if (fabs(global_variables.particles[i].fy * global_parameters.dt)>max_movement_y) max_movement_y = fabs(global_variables.particles[i].fy * global_parameters.dt);

	global_variables.particles[i].x += global_variables.particles[i].fx * global_parameters.mobility * global_parameters.dt;
	global_variables.particles[i].y += global_variables.particles[i].fy * global_parameters.mobility * global_parameters.dt;

	//PBC check will come here if implemented, for spin ice I don't need it
	if (global_variables.particles[i].x>global_parameters.SX) 	global_variables.particles[i].x -= global_parameters.SX;
	if (global_variables.particles[i].x<0.0) 					global_variables.particles[i].x += global_parameters.SX;
	if (global_variables.particles[i].y>global_parameters.SY) 	global_variables.particles[i].y -= global_parameters.SY;
	if (global_variables.particles[i].y<0.0) 					global_variables.particles[i].y += global_parameters.SY;

	//zero forces
	global_variables.particles[i].fx = 0.0;
	global_variables.particles[i].fy = 0.0;

    //zero energies
    global_variables.particles[i].E = 0.0;
	}

}

/* ==============================================================================================
calculate_pinning_forces()
//rewritten to real units-Barcelona
============================================================================================== */
void calculate_pinning_forces()
{
int i;
int pinning_i;
double dx,dy;
double forcex,forcey;
double lx2,ly2;
double sinfi,cosfi;

double k_spring, k_middle_spring;
double pinning_middle_height;

double dxprime,dyprime;
double forcexprime,forceyprime;

k_spring = global_parameters.pinning_k;
k_middle_spring = global_parameters.pinning_middle_k;

for(i=0;i<global_variables.N_particles;i++)
	{
	pinning_i = global_variables.particles[i].pinningsite_ID;

	dx = global_variables.particles[i].x - global_variables.pinningsites[pinning_i].x;
	dy = global_variables.particles[i].y - global_variables.pinningsites[pinning_i].y;

	//do a PBC check, just in case
	//none of the elongated pinning sites straddles a boundary so this is not necessary
	//they might straddle it when horizontal ... so I need this

	//I actually need to do this now that the particle never leaves the
	//pinning site and is always attracted back
	if (dx>(global_parameters.SX/2.0)) 	dx = dx - global_parameters.SX;
	if (dx<(-global_parameters.SX/2.0)) dx = dx + global_parameters.SX;
	if (dy>(global_parameters.SY/2.0)) 	dy = dy - global_parameters.SY;
	if (dy<(-global_parameters.SY/2.0)) dy = dy + global_parameters.SY;

	forcex = 0.0;
	forcey = 0.0;

	lx2 = global_variables.pinningsites[pinning_i].lx;
	ly2 = global_variables.pinningsites[pinning_i].ly;
	sinfi = global_variables.pinningsites[pinning_i].sinfi;
	cosfi = global_variables.pinningsites[pinning_i].cosfi;
	pinning_middle_height = global_variables.pinningsites[pinning_i].middle_height;


	//no need to check if the particle is in the pinning site
	//it should be there and always stay there
	//(remains permamently no matter how far it goes)

	//rotate all coordinates so the pin lies horizontally

	dxprime =  cosfi*dx + sinfi*dy;
	dyprime = -sinfi*dx + cosfi*dy;

	//spin related calculations
	if (global_variables.pinningsites[pinning_i].N_particles_in==1) calculate_spin_flipping(pinning_i,i,dxprime);

	//the particle is in the elongated middle part
	if (fabs(dxprime)<lx2)
		{
		//middle bump force calculation x direction
		//it pushes from the middle so + sign
		forcexprime = k_middle_spring * pinning_middle_height * dxprime;
		//y direction
		forceyprime = -k_spring*dyprime;
		}
	//the particle is in one one of the round ends
	else
		{
		//modifying dxprime for the circular end
		if (dxprime>0) 	dxprime -= lx2;
		else 		dxprime += lx2;
		//is it inside the circular end? we are not checking this any more
		//the pinning site goes out to infinity in all directions
		//if ((dxprime*dxprime+dyprime*dyprime)<ly2*ly2)
			{
			forcexprime = -k_spring*dxprime;
			forceyprime = -k_spring*dyprime;
			}
		}

	//rotate back
	//testing purposes printf
	//printf("typical pinning force = %15.13lf %15.13lf\n",forcexprime,forceyprime);

	forcex = cosfi*forcexprime - sinfi*forceyprime;
	forcey = sinfi*forcexprime + cosfi*forceyprime;


	if (fabs(forcex )>global_variables.max_movement_x) global_variables.max_movement_x = fabs(forcex );
	if (fabs(forcey )>global_variables.max_movement_y) global_variables.max_movement_y = fabs(forcey );

	global_variables.particles[i].fx += forcex;
	global_variables.particles[i].fy += forcey;
	}

}


/* ==============================================================================================
calculate_temperature_forces()
============================================================================================== */
void calculate_temperature_forces()
{
int i;
double forcex,forcey;

for(i=0;i<global_variables.N_particles;i++)
	{

	forcey = gasdev()*global_parameters.thermal_force_term;
	forcex = gasdev()*global_parameters.thermal_force_term;

	//printf("%lf %lf\n",forcex,forcey);
	//printf("typical temperature force = %15.13lf %15.13lf\n",forcex,forcey);

	global_variables.particles[i].fx += forcex;
	global_variables.particles[i].fy += forcey;
	}

}

/* ==============================================================================================
calculate_external_forces()
============================================================================================== */
void calculate_external_forces()
{
int i;
double forcex,forcey;

for(i=0;i<global_variables.N_particles;i++)
	{

	forcey = global_parameters.external_fx;
	forcex = global_parameters.external_fy;

	global_variables.particles[i].fx += forcex;
	global_variables.particles[i].fy += forcey;
	}

}



/* ==============================================================================================
build_Verlet_list()
============================================================================================== */
void build_Verlet_list()
{
int i,j;
double dx,dy;

global_variables.N_pairs = 0;
global_variables.Verlet_list_i = NULL;
global_variables.Verlet_list_j = NULL;

for(i=0;i<global_variables.N_particles;i++)
	for(j=i+1;j<global_variables.N_particles;j++)
		{
		dx = global_variables.particles[j].x - global_variables.particles[i].x;
		dy = global_variables.particles[j].y - global_variables.particles[i].y;
		//PBC fold back of distances
		if (dx>(global_parameters.SX/2.0)) 	dx = dx - global_parameters.SX;
		if (dx<(-global_parameters.SX/2.0)) dx = dx + global_parameters.SX;
		if (dy>(global_parameters.SY/2.0)) 	dy = dy - global_parameters.SY;
		if (dy<(-global_parameters.SY/2.0)) dy = dy + global_parameters.SY;
		//is it inside the cutoff radius
		if (dx*dx+dy*dy<global_parameters.r_cutoff_Verlet2)
			{
			global_variables.N_pairs++;
			global_variables.Verlet_list_i = (int *)realloc(global_variables.Verlet_list_i,global_variables.N_pairs*sizeof(int));
			global_variables.Verlet_list_j = (int *)realloc(global_variables.Verlet_list_j,global_variables.N_pairs*sizeof(int));
			global_variables.Verlet_list_i[global_variables.N_pairs-1] = i;
			global_variables.Verlet_list_j[global_variables.N_pairs-1] = j;
			}

		}

printf("Verlet list built. %d pairs found\n",global_variables.N_pairs);
fprintf(global_parameters.diagnostics_file,"Verlet list built. %d pairs found\n",global_variables.N_pairs);
printf("Verlet list: average neighbor/particle = %lf\n",global_variables.N_pairs/(double)global_variables.N_particles);
printf("Verlet list size [%ld bytes]\n",global_variables.N_pairs*sizeof(int));
fprintf(global_parameters.diagnostics_file,"Verlet list size [%ld bytes]\n",global_variables.N_pairs*sizeof(int));
fflush(global_parameters.diagnostics_file);
}


/* ==============================================================================================
calculate_pairwise_forces()
============================================================================================== */
void calculate_pairwise_forces()
{
int i,j,ii;
double dx,dy;
double dr2,dr;
double f;
double fx,fy;
int index;
double E;

for(ii=0;ii<global_variables.N_pairs;ii++)
	{
	i = global_variables.Verlet_list_i[ii];
	j = global_variables.Verlet_list_j[ii];

	dx = global_variables.particles[j].x - global_variables.particles[i].x;
	dy = global_variables.particles[j].y - global_variables.particles[i].y;
	//PBC fold back of distances
	if (dx>(global_parameters.SX/2.0)) 	dx = dx - global_parameters.SX;
	if (dx<(-global_parameters.SX/2.0)) dx = dx + global_parameters.SX;
	if (dy>(global_parameters.SY/2.0)) 	dy = dy - global_parameters.SY;
	if (dy<(-global_parameters.SY/2.0)) dy = dy + global_parameters.SY;

	dr2 = dx*dx+dy*dy;

	//this never happens as they stay in the non-overlapping pinning
	//if (dr2<global_parameters.r_cutoff_short2) f =

	if (dr2<global_parameters.r_cutoff_long2)
		{

		//direct calculation
		dr = sqrt(dr2);
		f = global_parameters.particle_particle_force_term/dr2/dr2;

        //energy calculations. we calculate E/B^2 as this is the quantity that
        //tells us something about how well the energy is minimized in a changing B field
        //E ~ B^2 and we want to see the change compared to that

        if (global_parameters.magnetic_field!=0.0)
            E = f*dr/3.0/global_parameters.magnetic_field/global_parameters.magnetic_field;
        else E = 0.0;

        fx = f*dx/dr;
		fy = f*dy/dr;

		/*
		dr = sqrt(dr2);
		f = global_parameters.multiplicator*exp(-global_parameters.particle_a0*dr);
		f = f*(1.0+global_parameters.particle_a0*dr)/dr2;
		fx = f*dx/dr;
		fy = f*dy/dr;
		*/

		//recall from tabulated
		/*
		index = (int)floor( (fabs(dr2)-global_parameters.r_cutoff_short2)/global_variables.step_tabulated);

		f = global_parameters.multiplicator * global_variables.f_tabulated[index];
		f = f * global_variables.particles[i].q * global_variables.particles[j].q;

		fx = f * dx;
		fy = f * dy;
		*/

		//printf("%lf %d\n",dr2,index);
		//fx = global_variables.f_tabulated[index] *dx;
		//fy = global_variables.f_tabulated[index] *dy;
		//double dr2check;
		//dr2check = global_parameters.r_cutoff_short2 + index*global_variables.step_tabulated;
		//printf("%lf %lf \n",f/dr,global_variables.f_tabulated[index]);

        global_variables.particles[i].E += E;
        global_variables.particles[j].E += E;

		global_variables.particles[i].fx -= fx;
		global_variables.particles[i].fy -= fy;

		global_variables.particles[j].fx += fx;
		global_variables.particles[j].fy += fy;
		}

	}

}

/* ==============================================================================================
void recalculate_vertex_types()
============================================================================================== */
void recalculate_vertex_types()
{
int i,j;
int neighbor_ID;

//this is where I clear all these flags
//this is error chekcing I fixed this
/*
for(i=0;i<global_variables.N_particles;i++)
    {
    global_variables.particles[i].wasclosetosomeone = 0;
    global_variables.particles[i].color = 2;
    }
*/
    
if ( (strcmp(global_parameters.particle_type,"spin_ice_particle")==0)
		&&(strcmp(global_parameters.particle_lattice_type,"in_spin_ice")==0)
		&&(strcmp(global_parameters.pinning_lattice_type,"square")==0) )
	for(i=0;i<global_variables.N_vertices;i++) recalculate_square_vertex_type(i);

// printf(">>%d\n",global_variables.N_vertices);


if ( (strcmp(global_parameters.particle_type,"spin_ice_particle")==0)
		&&(strcmp(global_parameters.particle_lattice_type,"in_spin_ice")==0)
		&&(strcmp(global_parameters.pinning_lattice_type,"hexagonal")==0) )
	for(i=0;i<global_variables.N_vertices;i++) recalculate_hexagonal_vertex_type(i);


//calculate neighbor screening
for(i=0;i<global_variables.N_vertices;i++)
    {
    global_variables.vertices[i].q_neighbor_z3 = 0;
    global_variables.vertices[i].q_neighbor_z4 = 0;
    
    for(j=0;j<global_variables.vertices[i].N_neighbor_vertices;j++)
        {
        neighbor_ID = global_variables.vertices[i].neighbor_vertex_ID[j];
        if (global_variables.vertices[neighbor_ID].z==3)
            global_variables.vertices[i].q_neighbor_z3 += global_variables.vertices[neighbor_ID].q;
        else if (global_variables.vertices[neighbor_ID].z==4)
            global_variables.vertices[i].q_neighbor_z4 += global_variables.vertices[neighbor_ID].q;
        }
    
    }

}

/* ==============================================================================================
void recalculate_square_vertex_type(int vertex_ID)
recalculates the vertex types for a square lattice
takes care of GS for the case of 2in 2 out
============================================================================================== */
void recalculate_square_vertex_type(int vertex_ID)
{
int i;
double dx,dy;
int n_close;

int i_pinningsite;
int pinningsite_direction;

n_close = 0;

for(i=0;i<global_variables.vertices[vertex_ID].N_particles_in;i++)
	{
    
	//this is the pinningsite we are in
    	i_pinningsite = global_variables.particles[global_variables.vertices[vertex_ID].particle_in_ID[i]].pinningsite_ID;
    
    //it knows its direction, cosfi is 1 for horizontal and 0 for vertical pinning sites
    pinningsite_direction = global_variables.pinningsites[i_pinningsite].cosfi;
    
    if (pinningsite_direction==1)
        {
        dx = global_variables.particles[global_variables.vertices[vertex_ID].particle_in_ID[i]].x - global_variables.vertices[vertex_ID].x;
        //PBC fold back
        if (dx>(global_parameters.SX/2.0)) 	dx = dx - global_parameters.SX;
        if (dx<(-global_parameters.SX/2.0)) dx = dx + global_parameters.SX;
        if ( (dx*dx) <= (global_parameters.pinning_lattice_ax*global_parameters.pinning_lattice_ax/4.0))
            {
            n_close++;
            //global_variables.particles[global_variables.vertices[vertex_ID].particle_in_ID[i]].wasclosetosomeone = 1;
            }
        }
    else
        {
        dy = global_variables.particles[global_variables.vertices[vertex_ID].particle_in_ID[i]].y - global_variables.vertices[vertex_ID].y;
        //PBC fold back
        if (dy>(global_parameters.SY/2.0)) 	dy = dy - global_parameters.SY;
        if (dy<(-global_parameters.SY/2.0)) dy = dy + global_parameters.SY;
        if ( (dy*dy) <= (global_parameters.pinning_lattice_ax*global_parameters.pinning_lattice_ax/4.0))
            {
            n_close++;
            //global_variables.particles[global_variables.vertices[vertex_ID].particle_in_ID[i]].wasclosetosomeone = 1;
            }
        }

	}

global_variables.vertices[vertex_ID].type = n_close;
if (n_close==2) recalculate_square_gs(vertex_ID);

//topological charge calculation
recalculate_square_topological_charge_on_vertex(vertex_ID);
}

/* ==============================================================================================
void recalculate_square_topological_charge_on_vertex(vertex_ID)
recalculates the topological charge on the vertex
this is not used in the write_statistics yet that calculates it separately
but it is used in the calcultation of the screening of topological charge
============================================================================================== */
void recalculate_square_topological_charge_on_vertex(int vertex_ID)
{

if (global_variables.vertices[vertex_ID].z==3)
switch (global_variables.vertices[vertex_ID].type)
    {
    case 0: global_variables.vertices[vertex_ID].q = -3;break;
    case 1: global_variables.vertices[vertex_ID].q = -1;break;
    case 2: global_variables.vertices[vertex_ID].q = 1;break;
    case 3: global_variables.vertices[vertex_ID].q = 3;break;
    case 4: printf("error in topological charge calculation\n");break;
    case 5: global_variables.vertices[vertex_ID].q = 1;break;
    case 6: global_variables.vertices[vertex_ID].q = 1;break;
    }
else if (global_variables.vertices[vertex_ID].z==4)
switch (global_variables.vertices[vertex_ID].type)
    {
    case 0: global_variables.vertices[vertex_ID].q = -4;break;
    case 1: global_variables.vertices[vertex_ID].q = -2;break;
    case 2: global_variables.vertices[vertex_ID].q = 0;break;
    case 3: global_variables.vertices[vertex_ID].q = 2;break;
    case 4: global_variables.vertices[vertex_ID].q = 4;break;
    case 5: global_variables.vertices[vertex_ID].q = 0;break;
    case 6: global_variables.vertices[vertex_ID].q = 0;break;
    }

}


/* ==============================================================================================
void recalculate_square_vertex_type(int vertex_ID)
if a vertex has 2 particles close if could be a GS (gorund state) vertex
this is the function to decide if it is, and if it's GS it also differentiates GS1/GS2
============================================================================================== */
void recalculate_square_gs(int vertex_ID)
{
int i;
double dx,dy;
int particle_in_given_direction[4];
int n_close;
int sum_close;

int particle_ID;
int pinningsite_ID;

for(i=0;i<4;i++)
	particle_in_given_direction[i]=0;

n_close = 0;
for(i=0;i<global_variables.vertices[vertex_ID].N_particles_in;i++)
	{
	particle_ID = global_variables.vertices[vertex_ID].particle_in_ID[i];

	dx = global_variables.particles[particle_ID].x - global_variables.vertices[vertex_ID].x;
	dy = global_variables.particles[particle_ID].y - global_variables.vertices[vertex_ID].y;
	//PBC fold back
	if (dx>(global_parameters.SX/2.0)) 	dx = dx - global_parameters.SX;
	if (dx<(-global_parameters.SX/2.0)) dx = dx + global_parameters.SX;
	if (dy>(global_parameters.SY/2.0)) 	dy = dy - global_parameters.SY;
	if (dy<(-global_parameters.SY/2.0)) dy = dy + global_parameters.SY;
	//printf("%d %lf %lf\n",vertex_ID,dx,dy);
	//find the two particles that are close
	if ( (dx*dx+dy*dy) < (global_parameters.pinning_lattice_ax*global_parameters.pinning_lattice_ax/4.0))
		{
		pinningsite_ID = global_variables.particles[particle_ID].pinningsite_ID;

		dx = global_variables.pinningsites[pinningsite_ID].x - global_variables.vertices[vertex_ID].x;
		dy = global_variables.pinningsites[pinningsite_ID].y - global_variables.vertices[vertex_ID].y;
		//PBC fold back
		if (dx>(global_parameters.SX/2.0)) 	dx = dx - global_parameters.SX;
		if (dx<(-global_parameters.SX/2.0)) dx = dx + global_parameters.SX;
		if (dy>(global_parameters.SY/2.0)) 	dy = dy - global_parameters.SY;
		if (dy<(-global_parameters.SY/2.0)) dy = dy + global_parameters.SY;

		//printf("close %lf %lf\n",dx,dy);
		//decide what direction it's in
		//it should be father than lattice_ax/2.0
		if (dx>(global_parameters.pinning_lattice_ax/3.0)) particle_in_given_direction[0]++;
		if (dx<-(global_parameters.pinning_lattice_ax/3.0)) particle_in_given_direction[2]++;
		if (dy>(global_parameters.pinning_lattice_ay/3.0)) particle_in_given_direction[1]++;
		if (dy<-(global_parameters.pinning_lattice_ay/3.0)) particle_in_given_direction[3]++;
		n_close++;
		//trying to speed up

		}
	}

sum_close = 0;
for(i=0;i<4;i++)
	sum_close += i*particle_in_given_direction[i];

//basically if chess is 0 and type is 2, then GS =5
//corresponding to this chess 1 and type 4 is also GS = 5
// the other version is chess=1 type =2 and chess=0 type=4 which are the GS=6
if (sum_close==2) global_variables.vertices[vertex_ID].type = 5 + global_variables.vertices[vertex_ID].chess_table_color;
else if (sum_close==4) global_variables.vertices[vertex_ID].type = 6 - global_variables.vertices[vertex_ID].chess_table_color;
/*
printf("recalculating %d to %d [%d %d %d %d]\n",vertex_ID,global_variables.vertices[vertex_ID].type,
	particle_in_given_direction[0],
	particle_in_given_direction[1],
	particle_in_given_direction[2],
	particle_in_given_direction[3]);
*/
}

/* ==============================================================================================
void recalculate_hexagonal_vertex_type(int vertex_ID)
============================================================================================== */
void recalculate_hexagonal_vertex_type(int vertex_ID)
{

int i;
double dx,dy;
int n_close;

n_close = 0;
for(i=0;i<global_variables.vertices[vertex_ID].N_particles_in;i++)
	{
	dx = global_variables.particles[global_variables.vertices[vertex_ID].particle_in_ID[i]].x - global_variables.vertices[vertex_ID].x;
	dy = global_variables.particles[global_variables.vertices[vertex_ID].particle_in_ID[i]].y - global_variables.vertices[vertex_ID].y;
	//PBC fold back
	if (dx>(global_parameters.SX/2.0)) 	dx = dx - global_parameters.SX;
	if (dx<(-global_parameters.SX/2.0)) dx = dx + global_parameters.SX;
	if (dy>(global_parameters.SY/2.0)) 	dy = dy - global_parameters.SY;
	if (dy<(-global_parameters.SY/2.0)) dy = dy + global_parameters.SY;
	//find how many are close
	if ( (dx*dx+dy*dy) < (global_parameters.pinning_lattice_ax*global_parameters.pinning_lattice_ax/4.0)) n_close++;
	}

global_variables.vertices[vertex_ID].type = n_close;

}


/* ==============================================================================================
void initialize_statistics()
============================================================================================== */
void initialize_statistics()
{
int i;
global_statistics.N_average = 0;

for(i=0;i<7;i++)
    {
	global_statistics.N_vertextype[i] = 0;
	global_statistics.N_vertextype3[i] = 0;
    }

//topological charge on the neighbors
global_statistics.q33 = 0;
global_statistics.q34 = 0;
global_statistics.q43 = 0;
global_statistics.q44 = 0;

global_statistics.N3_charged = 0;
global_statistics.N4_charged = 0;

global_statistics.total_energy = 0;
}

/* ==============================================================================================
void calculate_statistics_per_step()
============================================================================================== */
void calculate_statistics_per_step()
{
int i,j;
int i_neighbor;
int q_neighbor;
int z_neighbor;

recalculate_vertex_types();
for(i=0;i<global_variables.N_vertices;i++)
	{
	if (global_variables.vertices[i].z==3)
        global_statistics.N_vertextype3[global_variables.vertices[i].type]++;
	else if (global_variables.vertices[i].z==4)
        global_statistics.N_vertextype[global_variables.vertices[i].type]++;
	}

for(i=0;i<global_variables.N_particles;i++)
    global_statistics.total_energy += global_variables.particles[i].E;


//add the charge together for all vertices
for(i=0;i<global_variables.N_vertices;i++)
	{

	//only calculate screening for vertices that have a charge
	if (global_variables.vertices[i].q!=0)
		{
		
		//counting the number that is charged
		//this will replace N_average for these quantities
		if (global_variables.vertices[i].z==3) global_statistics.N3_charged++;
		if (global_variables.vertices[i].z==4) global_statistics.N4_charged++;

		//go over all the neighboring vertices
		for(j=0;j<global_variables.vertices[i].N_neighbor_vertices;j++)
			{
			i_neighbor = global_variables.vertices[i].neighbor_vertex_ID[j];
			q_neighbor = global_variables.vertices[i_neighbor].q;
			z_neighbor = global_variables.vertices[i_neighbor].z;
			//add the charge on them to the appropriate total
			if ((z_neighbor==3)&&(global_variables.vertices[i].z==3)) global_statistics.q33 += q_neighbor;
			else if ((z_neighbor==4)&&(global_variables.vertices[i].z==3)) global_statistics.q43 += q_neighbor;
			else if ((z_neighbor==3)&&(global_variables.vertices[i].z==4)) global_statistics.q34 += q_neighbor;
			else if ((z_neighbor==4)&&(global_variables.vertices[i].z==4)) global_statistics.q44 += q_neighbor;
			}
		}
	}

global_statistics.N_average++;
}

/* ==============================================================================================
void write_statistics()
write out the statistics at each step
============================================================================================== */
void write_statistics()
{
int i;
int N_hopping[2];
int total_hopping[2];
double pos1dist,pos3dist;
double posdist_temp;

//fprintf(global_parameters.statistics_file,"%d ",global_parameters.time);
fprintf(global_parameters.statistics_file,"%lf ",global_parameters.magnetic_field);

//write out the real time
//fprintf(global_parameters.statistics_file,"%lf ",global_parameters.time * global_parameters.dt);


//just the vertices in the line
//this means type 1 type 3 and type 2 vertices
/*
fprintf(global_parameters.statistics_file,"%lf %lf %lf ",(global_statistics.N_vertextype[1]+global_statistics.N_vertextype[2]+global_statistics.N_vertextype[3]) /(double)global_statistics.N_average, global_statistics.N_vertextype[1]/(double)global_statistics.N_average, global_statistics.N_vertextype[3]/(double)global_statistics.N_average);
*/

//write out the position of the 1 end and the 3 end
//this only works if there is no spontaneous nucleation
//I need it to find optimal ratchet biasing forces

//sometimes it won't find the type 3 (when it's caught in switching)
//then it just gets the old value out

//fixed it so it tolerates spont. nucleation: find smallest 1 and largest 3
/*
    pos1dist = 0.0;
    pos3dist = 0.0;

    global_statistics.pos1x = 0.0;
    global_statistics.pos1y = 0.0;
    global_statistics.pos3x = 0.0;
    global_statistics.pos3y = 0.0;


    for(i=0;i<global_variables.N_vertices;i++)
        {
        if (global_variables.vertices[i].type==1)
            {
                posdist_temp = global_variables.vertices[i].x*global_variables.vertices[i].x
                                + global_variables.vertices[i].y*global_variables.vertices[i].y;

                if ((posdist_temp<pos1dist)||(pos1dist==0.0))
                    {
                        global_statistics.pos1x = global_variables.vertices[i].x;
                        global_statistics.pos1y = global_variables.vertices[i].y;
                        pos1dist = sqrt(global_statistics.pos1x * global_statistics.pos1x
                                    + global_statistics.pos1y * global_statistics.pos1y);

                    }

            }

        if (global_variables.vertices[i].type==3)
            {

                posdist_temp = global_variables.vertices[i].x*global_variables.vertices[i].x
                                + global_variables.vertices[i].y*global_variables.vertices[i].y;

                if ((posdist_temp>pos3dist)||(pos3dist==0.0))
                {
                    global_statistics.pos3x = global_variables.vertices[i].x;
                    global_statistics.pos3y = global_variables.vertices[i].y;
                    pos3dist = sqrt(global_statistics.pos3x * global_statistics.pos3x
                    + global_statistics.pos3y * global_statistics.pos3y);

                }
            }
        }
  */

  /*
    fprintf(global_parameters.statistics_file,"%.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n",
            pos1x,pos1y,sqrt(pos1x*pos1x+pos1y*pos1y),
            pos3x,pos3y,sqrt(pos3x*pos3x+pos3y*pos3y));
  */

/* DEFECT LINE PROJECT END OF LINE STATISTICS
    //get rid of the 0 lines that happen at transition
    if ((pos1dist>100.0)&&(pos3dist>100.0))
        {
            //printf("%lf %lf\n",pos1dist,pos3dist);
            fprintf(global_parameters.statistics_file,"%lf ",global_parameters.time * global_parameters.dt);
            fprintf(global_parameters.statistics_file,"%.2lf %.2lf ",pos1dist,pos3dist);

            //printf out only the biased vertices (these measure the line length)
            fprintf(global_parameters.statistics_file,"%lf\n",
                    global_statistics.N_vertextype[2]/(double)global_statistics.N_average);

        }
*/

//all vertex type statistics

//these are vertices with 4 ends
 for(i=0;i<7;i++)
	 fprintf(global_parameters.statistics_file,"%lf ",
		 global_statistics.N_vertextype[i]/(double)global_statistics.N_average/(double)global_variables.N_vertices4);

//same thing not divided by Navg and N4
//		 for(i=0;i<7;i++)
//			fprintf(global_parameters.statistics_file,"%ld ",
//				 global_statistics.N_vertextype[i]);


//these are vertices with 3 ends
 for(i=0;i<7;i++)
	 fprintf(global_parameters.statistics_file,"%lf ",
		 global_statistics.N_vertextype3[i]/(double)global_statistics.N_average/(double)global_variables.N_vertices3);

//same thing not divided by Navg and N3
//		for(i=0;i<7;i++)
//			fprintf(global_parameters.statistics_file,"%ld ",
//				 global_statistics.N_vertextype3[i]);


int charge_type_z4[5];
int charge_type_z3[5];

charge_type_z4[0] = global_statistics.N_vertextype[2] +
	global_statistics.N_vertextype[5] +
	global_statistics.N_vertextype[6];//charge 0
charge_type_z4[1] = global_statistics.N_vertextype[1]; //charge -2
charge_type_z4[2] = global_statistics.N_vertextype[3]; //charge +2
charge_type_z4[3] = global_statistics.N_vertextype[0]; //charge -4
charge_type_z4[4] = global_statistics.N_vertextype[4]; //charge +4

charge_type_z3[0] = global_statistics.N_vertextype3[2] +
	global_statistics.N_vertextype3[5] +
	global_statistics.N_vertextype3[6];//charge +1
charge_type_z3[1] = global_statistics.N_vertextype3[1]; //charge -1
charge_type_z3[2] = global_statistics.N_vertextype3[3]; //charge +3
charge_type_z3[3] = global_statistics.N_vertextype3[0]; //charge -3
charge_type_z3[4] = global_statistics.N_vertextype3[4]; //charge this never happens


for(i=0;i<5;i++)
	fprintf(global_parameters.statistics_file,"%lf ",
        charge_type_z4[i]/(double)global_statistics.N_average/(double)global_variables.N_vertices4);

for(i=0;i<5;i++)
	fprintf(global_parameters.statistics_file,"%lf ",
		charge_type_z3[i]/(double)global_statistics.N_average/(double)global_variables.N_vertices3);

//same thing not divided by Navg and N3 or N4
//		for(i=0;i<5;i++)
//			fprintf(global_parameters.statistics_file,"%d ", charge_type_z4[i]);
//		for(i=0;i<5;i++)
//			fprintf(global_parameters.statistics_file,"%d ",charge_type_z3[i]);

int total_charge_on_z3;
int total_charge_on_z4;

total_charge_on_z4 = charge_type_z4[1] * (-2) + charge_type_z4[2] * (+2) + charge_type_z4[3] * (-4) + charge_type_z4[4]*(+4);
total_charge_on_z3 = charge_type_z3[0] * 1 + charge_type_z3[1] * (-1) + charge_type_z3[2] * 3 + charge_type_z3[3] * (-3);

fprintf(global_parameters.statistics_file,"%lf %lf",
    total_charge_on_z4/(double)global_statistics.N_average/(double)global_variables.N_vertices4,
    total_charge_on_z3/(double)global_statistics.N_average/(double)global_variables.N_vertices3);

// printf("%d %d\n",total_charge_on_z4,total_charge_on_z3);

	//fprintf(global_parameters.statistics_file,"%lf ",global_statistics.N_vertextype[i]/(double)global_statistics.N_average(double)global_variables.N_vertices);

//total energy
//fprintf(global_parameters.statistics_file," %lf",global_statistics.total_energy/(double)global_statistics.N_average);

//Q3 and Q4 for the decimated lattice project

fprintf(global_parameters.statistics_file," %lf %lf",
	global_statistics.N3_charged/(double)global_statistics.N_average/(double)global_variables.N_vertices3,
	global_statistics.N4_charged/(double)global_statistics.N_average/(double)global_variables.N_vertices4);

fprintf(global_parameters.statistics_file," %lf %lf %lf %lf",
	global_statistics.q33/(double)global_statistics.N3_charged/3.0,
	global_statistics.q43/(double)global_statistics.N3_charged/3.0,
	global_statistics.q34/(double)global_statistics.N4_charged/4.0,
	global_statistics.q44/(double)global_statistics.N4_charged/4.0);


fprintf(global_parameters.statistics_file,"\n");


//hopping statistics
/*
for(i=0;i<3;i++)
	{
	total_hopping[i] = 0;
	N_hopping[i] = 0;
	}

for(i=0;i<global_variables.N_pinningsites;i++)
	{
	//counts how many times the spins at a given distance from the double defect have hopped
	total_hopping[global_variables.pinningsites[i].spin_distance_from_dd_bin] += global_variables.pinningsites[i].spin_N_hops;
	//this counts how many spins are at a given distance from the double defect
	N_hopping[global_variables.pinningsites[i].spin_distance_from_dd_bin]++;
	}

for(i=0;i<3;i++)
	{
	//printf("hopping: %d %d %d %lf\n",i,total_hopping[i],N_hopping[i],total_hopping[i]/(double)N_hopping[i]/(double)global_statistics.N_average);
	if (N_hopping[i]!=0)
	fprintf(global_parameters.statistics_file,"%lf ",total_hopping[i]/(double)N_hopping[i]/(double)global_statistics.N_average);
	else fprintf(global_parameters.statistics_file,"%lf ",0.0);
	}

fprintf(global_parameters.statistics_file,"\n");
*/

//if (global_parameters.time==500000) writeout_hopping_file_during_run();
//if (global_parameters.time==500000) writeout_hopping_heightmap();

/*for(i=0;i<global_variables.N_pinningsites;i++)
	global_variables.pinningsites[i].spin_N_hops = 0;
*/

//zero down the variables that accumulate the data
initialize_statistics();
}


/* ==============================================================================================
void calculate_temperature()
sets the temperature according to the setpoints
============================================================================================== */
void calculate_temperature()
{
int i;
int time_setpoint_before, time_setpoint_after;
double temp_setpoint_before, temp_setpoint_after;

int N_steps_total, N_steps;
double temp_step;

double new_temp;


for(i=0;i<global_parameters.N_temperature_set_points;i++)
	{
	if ((global_parameters.time>global_parameters.set_temperature_time[i])&&
		(global_parameters.time<=global_parameters.set_temperature_time[i+1]))
			{
			time_setpoint_before = global_parameters.set_temperature_time[i];
			time_setpoint_after  = global_parameters.set_temperature_time[i+1];
			temp_setpoint_before = global_parameters.set_temperature[i];
			temp_setpoint_after  = global_parameters.set_temperature[i+1];
			break;
			}
	}

N_steps_total = (time_setpoint_after - time_setpoint_before) / global_parameters.statistics_write_time;
temp_step = (temp_setpoint_after - temp_setpoint_before) / (double)N_steps_total;

N_steps = (global_parameters.time- time_setpoint_before) / global_parameters.statistics_write_time;
new_temp = temp_setpoint_before + N_steps*temp_step;

//printf(">>%d %d %lf %lf\n",N_steps_total,N_steps,temp_step,new_temp);
//printf("%d %d %d %lf %lf\n",global_parameters.time,time_setpoint_before,time_setpoint_after,temp_setpoint_before,temp_setpoint_after);

global_parameters.temperature = new_temp;

/*
for(i=0;i<global_parameters.N_temperature_set_points;i++)
	{


	printf("%d %d %lf\n",i,global_parameters.set_temperature_time[i],global_parameters.set_temperature[i]);

	}
*/
}


/* ==============================================================================================
void calculate_external()
sets the external forcing according to the setpoints
============================================================================================== */
void calculate_external()
{
int i;
int time_setpoint_before, time_setpoint_after;
double fx_setpoint_before, fx_setpoint_after;
double fy_setpoint_before, fy_setpoint_after;

int N_steps_total, N_steps;
double fx_step,fy_step;

double new_fx,new_fy;

/*
//setpoints calculation
for(i=0;i<global_parameters.N_external_set_points;i++)
	{
	if ((global_parameters.time>global_parameters.set_external_time[i])&&
		(global_parameters.time<=global_parameters.set_external_time[i+1]))
			{
			time_setpoint_before = global_parameters.set_external_time[i];
			time_setpoint_after  = global_parameters.set_external_time[i+1];

			fx_setpoint_before = global_parameters.set_external_fx[i];
			fx_setpoint_after  = global_parameters.set_external_fx[i+1];

			fy_setpoint_before = global_parameters.set_external_fy[i];
			fy_setpoint_after  = global_parameters.set_external_fy[i+1];

			break;
			}
	}

N_steps_total = (time_setpoint_after - time_setpoint_before) / global_parameters.statistics_write_time;
fx_step = (fx_setpoint_after - fx_setpoint_before) / (double)N_steps_total;
fy_step = (fy_setpoint_after - fy_setpoint_before) / (double)N_steps_total;

N_steps = (global_parameters.time- time_setpoint_before) / global_parameters.statistics_write_time;
new_fx = fx_setpoint_before + N_steps*fx_step;
new_fy = fy_setpoint_before + N_steps*fy_step;
*/


//doing nothing just keeps the initial set biasing field
//it was set in the initialization and it is only changed in this function


//this is for the ratcheting part

/*
if ( (global_parameters.time%300000)<50000)
    {
        new_fx = 1.12;
        new_fy = 1.12;
    }
else
    {
        new_fx = 0.30;
        new_fy = 0.30;
    }

global_parameters.external_fx = new_fx;
global_parameters.external_fy = new_fy;
*/

}

/* ==============================================================================================
void calculate_spin_flipping(int pinning_i,int particle_i,double dxprime)
particle particle_i is in pinningsite_i and its calculated distance from
the center of the pinning site in the elongated direction is dxprime
dxprime>0 points up,right dxprime<0 points down,left
this helps to calculate if a spin flip occured -> stored in the pinning site

we should only get into this routine if the site is singly occupied
it does not check whether it's doubly occupied or empty
============================================================================================== */
void calculate_spin_flipping(int pinning_i,int particle_i,double dxprime)
{
int current_spin;

if (dxprime<0) current_spin=-1;
else current_spin = 1;

if (current_spin != global_variables.pinningsites[pinning_i].spin_value)
	{
	global_variables.pinningsites[pinning_i].spin_old_value = global_variables.pinningsites[pinning_i].spin_value;
	global_variables.pinningsites[pinning_i].spin_value = current_spin;
	global_variables.pinningsites[pinning_i].spin_N_hops++;
	}

}


/* ==============================================================================================
void writeout_hopping_file_during_run()
============================================================================================== */
void writeout_hopping_file_during_run()
{
FILE *f;
int i;
int min_hop,max_hop;
int n_wroteout;
int j;
double hopratio;
char filename[80];
char tempstring[80];

min_hop = global_variables.pinningsites[0].spin_N_hops;
max_hop = global_variables.pinningsites[0].spin_N_hops;
for(i=0;i<global_variables.N_pinningsites;i++)
	{
	if (global_variables.pinningsites[i].spin_N_hops>max_hop) max_hop = global_variables.pinningsites[i].spin_N_hops;
	if (global_variables.pinningsites[i].spin_N_hops<min_hop) min_hop = global_variables.pinningsites[i].spin_N_hops;
	}
printf("min hop = %d max hop = %d\n",min_hop,max_hop);
fprintf(global_parameters.diagnostics_file,"min hop = %d max hop = %d\n",min_hop,max_hop);


for(j=0;j<=10;j++)
	{
	hopratio = j/10.0;
	strcpy(filename,"results/");
	strcat(filename,global_parameters.running_name);
	strcat(filename,"_hopsite_");
	sprintf(tempstring,"%d",j*10);
	strcat(filename,tempstring);
	strcat(filename,".txt");
	f = fopen(filename,"wt");

	printf("Writing %s\n",filename);
	fprintf(global_parameters.diagnostics_file,"Writing %s\n",filename);

	n_wroteout = 0;
	for(i=0;i<global_variables.N_pinningsites;i++)
		{
		if (global_variables.pinningsites[i].spin_N_hops>max_hop*hopratio)
			{
			fprintf(f,"%lf %lf\n",global_variables.pinningsites[i].x,global_variables.pinningsites[i].y);
			n_wroteout++;
			}
		}
	printf("Wrote out %d pinningsites\n",n_wroteout);
	fprintf(global_parameters.diagnostics_file,"Wrote out %d pinningsites\n",n_wroteout);
	fclose(f);
	n_wroteout = 0;
	}

strcpy(filename,"results/");
strcat(filename,global_parameters.running_name);
strcat(filename,"_ddsites.txt");
f = fopen(filename,"wt");
printf("Writing %s\n",filename);
fprintf(global_parameters.diagnostics_file,"Writing %s\n",filename);
n_wroteout = 0;
for(i=0;i<global_variables.N_pinningsites;i++)
	{
	if (global_variables.pinningsites[i].N_particles_in==2)
		{
		fprintf(f,"%lf %lf\n",global_variables.pinningsites[i].x,global_variables.pinningsites[i].y);
		n_wroteout++;
		}
	}
printf("Wrote out %d pinningsites\n",n_wroteout);
fprintf(global_parameters.diagnostics_file,"Wrote out %d pinningsites\n",n_wroteout);
fclose(f);
fflush(global_parameters.diagnostics_file);
}


/* ==============================================================================================
void writeout_hopping_heightmap()
============================================================================================== */
void writeout_hopping_heightmap()
{
FILE *f;
char filename[80];
char tempstring[80];
int i;

strcpy(filename,"results/");
strcat(filename,global_parameters.running_name);
strcat(filename,"_hopheightmap.txt");
f = fopen(filename,"wt");
printf("Writing %s\n",filename);
fprintf(global_parameters.diagnostics_file,"Writing %s\n",filename);
for(i=0;i<global_variables.N_pinningsites;i++)
	{
	fprintf(f,"%lf %lf %d\n",global_variables.pinningsites[i].x,global_variables.pinningsites[i].y,global_variables.pinningsites[i].spin_N_hops);
	}

fclose(f);
fflush(global_parameters.diagnostics_file);
}

/* ==============================================================================================
void writeout_hophistogram()
============================================================================================== */
void writeout_hophistogram()
{
FILE *f;
char filename[80];
char tempstring[80];
int i;

strcpy(filename,"results/");
strcat(filename,global_parameters.running_name);
strcat(filename,"_hophistogram.txt");

if (global_parameters.time==0)
	{
	f = fopen(filename,"wt");
	printf("Created %s\n",filename);
	fprintf(global_parameters.diagnostics_file,"Created %s\n",filename);
	}
else f = fopen(filename,"at");


for(i=0;i<global_variables.N_pinningsites;i++)
	{
	fprintf(f,"%d\n",global_variables.pinningsites[i].spin_N_hops);
	}

fclose(f);
fflush(global_parameters.diagnostics_file);
}

/* ==============================================================================================
 void calculate_magnetic_field()
 sets the external magnetic field according to the protocol we want to use
 for example for the Kibbel-Zurek I want to magnetize very slowly
 ============================================================================================== */
void calculate_magnetic_field()
{
int time_period;

//current time:             global_parameters.time
//total running time:       global_parameters.total_runtime
//initial magnetic field:   global_parameters.initial_magnetic_field
//magnetic field:   global_parameters.magnetic_field

    //this is to ramp the field once up
    //global_parameters.magnetic_field = global_parameters.initial_magnetic_field;

    global_parameters.magnetic_field = (double)global_parameters.time/(double)global_parameters.total_runtime * global_parameters.initial_magnetic_field;


    //cycling the field multiple times
    /*
    time_period = 1000000;

    if (global_parameters.time%time_period < time_period/2)
         global_parameters.magnetic_field = (double)( (time_period/2) - global_parameters.time%time_period) / (time_period/2.0) * global_parameters.initial_magnetic_field;
    else global_parameters.magnetic_field = (double)(global_parameters.time%time_period - (time_period/2)) / (time_period/2.0) * global_parameters.initial_magnetic_field;

    */



    //redo the related calculations

    global_parameters.particle_magnetization = 1000.0 * global_parameters.magnetic_susceptibility * global_parameters.particle_volume * global_parameters.magnetic_field / global_parameters.magnetic_permeability;

    global_parameters.particle_particle_force_term = global_parameters.particle_particle_force_pre_term * global_parameters.particle_magnetization * global_parameters.particle_magnetization;

}
