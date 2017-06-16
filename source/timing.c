

#include <stdio.h>
#include "timing.h"
#include "global_variables.h"

#ifdef __MACH__
#include <sys/time.h>
//clock_gettime is not implemented on OSX
int clock_gettime(int clk_id, struct timespec* t) {
    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    if (rv) return rv;
    t->tv_sec  = now.tv_sec;
    t->tv_nsec = now.tv_usec * 1000;
    return 0;
}
#define CLOCK_REALTIME 0 
#define CLOCK_MONOTONIC 0
#endif

struct timespec program_start_time;
struct timespec program_finish_time;

/* ============================================================================================== 
void get_starting_time()
============================================================================================== */
void get_starting_time()
{
clock_gettime(CLOCK_MONOTONIC, &program_start_time);
}

/* ============================================================================================== 
void get_finishing_time()
============================================================================================== */
void get_finishing_time()
{
clock_gettime(CLOCK_MONOTONIC, &program_finish_time);
}

/* ============================================================================================== 
void echo_running_time()
============================================================================================== */
void echo_running_time()
{
double time_expired;
time_t time_endtime;
struct tm *timeinfo;

time_expired  = program_finish_time.tv_sec - program_start_time.tv_sec + (program_finish_time.tv_nsec - program_start_time.tv_nsec)*1e-9;
time(&time_endtime);
timeinfo = localtime(&time_endtime);

printf("==============================================================================================\n");
printf("Run %s finished\n",global_parameters.running_name);
printf("Run finished on: %s",asctime(timeinfo));
printf("Program running time was: %.2lf seconds\n",time_expired);
printf("==============================================================================================\n");

fprintf(global_parameters.diagnostics_file,"==============================================================================================\n");
fprintf(global_parameters.diagnostics_file,"Run %s finished\n",global_parameters.running_name);
fprintf(global_parameters.diagnostics_file,"Run finished on: %s",asctime(timeinfo));
fprintf(global_parameters.diagnostics_file,"Program running time was: %.2lf seconds\n",time_expired);
fprintf(global_parameters.diagnostics_file,"==============================================================================================\n");

fflush(global_parameters.diagnostics_file);
}

/* ============================================================================================== 
void estimate_running_time()
============================================================================================== */
void estimate_running_time()
{
double time_expired;
struct timespec program_current_time;

double time_seconds_to_add;
double time_running_time_in_seconds;

time_t time_estimated_endtime;
struct tm *timeinfo;


clock_gettime(CLOCK_MONOTONIC, &program_current_time);
time(&time_estimated_endtime);

time_expired  = program_current_time.tv_sec - program_start_time.tv_sec + (program_current_time.tv_nsec - program_start_time.tv_nsec)*1e-9;

if (global_parameters.time!=0)
	{
	time_seconds_to_add = time_expired * (global_parameters.total_runtime - global_parameters.time) / (double)global_parameters.time;
	time_running_time_in_seconds = time_expired * global_parameters.total_runtime / (double)global_parameters.time;
	}
else time_seconds_to_add = 0.0;

time_estimated_endtime = time_estimated_endtime + time_seconds_to_add;
timeinfo = localtime(&time_estimated_endtime);
printf("\n");
printf("Program estimated total runtime: %.2lf seconds\n",time_running_time_in_seconds);
printf("Program estimated to end at: %s",asctime(timeinfo));
printf("\n");

fprintf(global_parameters.diagnostics_file,"\n");
fprintf(global_parameters.diagnostics_file,"Program estimated total runtime: %.2lf seconds\n",time_running_time_in_seconds);
fprintf(global_parameters.diagnostics_file,"Program estimated to end at: %s",asctime(timeinfo));
fprintf(global_parameters.diagnostics_file,"\n");
fflush(global_parameters.diagnostics_file);
}
