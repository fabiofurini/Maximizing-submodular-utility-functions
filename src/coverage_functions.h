#ifndef COVERAGE_FUNCTIONS_HEADER
#define COVERAGE_FUNCTIONS_HEADER

#include "global_variables.h"

#include "global_functions.h"

/*****************************************************************/
double compute_val_diff(instance *inst, double p);
/*****************************************************************/

/*****************************************************************/
double compute_val_funct(instance *inst, double p);
/*****************************************************************/

/*****************************************************************/
double compute_subset_coverage_utility(instance *inst,double *y);
/*****************************************************************/

/*****************************************************************/
double compute_subset_coverage_utility_mod_scenario(instance *inst,double *y,int i);
/*****************************************************************/

/*****************************************************************/
int compute_subset_coverage_number(instance *inst,double *y,bool *reached_items);
/*****************************************************************/

/*****************************************************************/
double compute_subset_coverage_utility_scenario(instance *inst,double *y,int i);
/*****************************************************************/

/*****************************************************************/
double compute_subset_coverage_utility_scenario_fract(instance *inst,double *y,int i);
/*****************************************************************/


#endif
