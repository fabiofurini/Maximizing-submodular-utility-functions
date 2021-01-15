#ifndef VARIABLE_local_HEADER
#define VARIABLE_local_HEADER


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sstream>
#include <vector>
#include <algorithm>

// normal_distribution
#include <iostream>
#include <string>
#include <random>

#include </home/fabio/ILOG/CPLEX_Studio_AcademicResearch129/cplex/include/ilcplex/cplex.h>
//#include <cplex.h>

using namespace std;

////////////////////////////////////GENERAL///////////////////////////////////
#define MAX_CONNECTIONS 100000000 //maximum number of arcs
//////////////////////////////////////////////////////////////////////////////

//#define PRINT_FIGURE_TEX

//#define PRINT_SINGLE_FILE_OUTPUT

typedef struct
{



	///////////////////////////////INPUT//////////////////////////////////
	char *input_file;
	char *param_file;
	int algorithm;
	int option;
	double timelimit;
	int n_items;
	int m_scenarios;
	int n_meta_items;
	double lambda; //parameter for the  utility function
	double delta;
	int seed;
	double TOLL_OPTIMALITY;
	double ProbScenario;//probability of an item in one scenario
	int  cardinality;//value of the cardinality RHS
	double conflict_perc;//probability of a conflict between meta-items
	double value_a;//it determines the value of the a_{ij} (or its upper bound)
	bool distribute_a;//set to 1, to have the a_{ij} randomly distributed in the interval [0,value_a] with integer values

	int KP_constraint;//set to > 0 to insert the KP constraint instead of the cardinality constraint

	//the profit for the generation are the marginal contribution to the empty set (global contribution -> rho single)

//	1 -> Uncorrelated: w j u.r. in [1, R].
//	2 -> Weakly correlated:  w j u.r. in [max{1, p j − R/10}, p j + R/10].
//	3 -> Strongly correlated:  w j = p j + R/10.
//	4 -> Almost strongly correlated: w j u.r. in [p j + R/10 − R/500, p j + R/10 + R/500].
//	5 -> Subset-sum:  w j = p j.

	//the profits are scaled as follow (int)(CONSTANT_SCALING*pp[i])

	int KP_constraint_R_VALUE;
	double KP_constraint_perc_cap;//it defines the capacity as the percentage of total item weight
	double KP_constraint_CAPACITY;
	double *KP_constraint_weights;//weights of the meta-items

	int USE_WORST_CASE_INSTANCE;//1 to create the worst case instances and 0 otherwise (the only parameter is cardinality to create these instances)

	int type_of_zed_function;

	//1 too use the exponential utility function with lambda
	//2 identity function

	int partition_constraints;//1 to activate partition constraints
	int meta_item_per_element; //number of meta items per element of the partition (they are taken in order) (the last element can have fewer meta-item!)
	int budget_per_element;//cardinality (budget) for each element of the partition


	int TEST_ID;

	double scale_factor_alpha;//it is used to determine alpha in case of non monotone utility functions
	//////////////////////////////////////////////////////////////////////


	double alpha;//use for the non monotone utility function

	//////////////////////////////////////////////////////////////////////
	//computed according to the partition input parameters
	int n_element_partition;//it tells how elements there are
	int *meta_item_element_partition;//for each meta-item it tells in which element of the partition is
	//////////////////////////////////////////////////////////////////////

	double lambda_orig;

	double curvature;

	double curvature_scenario_avg;
	double curvature_scenario_min;
	double curvature_scenario_max;

	////////////////////////////////PARAMETER/////////////////////////////
	int number_of_CPU;
	double TOLL_DERIVATIVE;
	double TOLL_VIOL;
	double TOLL_VIOL_FRAC_BEN;
	//////////////////////////////////////////////////////////////////////

	///////////////////////////////OBJECTIVE FUNCTION AND CONSTRAINTS////////////
	double **a; //value item per scenario
	double *weight_item_MP; //capital requirements of item (MP INSTANCES)
	/////////////////////////////////////////////////////////////////////////////


	double AVERAGE_DEMAND;

	/////////////////////////////////INSTANCE TYPE///////////////////////////////
	bool FLAG_INSTANCE_LOCATION;
	bool FLAG_INSTANCE_MP;
	bool FLAG_INSTANCE_WORST_CASE;
	int **CONF_MATRIX;//it store the conflicts between the meta-items
	int n_conflicts;
	/////////////////////////////////////////////////////////////////////////////



	/////////////////////////////////////DATA COVERAGE/////////////////////////////

	///////////////////////////////////////////////////////////////////////////////
	int *NFS;     // NFS[i]: begin of FS(i) in AFS
	int *AFS;     // item index for FS
	int *DP;  // number of neighbours delta+

	int *NBS;     // NBS[i]: begin of BS(i) in ABS
	int *ABS;     // meta item index  for BS
	int *DM;  // number of neighbours delta-
	///////////////////////////////////////////////////////////////////////////////

	double *x_location;
	double *y_location;

	double *x_client;
	double  *y_client;

	int counter_c;//number of backward arcs
	int counter_l;//number of farward arcs

	int item_not_covered;
	int item_single_covered;

	bool *item_OK;//boolean vector checking if an item is covered by at least one meta-item

	double *reached_item;//used in the coverage routines
	////////////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////CPLEX/////////////////////////////////////

	CPXENVptr env_MOD;
	CPXLPptr  lp_MOD;

	CPXENVptr env_MOD_OUTER;
	CPXLPptr  lp_MOD_OUTER;

	CPXENVptr env_MOD_LOCAL;
	CPXLPptr  lp_MOD_LOCAL;

	CPXENVptr env_BEN;
	CPXLPptr  lp_BEN;


	int status,ccnt,rcnt,nzcnt;
	int* rmatbeg, *rmatind,*cmatbeg, *cmatind;
	double* rmatval,*cmatval,*x,*pi,*obj, *lb, *ub,*rhs;
	char *c_type,* sense;
	char **colname;
	double objval,bestobjval;
	int lpstat,nodecount;
	///////////////////////////////////////////////////////////////////////////////


	int n_cuts_MOD_LOWER;//number of cuts submodular lower bound


	///////////////////////////////////////////////////////////////////////////////
	int n_cuts_BEN_1;
	int n_cuts_BEN_FRAC_1;
	int n_cuts_BEN_2;
	int n_cuts_BEN_FRAC_2;

	double _cut_BEN_obj_value;
	double _cut_BEN_RHS;
	double _cut_BEN_value;

	double *_cut_BEN_rmatval;
	double *_cut_BEN_Y;
	int    *_cut_BEN_rmatind;
	double *_cut_BEN_local_rho;

	double **_cut_BEN_super_rho;
	double **_cut_BEN_single_rho;

	double *AUX_SOL_BEN;
	double *I_TILDE_BEN;
	///////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////
	double *AUX_SOL_MOD_LOCAL;
	double *I_TILDE_MOD_LOCAL;

	int n_cuts_MOD_LOCAL_1;
	int n_cuts_MOD_LOCAL_FRAC_1;
	int n_cuts_MOD_LOCAL_2;
	int n_cuts_MOD_LOCAL_FRAC_2;


	double _cut_MOD_LOCAL_obj_value;
	double _cut_MOD_LOCAL_RHS;
	double _cut_MOD_LOCAL_value;

	double *_cut_MOD_LOCAL_rmatval;
	double *_cut_MOD_LOCAL_Y;
	int    *_cut_MOD_LOCAL_rmatind;
	double *_cut_MOD_LOCAL_local_rho;

	double **_cut_MOD_LOCAL_super_rho;
	double **_cut_MOD_LOCAL_single_rho;
	///////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////
	int n_cuts_MOD_OUTER_1;
	int n_cuts_MOD_OUTER_FRAC_1;
	int n_cuts_MOD_OUTER_2;
	int n_cuts_MOD_OUTER_FRAC_2;

	double _cut_MOD_OUTER_obj_value;
	double _cut_MOD_OUTER_RHS;
	double _cut_MOD_OUTER_value;

	double *_cut_MOD_OUTER_rmatval;
	double *_cut_MOD_OUTER_Y;
	int    *_cut_MOD_OUTER_rmatind;
	double *_cut_MOD_OUTER_local_rho;

	double **_cut_MOD_OUTER_super_rho;
	double **_cut_MOD_OUTER_single_rho;

	double *AUX_SOL_MOD_OUTER;
	double *I_TILDE_MOD_OUTER;
	///////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////
	double *_cut_rmatval;
	double *_cut_Y;
	double *_cut_SET;

	int    *_cut_rmatind;

	double *_cut_local_rho;
	double *_cut_super_rho;
	double *_cut_single_rho;

	double _cut_obj_value;
	double _cut_RHS;
	double _cut_value;

	int n_cuts_MOD_1;
	int n_cuts_MOD_2;
	///////////////////////////////////////////////////////////////////////////////

} instance;


#endif
