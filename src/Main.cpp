/*
 *		Main.cpp
 *		Created on: 01/12/2017
 *		Author: Fabio Furini
 */


#include "MOD_GLOBAL.h"
#include "MOD_OUTER.h"
#include "MOD_LOCAL.h"

#include "BEN_OUTER.h"


#include "global_functions.h"
#include "global_variables.h"

#include "GREEDY_ALG.h"

#include "instance_reader.h"
#include "instance_generator.h"

//#define COMPUTE_CORVATURE

/***************************************************************************/
int main(int argc, char** argv)
/***************************************************************************/
{

	instance inst;

	inst.input_file = (char *) calloc(1000, sizeof(char));
	inst.param_file = (char *) calloc(1000, sizeof(char));

	////////////////////////////////////////////////////////////////////////////////////////
	if (argc == 28)
	{
		/*Param1*/strcpy(inst.input_file, argv[1]);
		/*Param2*/strcpy(inst.param_file, argv[2]);
		/*Param3*/inst.algorithm=atoi(argv[3]);
		/*Param4*/inst.option=atoi(argv[4]);
		/*Param5*/inst.timelimit=atof(argv[5]);
		/*Param6*/inst.n_items=atoi(argv[6]);
		/*Param7*/inst.m_scenarios=atoi(argv[7]);
		/*Param8*/inst.n_meta_items=atoi(argv[8]);
		/*Param9*/inst.lambda=atof(argv[9]);
		/*Param10*/inst.delta=atof(argv[10]);
		/*Param11*/inst.seed=atoi(argv[11]);
		/*Param12*/inst.TOLL_OPTIMALITY=atof(argv[12]);
		/*Param13*/inst.ProbScenario=atof(argv[13]);
		/*Param14*/inst.cardinality=atoi(argv[14]);
		/*Param15*/inst.conflict_perc=atof(argv[15]);

		/*Param16*/inst.value_a=atof(argv[16]);
		/*Param17*/inst.distribute_a=atoi(argv[17]);


		/*Param17*/inst.KP_constraint=atoi(argv[18]);;
		/*Param18*/inst.KP_constraint_R_VALUE=atoi(argv[19]);;
		/*Param19*/inst.KP_constraint_perc_cap=atof(argv[20]);;

		/*Param20*/inst.USE_WORST_CASE_INSTANCE=atoi(argv[21]);;
		/*Param21*/inst.type_of_zed_function=atoi(argv[22]);;


		/*Param22*/inst.partition_constraints=atoi(argv[23]);;
		/*Param23*/inst.meta_item_per_element=atoi(argv[24]);;
		/*Param24*/inst.budget_per_element=atoi(argv[25]);;

		/*Param25*/inst.TEST_ID=atoi(argv[26]);;

		/*Param26*/inst.scale_factor_alpha=atof(argv[27]);;


		srand(inst.seed);
	}
	else
	{
		cout << "ERROR NUMBER STANDARD PARAMETERS" << endl;
		cout << "Param1:\t  instance name (NULL)\n";
		cout << "Param2:\t param file\n";
		cout << "Param3:\t  algorithm\n";
		cout << "Param4:\t  option\n";
		cout << "Param5:\t  time limit\n";
		cout << "Param6:\t  n_items\n";
		cout << "Param7:\t  m_scenarios\n";
		cout << "Param8:\t  n_meta_items\n";
		cout << "Param9:\t  lambda\n";
		cout << "Param10:\t  delta\n";
		cout << "Param11:\t seed\n";
		cout << "Param12:\t TOLL_OPTIMALITY\n";
		cout << "Param13:\t ProbScenario\n";
		cout << "Param14:\t cardinality\n";
		cout << "Param15:\t conflict_perc\n";
		cout << "Param16:\t value_a\n";
		cout << "Param17:\t distribute_a\n";

		cout << "Param18:\t KP_constraint\n";
		cout << "Param19:\t KP_constraint_R_VALUE\n";
		cout << "Param20:\t KP_constraint_perc_cap\n";

		cout << "Param21:\t USE_WORST_CASE_INSTANCE\n";
		cout << "Param22:\t type_of_zed_function\n";

		cout << "Param23:\t partition_constraints\n";
		cout << "Param24:\t meta_item_per_element\n";
		cout << "Param25:\t budget_per_element\n";

		cout << "Param26:\t TEST_ID\n";

		cout << "Param27:\t scale_factor_alpha\n";

		exit(-1);
	}
	////////////////////////////////////////////////////////////////////////////////////////


	cout << "\n*****************************\n";
	cout << "Param1:\tinput_file\t" << inst.input_file << endl;
	cout << "Param2:\tparam_file\t" << inst.param_file << endl;
	cout << "Param3:\talgorithm\t" << inst.algorithm<< endl;
	cout << "Param4:\toption\t" << inst.option<< endl;
	cout << "Param5:\ttimelimit\t" <<  inst.timelimit << endl;
	cout << "Param6:\tn_items\t" <<  inst.n_items << endl;
	cout << "Param7:\tm_scenarios\t" <<  inst.m_scenarios << endl;
	cout << "Param8:\tn_meta_items\t" <<  inst.n_meta_items << endl;
	cout << "Param9:\tlambda\t" <<  inst.lambda << endl;
	cout << "Param10:\tdelta\t" <<  inst.delta << endl;
	cout << "Param11:\tseed\t" << inst.seed << endl;
	cout << "Param12:\tTOLL_OPTIMALITY\t" << inst.TOLL_OPTIMALITY << endl;
	cout << "Param13:\tProbScenario\t" << inst.ProbScenario << endl;
	cout << "Param14:\tcardinality\t" << inst.cardinality << endl;
	cout << "Param15:\tconflict_perc\t" << inst.conflict_perc << endl;

	cout << "Param16:\tvalue_a\t" << inst.value_a << endl;
	cout << "Param17:\tdistribute_a\t" << inst.distribute_a << endl;

	cout << "Param18:\tKP_constraint\t" << inst.KP_constraint << endl;
	cout << "Param19:\tKP_constraint_R_VALUE\t" << inst.KP_constraint_R_VALUE << endl;
	cout << "Param20:\tKP_constraint_perc_cap\t" << inst.KP_constraint_perc_cap << endl;

	cout << "Param21:\tUSE_WORST_CASE_INSTANCE\t" << inst.USE_WORST_CASE_INSTANCE << endl;
	cout << "Param22:\ttype_of_zed_function\t" << inst.type_of_zed_function << endl;

	cout << "Param23:\tpartition_constraints\t" << inst.partition_constraints << endl;
	cout << "Param24:\tmeta_item_per_element\t" << inst.meta_item_per_element << endl;
	cout << "Param25:\tbudget_per_element\t" << inst.budget_per_element << endl;


	cout << "Param26:\tTEST_ID\t" << inst.TEST_ID << endl;

	cout << "Param27:\tscale_factor_alpha\t" << inst.scale_factor_alpha << endl;


	cout << "\n*****************************\n";


	////////////////////////////////////////////////////////
	//function to create the instance files of the MP paper
	//generate_data_MP(&inst);
	////////////////////////////////////////////////////////

	inst.FLAG_INSTANCE_LOCATION=false;
	inst.FLAG_INSTANCE_MP=false;
	inst.FLAG_INSTANCE_WORST_CASE=false;

	//	////////////////////////////////////////////////////////////////////////////////////////
	if(inst.USE_WORST_CASE_INSTANCE==1)
	{
		generate_WORT_CASE_INSTANCES(&inst,inst.cardinality);

		inst.FLAG_INSTANCE_WORST_CASE=true;
	}
	else
	{

		ifstream in(inst.input_file);
		if(in)
		{

			cout << "READING MP INSTANCES....\n";

			inst.FLAG_INSTANCE_MP=true;

			read_file_MP(&inst);

		}
		else
		{

			cout << "GENERATING LOCATION INSTANCE....\n";

			inst.FLAG_INSTANCE_LOCATION=true;

			generate_data_LOCATION(&inst);

		}

		in.close();
	}
	//	////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////
	read_param_file(&inst);

	cout << "\n*****************************\n";
	cout << "NUMBER OF CPUs\t" << inst.number_of_CPU << endl;
	cout << "TOLL_DERIVATIVE\t" << setprecision(9) << setw(9) << inst.TOLL_DERIVATIVE << endl;
	cout << "TOLL_VIOL\t" << setprecision(9) << setw(9) << inst.TOLL_VIOL   << endl;
	cout << "TOLL_VIOL_FRAC_BEN\t" << setprecision(9) << setw(9)  << inst.TOLL_VIOL_FRAC_BEN << endl;
	cout << "FLAG_GREEDY_SOLs\t" << inst.FLAG_GREEDY_SOL << endl;
	cout << "*****************************\n";
	////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////
	init_data(&inst);
	//////////////////////

	///////////////////////////////////////
	inst.n_conflicts=load_conflicts(&inst);
	cout << "n_conflicts\t" << inst.n_conflicts << endl;
	///////////////////////////////////////

	///////////////////////////////////////
	if(inst.KP_constraint>0)
	{
		cout << "\n\n+++LOAD WEIGHTS FOR THE KP CONSTRAINTS+++\n\n";
		load_weights(&inst);
	}
	///////////////////////////////////////



	///////////////////////////////////////
	if(inst.partition_constraints>0)
	{
		cout << "\n\n+++LOAD PARTITION+++\n\n";
		load_partition(&inst);
	}
	///////////////////////////////////////


	//////////////////////////////////////
	if(inst.FLAG_INSTANCE_LOCATION==true ||inst.FLAG_INSTANCE_WORST_CASE==true)
	{
		set_lambda_and_alpha(&inst);
	}
	else
	{
		inst.alpha = inst.scale_factor_alpha;
		cout << "lambda\t" << inst.lambda << endl;
		cout << "alpha\t" << inst.alpha << endl;
	}

	//////////////////////////////////////


	//	write_data_dat_file(&inst);
	//	cout << "PRINT FILE!!!!\n\n";
	//	exit(-1);

	/////////////////////////////////
#ifdef COMPUTE_CORVATURE
	compute_curvature(&inst);
	exit(-1);
#endif
	/////////////////////////////////

	if(inst.algorithm==1)
	{

		if(inst.FLAG_INSTANCE_LOCATION==true || inst.FLAG_INSTANCE_WORST_CASE==true)
		{
			cout << "\n\n----------->>>>>>GREEDY\n";

			if(inst.KP_constraint>0)
			{
				greedy_algorithm_KP_CONSTRAINT(&inst);
			}

			if(inst.cardinality!=0 && inst.partition_constraints==0)
			{
				greedy_algorithm_CARDINALITY(&inst);
			}

			if(inst.partition_constraints!=0)
			{
				greedy_algorithm_PARTITION(&inst);
			}
		}

	}


	if(inst.algorithm==11)
	{


		if(inst.type_of_zed_function==3)
		{
			cout << "METHOD NOT VALID FOR THESE UTILITY FUNCTIONs...\n";
			exit(-1);
		}

		cout << "\n\n----------->>>>>>MOD GLOBAL\n";

		///////////////////////////
		build_model_MOD(&inst);
		solve_model_MOD(&inst);
		clean_model_MOD(&inst);
		///////////////////////////

	}

	if(inst.algorithm==12)
	{

		if(inst.type_of_zed_function==3)
		{
			cout << "METHOD NOT VALID FOR THESE UTILITY FUNCTIONs...\n";
			exit(-1);
		}

		cout << "\n\n----------->>>>>>MOD LOCAL\n";

		///////////////////////////
		build_model_MOD_LOCAL(&inst);
		solve_model_MOD_LOCAL(&inst);
		clean_model_MOD_LOCAL(&inst);
		///////////////////////////

	}

	if(inst.algorithm==13)
	{

		cout << "\n\n----------->>>>>>MOD OUTER\n";

		///////////////////////////
		build_model_MOD_OUTER(&inst);
		solve_model_MOD_OUTER(&inst);
		clean_model_MOD_OUTER(&inst);
		///////////////////////////

	}

	if(inst.algorithm==14)
	{

		cout << "\n\n----------->>>>>>OUTER BENDER\n";

		///////////////////////////
		build_model_BEN(&inst);
		solve_model_BEN(&inst);
		clean_model_BEN(&inst);
		///////////////////////////

	}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	free(inst.input_file);
	free(inst.param_file);

	//////////////////////
	free_data(&inst);
	//////////////////////


	printf("\nDONE!");

	return 1;
}



