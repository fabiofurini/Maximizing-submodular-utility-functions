

#include "GREEDY_ALG.h"


#define EPSILON_TIES 1e-9

/*****************************************************************/
void greedy_algorithm_CARDINALITY(instance *inst)
/*****************************************************************/
{

	clock_t time_start=clock();

	double *CURRENT_SOL=new double[inst->n_meta_items];
	for ( int j = 0; j < inst->n_meta_items; j++){CURRENT_SOL[j]=0;}

	cout << "cardinality\t" << inst->cardinality << endl;


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for ( int k = 0; k < inst->n_meta_items; k++)
	{
		if( k == inst->cardinality)
		{
			cout << "**BREAK**\n";
			break;
		}

		cout << "iter\t" << k << endl;

		int best_j=-1;
		double best_val=-1;

		double val_1=compute_subset_coverage_utility(inst,CURRENT_SOL);

		for ( int j = 0; j < inst->n_meta_items; j++)
		{

			//////////////////////////////////////////////////////////
			bool OK_CONFLITS=true;
			if(inst->n_conflicts>0)
			{
				for ( int i = 0; i < inst->n_meta_items; i++)
				{
					if(i==j){continue;}

					if(CURRENT_SOL[i]>0.5 && inst->CONF_MATRIX[i][j]==1)
					{

						cout << "**CONFLICT**\t" << i << "\t" << j << endl;

						OK_CONFLITS=false;
						break;
					}
				}
			}
			//////////////////////////////////////////////////////////

			if(CURRENT_SOL[j]<0.5 && OK_CONFLITS)
			{

				CURRENT_SOL[j]=1;

				double val_2=compute_subset_coverage_utility(inst,CURRENT_SOL);

				CURRENT_SOL[j]=0;

				//cout << "marginal\t" << (val_2 - val_1)  << endl;

				if(best_val < (val_2 - val_1) - EPSILON_TIES )
				{
					best_val = (val_2 - val_1);
					best_j=j;

					//cout << "UPDATE\t" << best_j << endl;
				}
			}
		}

		if(best_j!=-1)
		{
			cout << "best_j\t" << best_j << "\t" << best_val << endl;
			CURRENT_SOL[best_j]=1;
			//cin.get();
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int item_counter=0;
	for ( int j = 0; j < inst->n_meta_items; j++)
	{
		cout << "METAITEM\t" << j << "\t->\t"<< CURRENT_SOL[j] << endl;
		item_counter+=CURRENT_SOL[j];
	}

	double val_sol=compute_subset_coverage_utility(inst,CURRENT_SOL);

	cout << "val_sol\t" << val_sol << "\t item_counter \t" << item_counter << endl;

	delete []CURRENT_SOL;

	clock_t time_end=clock();
	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;


	//////////////////////////////////////////////////////////////////////////////////
	ofstream compact_file;
	compact_file.open("info_GREEDY_CARDINALITY.txt", ios::app);
	compact_file << fixed

			<<  val_sol << "\t"
			<<  item_counter << "\t"
			<<  solution_time << "\t"

			<<  inst->input_file<< "\t"
			<<  inst->param_file<< "\t"
			<<  inst->algorithm<< "\t"
			<<  inst->option<< "\t"
			<<  inst->timelimit<< "\t"
			<<  inst->n_items<< "\t"
			<<  inst->m_scenarios<< "\t"
			<<  inst->n_meta_items<< "\t"
			<<  inst->lambda<< "\t"
			<<  inst->delta<< "\t"
			<<  inst->seed<< "\t"
			<<  inst->TOLL_OPTIMALITY<< "\t"
			<<  inst->ProbScenario<< "\t"
			<<  inst->cardinality<< "\t"
			<<  inst->conflict_perc<< "\t"

			<<   inst->value_a	<< "\t"
			<<   inst->distribute_a << "\t"

			<< inst->KP_constraint << "\t"
			<< inst->KP_constraint_R_VALUE << "\t"
			<< inst->KP_constraint_perc_cap << "\t"
			<< inst->KP_constraint_CAPACITY << "\t"

			<< inst->USE_WORST_CASE_INSTANCE << "\t"
			<< inst->type_of_zed_function << "\t"

			<< inst->partition_constraints << "\t"
			<< inst->meta_item_per_element << "\t"
			<< inst->budget_per_element << "\t"

			<<  inst->TEST_ID << "\t"

			<<   inst->scale_factor_alpha << "\t"

			<< endl;
	compact_file.close();
	//////////////////////////////////////////////////////////////////////////////////

}


/*****************************************************************/
void greedy_algorithm_PARTITION(instance *inst)
/*****************************************************************/
{

	clock_t time_start=clock();

	double *ITEM_ELEMENT=new double[inst->n_element_partition];
	for ( int j = 0; j < inst->n_element_partition; j++){ITEM_ELEMENT[j]=0;}

	double *CURRENT_SOL=new double[inst->n_meta_items];
	for ( int j = 0; j < inst->n_meta_items; j++){CURRENT_SOL[j]=0;}

	int total_budget=inst->budget_per_element*inst->n_element_partition;

	cout << "budget element\t" << inst->budget_per_element << "\t total_budget \t" << total_budget<<endl;

	int item_counter=0;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for ( int k = 0; k < inst->n_meta_items; k++)
	{
		if( item_counter == total_budget)
		{
			cout << "**BREAK**\n";
			break;
		}

		cout << "iter\t" << k << "\t item_counter \t " << item_counter << endl;

		int best_j=-1;
		double best_val=-1;

		double val_1=compute_subset_coverage_utility(inst,CURRENT_SOL);

		for ( int j = 0; j < inst->n_meta_items; j++)
		{

			//////////////////////////////////////////////////////////
			bool OK_CONFLITS=true;
			if(inst->n_conflicts>0)
			{
				for ( int i = 0; i < inst->n_meta_items; i++)
				{
					if(i==j){continue;}

					if(CURRENT_SOL[i]>0.5 && inst->CONF_MATRIX[i][j]==1)
					{

						cout << "**CONFLICT**\t" << i << "\t" << j << endl;

						OK_CONFLITS=false;
						break;
					}
				}
			}
			//////////////////////////////////////////////////////////

			//////////////////////////////////////////////////////////
			bool OK_PARTITION=true;
			if( ITEM_ELEMENT[inst->meta_item_element_partition[j]] == inst->budget_per_element)
			{
				OK_PARTITION=false;
			}
			//////////////////////////////////////////////////////////


			if(CURRENT_SOL[j]<0.5 && OK_CONFLITS && OK_PARTITION)
			{

				CURRENT_SOL[j]=1;

				double val_2=compute_subset_coverage_utility(inst,CURRENT_SOL);

				CURRENT_SOL[j]=0;

				//cout << "marginal\t" << (val_2 - val_1)  << endl;

				if(best_val < (val_2 - val_1) - EPSILON_TIES )
				{
					best_val = (val_2 - val_1);
					best_j=j;

					//cout << "UPDATE\t" << best_j << endl;
				}
			}
		}

		if(best_j!=-1)
		{
			cout << "best_j\t" << best_j << "\t" << best_val << "\t element \t" << inst->meta_item_element_partition[best_j] << endl;
			CURRENT_SOL[best_j]=1;
			item_counter++;

			ITEM_ELEMENT[inst->meta_item_element_partition[best_j]]=ITEM_ELEMENT[inst->meta_item_element_partition[best_j]]+1;

			//cin.get();
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	item_counter=0;
	for ( int j = 0; j < inst->n_meta_items; j++)
	{
		cout << "METAITEM\t" << j << "\t->\t"<< CURRENT_SOL[j] << endl;
		item_counter+=CURRENT_SOL[j];
	}

	for ( int j = 0; j < inst->n_element_partition; j++)
	{
		cout << "ELEMENT\t " << j << "\t items \t"<< ITEM_ELEMENT[j] << endl;
	}

	double val_sol=compute_subset_coverage_utility(inst,CURRENT_SOL);

	cout << "val_sol\t" << val_sol << "\t item_counter \t" << item_counter << endl;

	delete []CURRENT_SOL;

	clock_t time_end=clock();
	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;


	//////////////////////////////////////////////////////////////////////////////////
	ofstream compact_file;
	compact_file.open("info_GREEDY_PARTITION.txt", ios::app);
	compact_file << fixed

			<<  val_sol << "\t"
			<<  item_counter << "\t"
			<<  solution_time << "\t"

			<<  inst->input_file<< "\t"
			<<  inst->param_file<< "\t"
			<<  inst->algorithm<< "\t"
			<<  inst->option<< "\t"
			<<  inst->timelimit<< "\t"
			<<  inst->n_items<< "\t"
			<<  inst->m_scenarios<< "\t"
			<<  inst->n_meta_items<< "\t"
			<<  inst->lambda<< "\t"
			<<  inst->delta<< "\t"
			<<  inst->seed<< "\t"
			<<  inst->TOLL_OPTIMALITY<< "\t"
			<<  inst->ProbScenario<< "\t"
			<<  inst->cardinality<< "\t"
			<<  inst->conflict_perc<< "\t"

			<<   inst->value_a	<< "\t"
			<<   inst->distribute_a << "\t"

			<< inst->KP_constraint << "\t"
			<< inst->KP_constraint_R_VALUE << "\t"
			<< inst->KP_constraint_perc_cap << "\t"
			<< inst->KP_constraint_CAPACITY << "\t"

			<< inst->USE_WORST_CASE_INSTANCE << "\t"
			<< inst->type_of_zed_function << "\t"

			<< inst->partition_constraints << "\t"
			<< inst->meta_item_per_element << "\t"
			<< inst->budget_per_element << "\t"

			<<  inst->TEST_ID << "\t"

			<<   inst->scale_factor_alpha << "\t"

			<< endl;
	compact_file.close();
	//////////////////////////////////////////////////////////////////////////////////

	//cin.get();

}



/*****************************************************************/
void greedy_algorithm_KP_CONSTRAINT(instance *inst)
/*****************************************************************/
{

	clock_t time_start=clock();

	double *CURRENT_SOL=new double[inst->n_meta_items];
	for ( int j = 0; j < inst->n_meta_items; j++){CURRENT_SOL[j]=0;}

	cout << "KP CAPACITY\t" << inst->KP_constraint_CAPACITY << endl;

	int USED_CAPACITY=0;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for ( int k = 0; k < inst->n_meta_items; k++)
	{

		cout << "iter\t" << k << "\t USED_CAPACITY \t" << USED_CAPACITY << endl;

		int best_j=-1;
		double best_val=-1;

		double val_1=compute_subset_coverage_utility(inst,CURRENT_SOL);

		for ( int j = 0; j < inst->n_meta_items; j++)
		{

			//////////////////////////////////////////////////////////
			bool OK_CONFLITS=true;
			if(inst->n_conflicts>0)
			{
				for ( int i = 0; i < inst->n_meta_items; i++)
				{
					if(i==j){continue;}

					if(CURRENT_SOL[i]>0.5 && inst->CONF_MATRIX[i][j]==1)
					{

						cout << "**CONFLICT**\t" << i << "\t" << j << endl;

						OK_CONFLITS=false;
						break;
					}
				}
			}
			//////////////////////////////////////////////////////////

			//////////////////////////////////////////////////////////
			bool OK_KP_CAP=true;

			if( (inst->KP_constraint_weights[j] + USED_CAPACITY ) > inst->KP_constraint_CAPACITY)
			{
				OK_KP_CAP=false;
			}
			//////////////////////////////////////////////////////////

			if(CURRENT_SOL[j]<0.5 && OK_CONFLITS && OK_KP_CAP)
			{

				CURRENT_SOL[j]=1;

				double val_2=compute_subset_coverage_utility(inst,CURRENT_SOL);

				CURRENT_SOL[j]=0;

				//cout << "marginal\t" << (val_2 - val_1)  << endl;

				if(best_val < (val_2 - val_1) - EPSILON_TIES )
				{
					best_val = (val_2 - val_1);
					best_j=j;

					//cout << "UPDATE\t" << best_j << endl;
				}
			}
		}

		if(best_j!=-1)
		{
			cout << "best_j\t" << best_j << "\t" << best_val << endl;
			CURRENT_SOL[best_j]=1;

			USED_CAPACITY+=inst->KP_constraint_weights[best_j];
			//cin.get();
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int item_counter=0;
	for ( int j = 0; j < inst->n_meta_items; j++)
	{
		cout << "METAITEM\t" << j << "\t->\t"<< CURRENT_SOL[j] << endl;
		item_counter+=CURRENT_SOL[j];
	}

	double val_sol=compute_subset_coverage_utility(inst,CURRENT_SOL);

	cout << "val_sol\t" << val_sol << "\t item_counter \t" << item_counter << endl;

	delete []CURRENT_SOL;

	clock_t time_end=clock();
	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;


	//////////////////////////////////////////////////////////////////////////////////
	ofstream compact_file;
	compact_file.open("info_GREEDY_KP.txt", ios::app);
	compact_file << fixed

			<<  val_sol << "\t"
			<<  item_counter << "\t"
			<<  solution_time << "\t"

			<<  inst->input_file<< "\t"
			<<  inst->param_file<< "\t"
			<<  inst->algorithm<< "\t"
			<<  inst->option<< "\t"
			<<  inst->timelimit<< "\t"
			<<  inst->n_items<< "\t"
			<<  inst->m_scenarios<< "\t"
			<<  inst->n_meta_items<< "\t"
			<<  inst->lambda<< "\t"
			<<  inst->delta<< "\t"
			<<  inst->seed<< "\t"
			<<  inst->TOLL_OPTIMALITY<< "\t"
			<<  inst->ProbScenario<< "\t"
			<<  inst->cardinality<< "\t"
			<<  inst->conflict_perc<< "\t"

			<<   inst->value_a	<< "\t"
			<<   inst->distribute_a << "\t"

			<< inst->KP_constraint << "\t"
			<< inst->KP_constraint_R_VALUE << "\t"
			<< inst->KP_constraint_perc_cap << "\t"
			<< inst->KP_constraint_CAPACITY << "\t"

			<< inst->USE_WORST_CASE_INSTANCE << "\t"
			<< inst->type_of_zed_function << "\t"

			<< inst->partition_constraints << "\t"
			<< inst->meta_item_per_element << "\t"
			<< inst->budget_per_element << "\t"

			<<  inst->TEST_ID << "\t"

			<<   inst->scale_factor_alpha << "\t"

			<< endl;
	compact_file.close();
	//////////////////////////////////////////////////////////////////////////////////

}
