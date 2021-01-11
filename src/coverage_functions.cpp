

#include "coverage_functions.h"


/***********************************************************************************/
double compute_val_diff(instance *inst, double p)
/***********************************************************************************/
{

	switch (inst->type_of_zed_function)
	{

	case 1 :
		return   ( exp ( -p/inst->lambda ) ) / inst->lambda ;
		break;

	case 2 :
		return   1;
		break;

	case 3 :
		return     ( exp ( -p/inst->lambda ) ) / inst->lambda - inst->alpha;
		break;

	default :
		cout << "wrong value for type_of_zed_function";
		exit(-1);

	}

}


/***********************************************************************************/
double compute_val_funct(instance *inst, double p)
/***********************************************************************************/
{

	switch (inst->type_of_zed_function)
	{

	case 1 :
		return  (1- exp ( -p/inst->lambda ) ) ;
		break;

	case 2 :
		return   p;
		break;

	case 3 :
		return    (1- exp ( -p/inst->lambda ) ) - inst->alpha * p;
		break;

	default :
		cout << "wrong value for type_of_zed_function";
		exit(-1);

	}

}

/*****************************************************************/
double compute_subset_coverage_utility(instance *inst,double *y)
/*****************************************************************/
{


	bool *reached_item=new bool[inst->n_items];

	for ( int j = 0; j < inst->n_items; j++)
	{
		reached_item[j]=false;
	}

	for ( int k = 0; k < inst->n_meta_items; k++)
	{

		if(y[k]< 0.5){continue;}

		for (int kkk = inst->NFS[k]; kkk < inst->NFS[k+1]; kkk++ )
		{
			reached_item[inst->AFS[kkk]]=true;
		}
	}


	//	cout << "reached_item\t:";
	//	for ( int j = 0; j < inst->n_items; j++)
	//	{
	//		cout << reached_item[j] ;
	//	}
	//	cout << endl;

	double value=0;

	for ( int i = 0; i < inst->m_scenarios; i++)
	{
		double scenario_value=0;
		for ( int j = 0; j < inst->n_items; j++)
		{
			if (reached_item[j] == false) continue;
			scenario_value += inst->a[i][j];
		}


		switch (inst->type_of_zed_function)
		{

		case 1 :
			value+=(1-exp(-scenario_value/inst->lambda));
			break;

		case 2 :
			value+=scenario_value;
			break;

		case 3 :
			value+= ( ( 1-exp(-scenario_value/inst->lambda) ) - inst->alpha * scenario_value);
			break;

		default :
			cout << "wrong value for type_of_zed_function";
			exit(-1);

		}


	}

	delete []reached_item;


	return value;
}


/*****************************************************************/
double compute_subset_coverage_utility_mod_scenario(instance *inst,double *y,int i)
/*****************************************************************/
{



	bool *reached_item=new bool[inst->n_items];

	for ( int j = 0; j < inst->n_items; j++)
	{
		reached_item[j]=false;
	}

	for ( int k = 0; k < inst->n_meta_items; k++){

		if(y[k]< 0.5){continue;}

		for (int kkk = inst->NFS[k]; kkk < inst->NFS[k+1]; kkk++ )
		{
			reached_item[inst->AFS[kkk]]=true;
		}
	}


	double value=0;

	//	for ( int i = 0; i < inst->m_scenarios; i++)
	//	{
	double scenario_value=0;
	for ( int j = 0; j < inst->n_items; j++)
	{
		if (reached_item[j] == false) continue;
		scenario_value += inst->a[i][j];
	}



	switch (inst->type_of_zed_function)
	{

	case 1 :
		value+=(1-exp(-scenario_value/inst->lambda));
		break;

	case 2 :
		value+=scenario_value;
		break;

	case 3 :
		value+= ( (1-exp(-scenario_value/inst->lambda)) - inst->alpha * scenario_value);
		break;

	default :
		cout << "wrong value for type_of_zed_function";
		exit(-1);

	}


	delete []reached_item;


	return value;
}






/*****************************************************************/
int compute_subset_coverage_number(instance *inst,double *y,bool *reached_items)
/*****************************************************************/
{

	int number=0;

	for ( int j = 0; j < inst->n_items; j++)
	{
		reached_items[j]=false;
	}

	for ( int j = 0; j < inst->n_meta_items; j++)
	{

		if(y[j]< 0.5){continue;}

		for (int  k = inst->NFS[j]; k < inst->NFS[j+1]; k++ )
		{
			if(reached_items[inst->AFS[k]]==false)
			{
				reached_items[inst->AFS[k]]=true;
				number++;
			}
		}
	}

	return number;
}


/*****************************************************************/
double compute_subset_coverage_utility_scenario(instance *inst,double *y,int i)
/*****************************************************************/
{



	for ( int j = 0; j < inst->n_items; j++)
	{
		inst->reached_item[j]=false;
	}

	for ( int k = 0; k < inst->n_meta_items; k++){

		if(y[k]< 0.5){continue;}

		for (int kkk = inst->NFS[k]; kkk < inst->NFS[k+1]; kkk++ )
		{
			inst->reached_item[inst->AFS[kkk]]=true;
		}
	}


	double value=0;


	double scenario_value=0;
	for ( int j = 0; j < inst->n_items; j++)
	{
		if (inst->reached_item[j] == false) continue;
		scenario_value += inst->a[i][j];
	}


	return scenario_value;

}

/*****************************************************************/
double compute_subset_coverage_utility_scenario_fract(instance *inst,double *y,int i)
/*****************************************************************/
{



	for ( int j = 0; j < inst->n_items; j++)
	{
		inst->reached_item[j]=0.0;
	}

	for ( int k = 0; k < inst->n_meta_items; k++)
	{

		for (int kkk = inst->NFS[k]; kkk < inst->NFS[k+1]; kkk++ )
		{
			inst->reached_item[inst->AFS[kkk]]+=y[k];
		}
	}


	for ( int k = 0; k < inst->n_items; k++)
	{
		if(inst->reached_item[k]>= 1 - 0.00001)
		{
			inst->reached_item[k]=1.0;
		}

		if(inst->reached_item[k]<=  0.00001)
		{
			inst->reached_item[k]=0.0;
		}
	}


	double scenario_value=0;

	for ( int j = 0; j < inst->n_items; j++)
	{
		scenario_value += (inst->reached_item[j] * inst->a[i][j]);
	}

	return scenario_value;

}

