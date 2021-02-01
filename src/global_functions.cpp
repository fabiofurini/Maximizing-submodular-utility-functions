

#include "global_functions.h"


/*****************************************************************/
void set_lambda_and_alpha(instance *inst)
/*****************************************************************/
{
	double D_tot=0.0;

	for ( int i = 0; i < inst->m_scenarios; i++)
	{
		double D=0.0;

		for(int j=0;j<inst->n_items;j++)
		{
			D+=inst->a[i][j];
		}

		cout << "scenario\t" << i << "\t total demand\t" << D << "\t total \t" << D_tot << endl;

		D_tot+=D;
	}

	inst->AVERAGE_DEMAND=(double) D_tot / (double)inst->m_scenarios;

	cout << fixed << "average demand\t" << (double) D_tot / (double)inst->m_scenarios << endl;

	inst->lambda_orig=inst->lambda;

	inst->lambda= inst->lambda * inst->AVERAGE_DEMAND;

	cout << "SCALED value of LAMBDA\t" << inst->lambda << endl;

	cout << "scale_factor_alpha\t" << inst->scale_factor_alpha << endl;

	//	inst->alpha = 1.0 / (inst->scale_factor_alpha * inst->AVERAGE_DEMAND);

	inst->alpha = inst->scale_factor_alpha / (inst->lambda);

	cout << "SCALED value of alpha\t" << inst->alpha << endl;
}

/*****************************************************************/
void scale_obj_coefficient(instance *inst)
/*****************************************************************/
{
	double D_tot=0.0;

	for ( int i = 0; i < inst->m_scenarios; i++)
	{
		double D=0.0;

		for(int j=0;j<inst->n_items;j++)
		{
			D+=inst->a[i][j];
		}

		cout << "scenario\t" << i << "\t total demand\t" << D << "\t total \t" << D_tot << endl;

		D_tot+=D;
	}

	inst->AVERAGE_DEMAND=(double) D_tot / (double)inst->m_scenarios;

	cout << fixed << "average demand\t" << (double) D_tot / (double)inst->m_scenarios << endl;

	for ( int i = 0; i < inst->m_scenarios; i++)
	{

		for(int j=0;j<inst->n_items;j++)
		{
			inst->a[i][j]=inst->a[i][j]/inst->AVERAGE_DEMAND;
		}

	}


}


/*****************************************************************/
void set_alpha(instance *inst)
/*****************************************************************/
{

	inst->alpha = inst->scale_factor_alpha/(inst->lambda);
}

/*****************************************************************/
void compute_curvature(instance *inst)
/*****************************************************************/
{

	clock_t time_start=clock();

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "***COMPUTE CURVATURE***\n";

	inst->curvature=CPX_INFBOUND;

	double *YY=new double[inst->n_meta_items];

	for ( int j = 0; j < inst->n_meta_items; j++)
	{
		YY[j]=1;
	}

	double val_all=compute_subset_coverage_utility(inst,YY);

	for ( int j = 0; j < inst->n_meta_items; j++)
	{
		YY[j]=0;
		double rho_N_minus_j= val_all - compute_subset_coverage_utility(inst,YY);
		YY[j]=1;


		double rho_singleton=0;

		for ( int k = 0; k < inst->m_scenarios; k++ )
		{

			double dummy=0;
			for (int i = inst->NFS[j]; i < inst->NFS[j+1]; i++ )
			{
				dummy+=inst->a[k][inst->AFS[i]];
			}


			switch (inst->type_of_zed_function)
			{

			case 1 :
				rho_singleton+=(1-exp(-dummy/inst->lambda));
				break;

			case 2 :
				rho_singleton+=dummy;
				break;

			case 3 :
				rho_singleton+= ( (1-exp(-dummy/inst->lambda))- inst->alpha * dummy );
				break;

			default :
				cout << "wrong value for type_of_zed_function";
				exit(-1);

			}



		}

		//cout << "meta item\t" << j << "\t rho_N_minus_j \t" << rho_N_minus_j << "\t rho_singleton \t" << rho_singleton << "\t ratio \t" << rho_N_minus_j/ rho_singleton << endl;

		if( inst->curvature > (rho_N_minus_j / rho_singleton) )
		{
			inst->curvature = (rho_N_minus_j / rho_singleton);
		}

	}

	inst->curvature=1-inst->curvature;

	cout << "curvature\t" << inst->curvature << endl;

	delete[]YY;

	inst->curvature_scenario_avg=0;
	inst->curvature_scenario_min=CPX_INFBOUND;
	inst->curvature_scenario_max=0;


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for ( int k = 0; k < inst->m_scenarios; k++)
	{
		double curvature_scenario=CPX_INFBOUND;

		double rho_N_minus_j= 0;
		double rho_singleton = 0;

		for ( int j = 0; j < inst->n_meta_items; j++)
		{

			double dummy=0;
			for (int i = inst->NFS[j]; i < inst->NFS[j+1]; i++ )
			{
				if(inst->DM[inst->AFS[i]]==1)
				{
					dummy+=inst->a[k][inst->AFS[i]];
				}
			}

			rho_N_minus_j=dummy;

			dummy=0;
			for (int i = inst->NFS[j]; i < inst->NFS[j+1]; i++ )
			{
				dummy+=inst->a[k][inst->AFS[i]];
			}

			rho_singleton=dummy;

			if( curvature_scenario > (rho_N_minus_j / rho_singleton) )
			{
				curvature_scenario = (rho_N_minus_j / rho_singleton);
			}

		}

		curvature_scenario=1-curvature_scenario;

		inst->curvature_scenario_avg+=curvature_scenario;

		if(inst->curvature_scenario_min>curvature_scenario)
		{
			inst->curvature_scenario_min=curvature_scenario;
		}

		if(inst->curvature_scenario_max<curvature_scenario)
		{
			inst->curvature_scenario_max=curvature_scenario;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	inst->curvature_scenario_avg=inst->curvature_scenario_avg/inst->m_scenarios;

	cout << "curvature_scenario_avg\t" << inst->curvature_scenario_avg << endl;
	cout << "curvature_scenario_max\t" << inst->curvature_scenario_max << endl;
	cout << "curvature_scenario_min\t" << inst->curvature_scenario_min << endl;



	clock_t time_end=clock();
	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;

	//////////////////////////////////////////////////////////////////////////////////
	ofstream compact_file;
	compact_file.open("info_CURVATURE.txt", ios::app);
	compact_file << fixed

			<<  inst->curvature << "\t"
			<<  inst->curvature_scenario_avg << "\t"
			<<  inst->curvature_scenario_max << "\t"
			<<  inst->curvature_scenario_min << "\t"

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


	exit(-1);
}

/*****************************************************************/
void read_param_file(instance *inst)
/*****************************************************************/
{

	ifstream in(inst->param_file);
	if ( !in )
	{
		cout << "File with PARAMATERS could not be opened. " << endl;
		exit(1);
	}

	char DUMMY[1000];


	in >> DUMMY; in >> inst->number_of_CPU;

	in >> DUMMY; in >> inst->TOLL_DERIVATIVE;

	in >> DUMMY; in >> inst->TOLL_VIOL;

	in >> DUMMY; in >> inst->TOLL_VIOL_FRAC_BEN;

	in >> DUMMY; in >> inst->FLAG_GREEDY_SOL;



	in.close();

}

// return a integer random value in range min-max
/*****************************************************************/
int randomBETWEEN(int min,int max)
/*****************************************************************/
{
	return rand() % (max - min +1) +min;
}

// return random value in range min-max
/*****************************************************************/
double randomBETWEEN_double(int min,int max)
/*****************************************************************/
{
	return (rand()/(double)RAND_MAX)*(max-min) + min;
}


// return a random value in range 0.0-1.0
/*****************************************************************/
double random01()
/*****************************************************************/
{
	return ((double) rand() / RAND_MAX);
}


/*****************************************************************/
void from_int_to_set(int *set,int number,int dim)
/*****************************************************************/
{

	for(int i=0;i<dim;i++)
	{
		set[i]=0;
	}

	int bit;
	double num1=number;
	double num2=2;

	if (number!=0)
	{

		bit=(log(num1)/log(num2));
		bit++;

		for(int i=0;i<bit;i++)
		{
			set[i]=(number%2);
			number=number/2;
		}
	}
}


/*****************************************************************/
void init_data(instance *inst)
/*****************************************************************/
{


	inst->KP_constraint_CAPACITY=0;

	inst->KP_constraint_weights=new double[inst->n_meta_items];

	inst->reached_item=new double[inst->n_items];

	inst->meta_item_element_partition=new int[inst->n_meta_items];

	inst->GREEDY_SOL=new double[inst->n_meta_items];

}

/*****************************************************************/
void free_data(instance *inst)
/*****************************************************************/
{


	cout << "FREEEEE\n\n";

	delete []inst->GREEDY_SOL;

	delete []inst->KP_constraint_weights;


	for ( int j = 0; j < inst->n_meta_items; j++ )
	{
		delete []inst->CONF_MATRIX[j];
	}
	delete []inst->CONF_MATRIX;

	delete []inst->reached_item;

	delete []inst->NFS;
	delete []inst->AFS;
	delete []inst->DP;
	delete []inst->NBS;
	delete []inst->ABS;
	delete []inst->DM;

	if(inst->FLAG_INSTANCE_MP==true)
	{
		delete []inst->weight_item_MP;
	}

	for ( int i = 0; i < inst->m_scenarios; i++)
	{
		delete []inst->a[i];
	}
	delete []inst->a;

	delete []inst->meta_item_element_partition;

	if(inst->FLAG_INSTANCE_LOCATION==true)
	{

		delete []inst->item_OK;

		free(inst->x_location);
		free(inst->y_location);

		free(inst->x_client);
		free(inst->y_client);

	}

}










/*****************************************************************/
void draw_grid_sol_market(instance *inst,int *yy,int *zz)
/*****************************************************************/
{

	double scale_factor=2.5;

	cout << "\n\nWRITING SOLUTION FIGURE\n\n";
	//cin.get();

	ofstream compact_file;
	compact_file.open("./figures/gridSol_market.tex");
	compact_file
	<< "\\documentclass[]{article}" << "\n"
	<< "\\usepackage{tikz}" << "\n"
	<< "\\begin{document}" << "\n"
	<< "\\begin{tikzpicture}" << "\n"

	<< endl;

	int **coverage_clients=new int*[inst->n_items];
	for ( int i = 0; i < inst->n_items; i++ )
	{
		coverage_clients[i]=new int[inst->m_scenarios];

		for ( int k = 0; k < inst->m_scenarios; k++ )
		{

			coverage_clients[i][k]=0;

			for (int  j = inst->NBS[i]; j < inst->NBS[i+1]; j++ )
			{
				if(yy[inst->ABS[j]]==1)
				{
					coverage_clients[i][k]=coverage_clients[i][k]+1;
				}
			}

		}
	}


	int **color=new int*[200];
	for ( int k = 0; k < inst->m_scenarios; k++ )
	{
		color[k]=new int[3];

		color[k][0]=1;
		color[k][1]=1;
		color[k][2]=1;

	}

	color[0][0]=1;
	color[0][1]=1;
	color[0][2]=1;

	color[1][0]=255;
	color[1][1]=255;
	color[1][2]=255;

	color[2][0]=255;
	color[2][1]=0;
	color[2][2]=0;

	color[3][0]=0;
	color[3][1]=255;
	color[3][2]=0;

	color[5][0]=0;
	color[5][1]=0;
	color[5][2]=255;

	color[6][0]=255;
	color[6][1]=255;
	color[6][2]=0;

	color[7][0]=0;
	color[7][1]=255;
	color[7][2]=255;

	color[8][0]=255;
	color[8][1]=0;
	color[8][2]=255;

	color[9][0]=255;
	color[9][1]=215;
	color[9][2]=0;

	double **perturb=new double*[200];
	for ( int k = 0; k < inst->m_scenarios; k++ )
	{

		double angolo=(2*3.1428/inst->m_scenarios)*k;

		perturb[k]=new double[2];

		perturb[k][0]=0.5*cos(angolo);
		perturb[k][1]=0.5*sin(angolo);

	}


	for ( int i = 0; i < inst->n_items; i++ )
	{
		for ( int k = 0; k < inst->m_scenarios; k++ )
		{

			//			if(zz[i]==1){
			//				compact_file <<
			//						"\\draw[fill=black] (" <<
			//						(inst->x_client[i] )/scale_factor  <<
			//						"," <<
			//						(inst->y_client[i])/scale_factor <<
			//						") circle (0.1);\n"
			//						<< endl;
			//			}

			if(zz[i]==1 && inst->a[k][i]==1)
			{
				if(coverage_clients[i][k]==1)
				{
					compact_file <<
							"\\draw[line width=0.01mm,fill={rgb:red," << color[k][0] << " ;green," << color[k][1] << ";blue," <<  color[k][2] << "}] (" <<
							(inst->x_client[i] + perturb[k][0])/scale_factor  <<
							"," <<
							(inst->y_client[i]+ perturb[k][1] )/scale_factor <<
							") circle (0.06);\n"
							<< endl;
				}
				else
				{
					compact_file <<
							"\\draw[line width=0.04mm,fill={rgb:red," << color[k][0] << " ;green," << color[k][1] << ";blue," <<  color[k][2] << "}] (" <<
							(inst->x_client[i]+ perturb[k][0])/scale_factor  <<
							"," <<
							(inst->y_client[i]+ perturb[k][1])/scale_factor <<
							") circle (0.06);\n"
							<< endl;
				}
			}
		}
	}

	for ( int i = 0; i < inst->n_items; i++ )
	{
		for ( int k = 0; k < inst->m_scenarios; k++ )
		{

			if(zz[i]==1){
				compact_file <<
						"\\draw[fill=black] (" <<
						(inst->x_client[i] )/scale_factor  <<
						"," <<
						(inst->y_client[i])/scale_factor <<
						") circle (0.1);\n"
						<< endl;
			}
		}
	}



	for ( int j = 0; j < inst->n_meta_items; j++ )
	{
		if(yy[j]==1){
			compact_file <<
					"\\draw[line width=0.7mm,circle,fill=white] (" <<
					inst->x_location[j]/scale_factor  <<
					"," <<
					inst->y_location[j]/scale_factor <<
					") circle (0.2);"
					<< endl;
		}
	}

	for ( int j = 0; j < inst->n_meta_items; j++ )
	{
		if(yy[j]==1){
			compact_file <<
					"\\node[circle, draw, ultra thick, dotted,line width=2.6pt,fill=none,minimum size=" <<
					inst->delta <<
					"cm]at(" <<
					inst->x_location[j]/scale_factor  <<
					"," <<
					inst->y_location[j]/scale_factor <<
					"){};"
					<< endl;
		}
	}

	int size_grid=30;

	double point=size_grid/scale_factor;

	compact_file << "\\draw[ line width=0.7mm ] (0,0) -- (0,"<< point << ");" << endl;
	compact_file << "\\draw[ line width=0.7mm ] (0,0) -- (" << point << ",0);" << endl;
	compact_file << "\\draw[ line width=0.7mm ] ("<< point << ","<<  point <<") -- (" << point << ",0);" << endl;
	compact_file << "\\draw[ line width=0.7mm ] ("<< point << ","<<  point <<") -- (0," << point <<  ");" << endl;
	compact_file
	<< "\\end{tikzpicture}" << "\n"
	<< "\\end{document}" << "\n"
	<< "\\begin{document}" << "\n"
	<< endl;

	//delete []coverage_clients;


	compact_file.close();

}

///***************************************************************************/
void write_data_dat_file(instance *inst)
/***************************************************************************/
{

	double val=1;

	cout << "**********ATTENTION FIX THE CORRECT VALUE OF VAL***********\n";
	exit(-1);


	char DUMMY[1000];

	ofstream compact_names;

	if(inst->FLAG_INSTANCE_LOCATION==false){
		//sprintf(DUMMY,"../ampl_code/data.dat");
		sprintf(DUMMY,"ForBARON_NEW_FINANCE/FINANCEu%d_v%d_k%d_l%.2f_b%.2f_r%.2f_s%d_data.dat",
				inst->n_items,inst->n_meta_items,inst->m_scenarios,inst->lambda,val,inst->delta,inst->seed);

		printf("ForBARON_NEW_FINANCE/FINANCEu%d_v%d_k%d_l%.2f_b%.2f_r%.2f_s%d_data.dat.log\n",
				inst->n_items,inst->n_meta_items,inst->m_scenarios,inst->lambda,val,inst->delta,inst->seed);

		compact_names.open("ForBARON_NEW_FINANCE/names.txt", ios::app);

	}else{

		//sprintf(DUMMY,"../ampl_code/data.dat");
		sprintf(DUMMY,"ForBARON_NEW_LOCATION/LOCATIONu%d_v%d_k%d_l%.2f_b%.2f_r%.2f_s%d_data.dat",
				inst->n_items,inst->n_meta_items,inst->m_scenarios,inst->lambda,val,inst->delta,inst->seed);

		printf("ForBARON_NEW_LOCATION/LOCATIONu%d_v%d_k%d_l%.2f_b%.2f_r%.2f_s%d_data.dat.log\n",
				inst->n_items,inst->n_meta_items,inst->m_scenarios,inst->lambda,val,inst->delta,inst->seed);

		compact_names.open("ForBARON_NEW_LOCATION/names.txt", ios::app);

	}


	compact_names <<  DUMMY << endl;
	compact_names.close();

	ofstream compact_file;
	compact_file.open(DUMMY);


	compact_file << "param vv := " << inst->n_items << ";"<<   endl;
	compact_file << "param uu := " << inst->n_meta_items<< ";"<< endl;
	compact_file << "param kk := " << inst->m_scenarios << ";" << endl;

	compact_file << "param lambda := " << inst->lambda << ";" << endl;


	compact_file << "param b := " << val << ";" << endl;

	cout << "************TO BE FIXED*************";

	//	compact_file << "param a := " << endl;
	//	for ( int i = 0; i < inst->n_meta_items; i++ )
	//	{
	//		compact_file << i+1 << " "<< inst->b[i]  << endl;
	//	}
	//	compact_file << ";\n";

	cout << "************TO BE FIXED*************";


	compact_file << "set UV := " << endl;
	for ( int j = 0; j < inst->n_meta_items; j++){

		for (int  k = inst->NFS[j]; k < inst->NFS[j+1]; k++ )
		{
			compact_file << "(" << j+1 << "," << inst->AFS[k]+1 << ")" << endl;
		}
	}
	compact_file << ";\n";

	compact_file << "param c (tr):= " << endl;
	for ( int j = 0; j< inst->n_items; j++ )
	{
		compact_file << j+1 << " ";
	}
	compact_file << ":=\n";

	for(int k=0;k<inst->m_scenarios;k++)
	{
		compact_file << k+1 << "  ";
		for(int j=0;j<inst->n_items;j++)
		{
			compact_file << inst->a[k][j] << " ";
		}
		compact_file <<  endl;
	}

	compact_file << ";\n";


	compact_file.close();


	cout << "ESCO FURTIVO DOPO AVER SCRITTO L'istanza....\n\n";
	exit(-1);



}


