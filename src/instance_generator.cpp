

#include "instance_generator.h"

//#define print_stars
//#define print_a_coefficients

//#define print_weights



//#define print_neighbourhoods

/***************************************************************************/
void generate_WORT_CASE_INSTANCES(instance *inst, int k)
/***************************************************************************/
{

	cout << "INPUT PARAMETER k->\t" << k << endl;

	inst->n_items=k*k+k;
	inst->n_meta_items=2*k;
	inst->m_scenarios=1;

	cout << "n_items\t" << inst->n_items << endl;
	cout << "n_meta_items\t" << inst->n_meta_items << endl;
	cout << "m_scenarios\t" << inst->m_scenarios << endl;

	inst->a=new double*[inst->m_scenarios];
	for ( int i = 0; i < inst->m_scenarios; i++)
	{
		inst->a[i]=new double[inst->n_items];
	}


	int *father=new int[inst->n_items];
	int *mather=new int[inst->n_items];


	int counter_item=0;
	for ( int i = 0; i <= k; i++)
	{

		//////////////////////////////////////////////////////////
		if(i==0)
		{
			for ( int j = 1; j <= k; j++)
			{


				if(j==1)
				{
					cout << k-1 << "\t";
					inst->a[0][counter_item]= k - 1;

					mather[counter_item]= j + k -1;
					father[counter_item++]= -1;

				}
				else
				{
					cout << k-2 << "\t";
					inst->a[0][counter_item]= k - 2;

					mather[counter_item]= j + k -1;
					father[counter_item++]= -1;

				}
			}
			cout << endl;
		}
		//////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////
		if(i==1)
		{
			for ( int j = 1; j <= k; j++)
			{


				if(j==1)
				{
					cout << 0 << "\t";
					inst->a[0][counter_item]=0;

					mather[counter_item]= j + k -1;
					father[counter_item++]= i -1;


				}
				else
				{
					cout << 1 << "\t";
					inst->a[0][counter_item]= 1;

					mather[counter_item]= j + k -1;
					father[counter_item++]= i -1;

				}
			}
			cout << endl;
		}
		//////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////
		if(i>1)
		{
			for ( int j = 1; j <= k; j++)
			{

				cout << pow  ( (double) k / (double) (k-1) , (double) i-2  ) << "\t";

				inst->a[0][counter_item]= pow  ( (double) k / (double) (k-1) , (double) i-2  );

				mather[counter_item]= j + k - 1;
				father[counter_item++]= i - 1;

			}
			cout << endl;
		}
		//////////////////////////////////////////////////////////

	}


	////////////////////////////////////////////////////////////////////////////////////////////
	cout << "BUILDING neighbourhoods\n";
	vector < vector < int > > dummy;
	for ( int i = 0; i < inst->n_items; i++ )
	{

		vector < int > local_dummy;

		if(father[i]!=-1)
		{
			local_dummy.push_back(father[i]);
		}

		local_dummy.push_back(mather[i]);

		dummy.push_back(local_dummy);
	}
	////////////////////////////////////////////////////////////////////////////////////////////



	delete []father;
	delete []mather;



#ifdef print_neighbourhoods
	for ( int i = 0; i < inst->n_items; i++ ){

		cout << "ITEM\t" << i << "\t number of META ITEMS \t " << dummy[i].size() << "\t\t";

		for ( int j = 0; j <  dummy[i].size(); j++)
		{
			cout << dummy[i][j] << " ";
		}

		cout << endl;
	}
	cout << endl;
#endif




	//////////////////////////////////////////////
	for ( int j = 0; j < inst->n_items; j++ )
	{
		if((int)(dummy[j].size())==0)
		{
			cout << "ITEM NOT COVERED BY ANY META-ITEM....\n\n";
			exit(-1);
		}
	}
	//////////////////////////////////////////////

	inst->item_not_covered=0;
	inst->item_single_covered=0;
	for ( int j = 0; j < inst->n_items; j++ )
	{
		if((int)(dummy[j].size())==0)
		{
			inst->item_not_covered++;
		}
		if((int)(dummy[j].size())==1)
		{
			inst->item_single_covered++;
		}
	}
	cout << "item_not_covered\t" << inst->item_not_covered << endl;
	cout << "item_single_covered\t" << inst->item_single_covered << endl;

	inst->NFS=new int[inst->n_meta_items+1];
	inst->AFS=new int[MAX_CONNECTIONS];
	inst->DP=new int[inst->n_meta_items];
	inst->NBS=new int[inst->n_items+1];
	inst->ABS=new int[MAX_CONNECTIONS];
	inst->DM=new int[inst->n_items];

	inst->counter_c=0;


	cout << "\n***BUILDING BACKWARD STARS\t" << endl;

	inst->NBS[0]=0;
	for ( int j = 0; j < inst->n_items; j++ ){

		inst->DM[j]=dummy[j].size();

		inst->NBS[j+1]=inst->NBS[j]+dummy[j].size();

		for ( int jjj = 0; jjj <  dummy[j].size(); jjj++){

			inst->ABS[inst->counter_c++]=dummy[j][jjj];
		}
	}

	cout << "\n***REAL_CONNECTIONS\t" << inst->counter_c << endl;


	cout << "\n***BUILDING FORWARD STARS\t" << endl;
	inst->counter_l=0;

	inst->NFS[0]=0;
	for ( int k = 0; k < inst->n_meta_items; k++)
	{

		int local_counter=0;

		for ( int j = 0; j < inst->n_items; j++ )
		{
			bool found=false;
			for ( int jj = 0; jj <  dummy[j].size() && !found; jj++)
			{
				if(k==dummy[j][jj])
				{
					inst->AFS[inst->counter_l++]=j;
					local_counter++;
					found=true;
				}

			}
		}

		inst->DP[k]=local_counter;

		inst->NFS[k+1]=inst->NFS[k]+local_counter;
	}

	cout << "\n***REAL_CONNECTIONS\t" << inst->counter_l << endl;


#ifdef print_stars

	for(int j=0;j<inst->n_items;j++)
	{
		cout << "Backward star of item \t" << j << "\t size of his neighbourhood\t" << inst->DM[j] << endl;
		for (int  k = inst->NBS[j]; k < inst->NBS[j+1]; k++ )
		{
			cout << "meta item\t" << inst->ABS[k]  << endl;
		}

	}
	cin.get();

	for ( int k = 0; k < inst->n_meta_items; k++){
		cout << "Forward star of meta item \t" << k << "\t size of his neighbourhood\t" << inst->DP[k] << endl;
		for (int  kk = inst->NFS[k]; kk < inst->NFS[k+1]; kk++ )
		{
			cout << "item\t" << inst->AFS[kk] << endl;
		}
	}
	cin.get();

#endif




}

/***************************************************************************/
void load_partition(instance *inst)
/***************************************************************************/
{

	if(inst->budget_per_element>=inst->meta_item_per_element)
	{
		cout << "budget too large!!";
		exit(-1);
	}


	if(inst->meta_item_per_element <=0 || inst->meta_item_per_element >= inst->n_meta_items)
	{
		cout << "meta_item_per_element -> not valid";
		exit(-1);
	}

	if(inst->n_meta_items%inst->meta_item_per_element!=0)
	{
		cout << "\n***ATTENTION THE LAST ELEMENT HAS LESS META ITEMS****\n";

		cout << "The number of meta-items is not a multiple of meta_item_per_element....\n";

		exit(-1);
	}

	inst->n_element_partition=ceil((double)inst->n_meta_items/(double)inst->meta_item_per_element);

	cout << "n_element_partition\t" << inst->n_element_partition << endl;

	int counter_meta=0;
	for ( int p = 0; p < inst->n_element_partition; p++ )
	{

		for ( int kk = 0; kk < inst->meta_item_per_element; kk++ )
		{
			inst->meta_item_element_partition[counter_meta++]=p;

			if(counter_meta==inst->n_meta_items)
			{
				break;
			}
		}
	}

	for ( int kk = 0; kk < inst->n_meta_items; kk++ )
	{
		cout << "meta_item \t" << kk << "\tin element ->\t" << inst->meta_item_element_partition[kk] << endl;
	}

	//exit(-1);
}

/***************************************************************************/
void load_weights(instance *inst)
/***************************************************************************/
{


	////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "compute meta-item profits\n";


	cout << "KP_constraint_R_VALUE\t" << inst->KP_constraint_R_VALUE << endl;

	double *pp=new double [inst->n_meta_items];

	for ( int j = 0; j < inst->n_meta_items; j++)
	{
		pp[j]=0.0;

		for ( int k = 0; k < inst->m_scenarios; k++ )
		{

			double dummy=0;
			for (int i = inst->NFS[j]; i < inst->NFS[j+1]; i++ )
			{
				dummy+=inst->a[k][inst->AFS[i]];
			}

			pp[j]+=dummy;

		}

		pp[j]/=inst->m_scenarios;

	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////


	double sum_weight=0;



	switch (inst->KP_constraint)
	{

	case 1 :
		cout << "Uncorrelated\n";

		for(int i=0; i<inst->n_meta_items; i++)
		{

			inst->KP_constraint_weights[i]=randomBETWEEN(1,inst->KP_constraint_R_VALUE);

			sum_weight+=inst->KP_constraint_weights[i];
		}

		break;

	case 2 :
		cout << "Weakly correlated\n";

		for(int i=0; i<inst->n_meta_items; i++)
		{

			int LB_P= max( 1.0 , (int)(pp[i]) - (inst->KP_constraint_R_VALUE / 10.0) );

			int UB_P= (int)(pp[i]) + (inst->KP_constraint_R_VALUE / 10.0);

			inst->KP_constraint_weights[i]=randomBETWEEN(LB_P,UB_P);

			sum_weight+=inst->KP_constraint_weights[i];
		}


		break;

	case 3 :
		cout << "Strongly correlated\n";

		for(int i=0; i<inst->n_meta_items; i++)
		{

			inst->KP_constraint_weights[i] = (int)(pp[i]) + (inst->KP_constraint_R_VALUE / 10);

			sum_weight+=inst->KP_constraint_weights[i];
		}


		break;

	case 4 :
		cout << "Almost strongly correlated\n";

		for(int i=0; i<inst->n_meta_items; i++)
		{

			int LB_P= (int)(pp[i]) + (inst->KP_constraint_R_VALUE / 10.0) - (inst->KP_constraint_R_VALUE / 500.0);

			int UB_P= (int)(pp[i]) + (inst->KP_constraint_R_VALUE / 10.0) + (inst->KP_constraint_R_VALUE / 500.0);

			inst->KP_constraint_weights[i] =randomBETWEEN(LB_P,UB_P);

			sum_weight+=inst->KP_constraint_weights[i];
		}


		break;

	case 5 :
		cout << "Subset-sum\n";


		for(int i=0; i<inst->n_meta_items; i++)
		{

			inst->KP_constraint_weights[i] = pp[i];

			sum_weight+=inst->KP_constraint_weights[i];
		}


		break;

	default :
		cout << "wrong value for KP_constraint";
		exit(-1);

	}

	inst->KP_constraint_CAPACITY=  (int)(sum_weight * inst->KP_constraint_perc_cap/100.0);

	cout << "KP_constraint_CAPACITY\t" <<  inst->KP_constraint_CAPACITY << endl;

#ifdef 	print_weights

	for(int i=0; i<inst->n_meta_items; i++)
	{
		cout << "meta item\t" << pp[i] << "\t" << inst->KP_constraint_weights[i] << endl;
	}

	exit(-1);
#endif

	delete[]pp;
}




/***************************************************************************/
int load_conflicts(instance *inst)
/***************************************************************************/
{

	int n_conflicts=0;

	inst->CONF_MATRIX=new int*[inst->n_meta_items];
	for ( int j = 0; j < inst->n_meta_items; j++ )
	{
		inst->CONF_MATRIX[j]=new int[inst->n_meta_items];
	}

	for ( int j = 0; j < inst->n_meta_items; j++ )
	{
		for ( int i = 0; i < inst->n_meta_items; i++ )
		{
			inst->CONF_MATRIX[j][i]=0;
		}
	}

	for(int i=0; i<inst->n_meta_items; i++)
	{

		for(int j=i+1; j<inst->n_meta_items; j++)
		{

			double n_random=randomBETWEEN_double(0,1);

			if(n_random<=inst->conflict_perc)
			{

				//					cout << "n_random\t" << n_random << endl;
				//					cout << "conflict\t" << i << "\t" << j << endl;

				inst->CONF_MATRIX[i][j]=1;
				inst->CONF_MATRIX[j][i]=1;

				n_conflicts++;
			}
		}
	}

	return n_conflicts;
}


/***************************************************************************/
double distance(double x1,double y1,double x2, double y2)
/***************************************************************************/
{


	double distancex = (x2 - x1)*(x2 - x1);
	double distancey = (y2 - y1)*(y2 - y1);

	return   sqrt(distancex + distancey);
}

/***************************************************************************/
void generate_data_LOCATION(instance *inst)
/***************************************************************************/
{


	double size_grid=30;

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	inst->x_location= (double *) calloc(inst->n_meta_items, sizeof(double));
	inst->y_location= (double *) calloc(inst->n_meta_items, sizeof(double));

	inst->x_client= (double *) calloc(inst->n_items, sizeof(double));
	inst->y_client= (double *) calloc(inst->n_items, sizeof(double));


	cout << "SAMPLING POINT IN THE GRID\t" << 0 << "\t" << size_grid << endl;

	for ( int j = 0; j < inst->n_meta_items; j++ )
	{


		double x=randomBETWEEN_double(0,size_grid);
		double y=randomBETWEEN_double(0,size_grid);
		inst->x_location[j]=x;
		inst->y_location[j]=y;

	}

	inst->item_OK=new bool[inst->n_items];

	for ( int j = 0; j < inst->n_items; j++ )
	{
		inst->item_OK[j]=false;
	}


	for ( int i = 0; i < inst->n_items; i++ )
	{

		double x;
		double y;

		//sampling the location of the item until it is covered by at least one meta-item
		while(inst->item_OK[i]==false)
		{
			x=randomBETWEEN_double(0,size_grid);
			y=randomBETWEEN_double(0,size_grid);

			for ( int j = 0; j < inst->n_meta_items; j++)
			{

				if(distance(x,y,inst->x_location[j],inst->y_location[j]) < inst->delta)
				{
					inst->item_OK[i]=true;
					break;
				}
			}
		}

		inst->x_client[i]=x;
		inst->y_client[i]=y;
	}


	cout << "BUILDING neighbourhoods\n";
	vector < vector < int > > dummy;
	for ( int i = 0; i < inst->n_items; i++ )
	{

		vector < int > local_dummy;

		for ( int j = 0; j < inst->n_meta_items; j++){

			if(distance(inst->x_client[i],inst->y_client[i],inst->x_location[j],inst->y_location[j]) < inst->delta)
			{
				local_dummy.push_back(j);
			}
		}
		dummy.push_back(local_dummy);
	}


	inst->a=new double*[inst->m_scenarios];
	for ( int i = 0; i < inst->m_scenarios; i++)
	{
		inst->a[i]=new double[inst->n_items];
	}


#ifdef print_neighbourhoods
	for ( int i = 0; i < inst->n_items; i++ ){

		cout << "ITEM\t" << i << "\t number of META ITEMS \t " << dummy[i].size() << "\t\t";

		for ( int j = 0; j <  dummy[i].size(); j++)
		{
			cout << dummy[i][j] << " ";
		}

		cout << endl;
	}
	cout << endl;
#endif




	//////////////////////////////////////////////
	for ( int j = 0; j < inst->n_items; j++ )
	{
		if((int)(dummy[j].size())==0)
		{
			cout << "ITEM NOT COVERED BY ANY META-ITEM....\n\n";
			exit(-1);
		}
	}
	//////////////////////////////////////////////

	inst->item_not_covered=0;
	inst->item_single_covered=0;
	for ( int j = 0; j < inst->n_items; j++ )
	{
		if((int)(dummy[j].size())==0)
		{
			inst->item_not_covered++;
		}
		if((int)(dummy[j].size())==1)
		{
			inst->item_single_covered++;
		}
	}
	cout << "item_not_covered\t" << inst->item_not_covered << endl;
	cout << "item_single_covered\t" << inst->item_single_covered << endl;

	inst->NFS=new int[inst->n_meta_items+1];
	inst->AFS=new int[MAX_CONNECTIONS];
	inst->DP=new int[inst->n_meta_items];
	inst->NBS=new int[inst->n_items+1];
	inst->ABS=new int[MAX_CONNECTIONS];
	inst->DM=new int[inst->n_items];

	inst->counter_c=0;


	cout << "\n***BUILDING BACKWARD STARS\t" << endl;

	inst->NBS[0]=0;
	for ( int j = 0; j < inst->n_items; j++ ){

		inst->DM[j]=dummy[j].size();

		inst->NBS[j+1]=inst->NBS[j]+dummy[j].size();

		for ( int jjj = 0; jjj <  dummy[j].size(); jjj++){

			inst->ABS[inst->counter_c++]=dummy[j][jjj];
		}
	}

	cout << "\n***REAL_CONNECTIONS\t" << inst->counter_c << endl;


	cout << "\n***BUILDING FORWARD STARS\t" << endl;
	inst->counter_l=0;

	inst->NFS[0]=0;
	for ( int k = 0; k < inst->n_meta_items; k++)
	{

		int local_counter=0;

		for ( int j = 0; j < inst->n_items; j++ )
		{
			bool found=false;
			for ( int jj = 0; jj <  dummy[j].size() && !found; jj++)
			{
				if(k==dummy[j][jj])
				{
					inst->AFS[inst->counter_l++]=j;
					local_counter++;
					found=true;
				}

			}
		}

		inst->DP[k]=local_counter;

		inst->NFS[k+1]=inst->NFS[k]+local_counter;
	}

	cout << "\n***REAL_CONNECTIONS\t" << inst->counter_l << endl;


#ifdef print_stars

	for(int j=0;j<inst->n_items;j++)
	{
		cout << "Backward star of item \t" << j << "\t size of his neighbourhood\t" << inst->DM[j] << endl;
		for (int  k = inst->NBS[j]; k < inst->NBS[j+1]; k++ )
		{
			cout << "meta item\t" << inst->ABS[k]  << endl;
		}

	}
	cin.get();

	for ( int k = 0; k < inst->n_meta_items; k++){
		cout << "Forward star of meta item \t" << k << "\t size of his neighbourhood\t" << inst->DP[k] << endl;
		for (int  kk = inst->NFS[k]; kk < inst->NFS[k+1]; kk++ )
		{
			cout << "item\t" << inst->AFS[kk] << endl;
		}
	}
	cin.get();

#endif





	for ( int i = 0; i < inst->m_scenarios; i++)
	{

		for(int j=0;j<inst->n_items;j++)
		{

			if(randomBETWEEN_double(0,1)<inst->ProbScenario)
			{
				if(inst->distribute_a==false)
				{
					inst->a[i][j]=inst->value_a;
				}
				else
				{
					//inst->a[i][j]=randomBETWEEN_double(0,inst->value_a);
					inst->a[i][j]=randomBETWEEN(1,inst->value_a);
				}
			}
			else
			{
				inst->a[i][j]=0.0;
			}
		}
	}


#ifdef print_a_coefficients

	cout << "OBJECTIVE FUNCTION COEFFICIENTS:\n";
	for(int j=0;j<inst->n_items;j++)
	{
		cout << "item\t" << j << "\t\t";
		for ( int i = 0; i < inst->m_scenarios; i++)
		{
			cout << inst->a[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	cin.get();

#endif

}

/***************************************************************************/
void generate_data_MP(instance *inst)
/***************************************************************************/
{

	/*
	The exputil problem instances
	Shabbir Ahmed (sahmed@isye.gatech.edu)
	06/17/2013

	This directory consists of 6 instances of an expected utility knapsack
	Instances generated as described in the paper:

	Shabbir Ahmed and Alper Atamturk. Maximizing a class of submodular utility functions.
	Mathematical Programming, 128(1-2):149ñ169, 2011.

	Formulation:
	------------
	m = scenarios
	n = variables
	L = risk coeff


	min sum_(i=1 to m) pi_i ( -1 + exp((1/L)*(sum_j -1*v_ij x_j)))
	s.t. sum_(j=1 to n) a_j x_j <= 1, x_j in {0,1}

	Data generation:
	----------------
	L = 4
	pi_i = 1/m
	a_j = Uniform(0,0.2)
	v_ij = r_ij*a_j
	where
	ln r_ij = alpha_j + beta_j*lnf_i + e_ij
	and
	alpha_j = Uniform(0.05,0.10)
	beta_j = Uniform(0,1)
	lnf_i = Normal(0.05,0.0025)
	e_ij = Normal(0,0.0025)


	Instances and optimal values:
	-----------------------------
	The instances are named as exputil_n_m with n being the number of variables and m the number of scenarios.

	------------------------------
	Instance	 Optimal value
	------------------------------
	exputil_25_100	 -0.242304
	exputil_25_500	 -0.246211
	exputil_25_1000	 -0.242750
	exputil_50_100	 -0.246343
	exputil_50_500	 -0.247751
	exputil_50_1000	 -0.247643
	exputil_100_100	 -0.247542
	exputil_100_500	 -0.248989
	exputil_100_1000 -0.249317
	------------------------------

	The data files:
	---------------
	The data file for each instances (with a .dat extension) is formatted as follows:
	The first line has m and n
	The next n lines list a_j
	The next m*n lines list v_ij

	 */


	cout << "n_items\t" << inst->n_items << endl;
	cout << "n_meta_items\t" << inst->n_meta_items << endl;
	cout << "m_scenarios\t" << inst->m_scenarios << endl;


	inst->weight_item_MP=new double[inst->n_items];


	inst->a=new double*[inst->m_scenarios];
	for ( int i = 0; i < inst->m_scenarios; i++)
	{
		inst->a[i]=new double[inst->n_items];
	}


	// ----->>>>> weight_item_MP = Uniform(0,0.2)

	//std::random_device rd_a;  //Will be used to obtain a seed for the random number engine
	//std::mt19937 gen_a(rd_a()); //Standard mersenne_twister_engine seeded with rd()
	std::mt19937 gen_a(inst->seed);

	std::uniform_real_distribution<> dis_a(0.00, 0.2);
	for ( int j = 0; j < inst->n_items; j++)
	{
		inst->weight_item_MP[j]=dis_a(gen_a);
		//cout << inst->weight_item_MP[j] << endl;
	}


	//	---->>> alpha_j = Uniform(0.05,0.10)
	//	---->>> beta_j = Uniform(0,1)

	double *alpha=new double[inst->n_items] ;
	double *beta =new double[inst->n_items] ;

	//	std::random_device rd_alfa;  //Will be used to obtain a seed for the random number engine
	//	std::mt19937 gen_alfa(rd_alfa()); //Standard mersenne_twister_engine seeded with rd()
	std::mt19937 gen_alfa(inst->seed);

	std::uniform_real_distribution<> dis_alfa(0.05,0.10);

	//	std::random_device rd_beta;  //Will be used to obtain a seed for the random number engine
	//	std::mt19937 gen_beta(rd_beta()); //Standard mersenne_twister_engine seeded with rd()
	std::mt19937 gen_beta(inst->seed);

	std::uniform_real_distribution<> dis_beta(0,1);

	for(int j=0;j<inst->n_items;j++)
	{
		alpha[j] = dis_alfa(gen_a);
		beta[j] = dis_beta(gen_a);
	}

	// -->>>>	lnf_i = Normal(0.05,0.0025)
	//-->>>> e_ij = Normal(0,0.0025)

	std::default_random_engine generator;
	std::normal_distribution<double> normal_lnf(0.05,0.0025);
	std::normal_distribution<double> normal_eps(0,0.0025);

	for ( int i = 0; i < inst->m_scenarios; i++)
	{
		double lnf = normal_lnf(generator);

		for(int j=0;j<inst->n_items;j++)
		{
			double epsilon = normal_eps(generator);

			//---->>>>	a_ij = r_ij*weight_item_MP_j

			//---->>>>	ln r_ij = alpha_j + beta_j*lnf_i + e_ij

			//			cout << "BETA\t" << beta[j] << endl;
			//			cout << "lnf\t" << lnf << endl;
			//			cout << "epsilon\t" << epsilon << endl;

			inst->a[i][j] = exp(alpha[j] + beta[j] * lnf + epsilon)*inst->weight_item_MP[j];
			//cout << inst->a[i][j] << endl;

			if(inst->a[i][j]<-0.0001)
			{
				cout << inst->a[i][j] << "NEGATIVO\n\n";
				exit(-1);
			}

		}
	}

	delete []alpha;
	delete []beta;


	/////////////////////////////////////////////////////////////////////////////////////////
	char FILE_NAME[1000];
	sprintf(FILE_NAME,"MP_INSTANCES/exputil_items%d_scenarios%d_seed%d.dat",inst->n_items,inst->m_scenarios,inst->seed);


	ofstream file_MP;
	file_MP.open(FILE_NAME);

	file_MP << inst->m_scenarios << "\t" << inst->n_items << endl;


	for ( int j = 0; j < inst->n_items; j++)
	{

		file_MP << setprecision(15) << setw(15) << inst->weight_item_MP[j] << endl;
	}

	for ( int i = 0; i < inst->m_scenarios; i++)
	{

		for(int j=0;j<inst->n_items;j++)
		{
			file_MP << setprecision(15) << setw(15)  << inst->a[i][j] << endl;
		}
	}


	ofstream compact_file;
	compact_file.open("info_names.txt", ios::app);
	compact_file << FILE_NAME << endl;
	compact_file.close();

	file_MP.close();
	exit(-1);
	/////////////////////////////////////////////////////////////////////////////////////////

}
