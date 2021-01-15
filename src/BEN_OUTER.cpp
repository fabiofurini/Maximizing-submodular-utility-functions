

#include "BEN_OUTER.h"

#define print_solution

//#define write_prob

//#define DEBUG
//#define DEBUG_CUTS

//#define print_rho

/*****************************************************************/
int CPXPUBLIC mycutcallback_BEN(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p)
/*****************************************************************/
{

	(*useraction_p)=CPX_CALLBACK_DEFAULT;

	instance *inst=(instance *) cbhandle;

	int num_variables=inst->n_meta_items+inst->m_scenarios;
	int status;


	status=CPXgetcallbacknodex(env,cbdata,wherefrom,inst->_cut_BEN_Y,0,num_variables-1);
	if(status!=0){
		printf("cannot get the x\n");
		exit(-1);
	}



	load_I_TILDE(inst,true,inst->I_TILDE_BEN,inst->_cut_BEN_Y);

#ifdef 	DEBUG
	cout << "metaitems\n";
	for (int i = 0; i < inst->n_meta_items; i++)
	{
		cout << i << "\t" << inst->_cut_BEN_Y[i] << "\n";
	}
	cout << endl;

	cout << "scenarios\n";
	for (int k = 0; k < inst->m_scenarios; k++)
	{

		cout << k << "\t" << inst->_cut_BEN_Y[inst->n_meta_items+k] << "\n";
	}
	cout << endl;

	cout << "I_TILDE\n";
	for ( int j = 0; j < inst->n_items; j++)
	{

		cout << j << "\t" << inst->I_TILDE_BEN[j] << "\n";
	}
	cout << endl;
	cin.get();
#endif

	for (int k = 0; k < inst->m_scenarios; k++)
	{

		double  punto=compute_subset_coverage_utility_scenario(inst,inst->_cut_BEN_Y,k);

		double denominator=compute_val_diff(inst,punto);

#ifdef 	DEBUG
		cout << "punto\t" << punto << "\t denominator \t" << denominator << endl;
#endif

		if( denominator>-inst->TOLL_DERIVATIVE && denominator<inst->TOLL_DERIVATIVE)
		{
			cout << "NULL DERIVATIVE!";

			int nzcnt=1;

			inst->_cut_BEN_rmatval[0]=1.0;

			inst->_cut_BEN_rmatind[0]=inst->n_meta_items+k;

			inst->_cut_BEN_RHS = punto;

			status=CPXcutcallbackadd (env,cbdata,wherefrom,nzcnt,inst->_cut_BEN_RHS,'L',inst->_cut_BEN_rmatind,inst->_cut_BEN_rmatval,0);
			if(status!=0){
				printf("CPXcutcallbackadd\n");
				exit(-1);
			}

			inst->n_cuts_MOD_LOWER++;

			(*useraction_p)=CPX_CALLBACK_SET;

		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(denominator<-inst->TOLL_DERIVATIVE)
		{



			double nu=inst->_cut_BEN_Y[inst->n_meta_items+k];

			double magic_value= (compute_val_funct(inst,punto) - compute_val_diff(inst,punto) * punto) / compute_val_diff(inst,punto);

			inst->_cut_BEN_RHS=punto;

			if( nu/denominator < punto + magic_value - inst->TOLL_VIOL)
			{

				int nzcnt=inst->n_meta_items+1;

				for (int i = 0; i < inst->n_meta_items; i++)
				{
					if(inst->_cut_BEN_Y[i]<0.5)
					{
						inst->_cut_BEN_rmatval[i]=-inst->_cut_BEN_super_rho[k][i];
						inst->_cut_BEN_rmatind[i]=i;
					}
					else
					{
						inst->_cut_BEN_RHS-=inst->_cut_BEN_single_rho[k][i];
						inst->_cut_BEN_rmatval[i]=-inst->_cut_BEN_single_rho[k][i];
						inst->_cut_BEN_rmatind[i]=i;
					}
				}

				inst->_cut_BEN_rmatval[inst->n_meta_items]=1.0/denominator;

				inst->_cut_BEN_rmatind[inst->n_meta_items]=inst->n_meta_items+k;

				inst->_cut_BEN_RHS += magic_value;

#ifdef 	DEBUG_CUTS
		cout << "***MOD_LOWER***\n";
		for (int i = 0; i < inst->n_meta_items; i++)
		{
			cout << "META-ITEM\t" << i << "\t _cut_BEN_rmatind \t" << inst->_cut_BEN_rmatind[i] << "\t \t"<< inst->_cut_BEN_rmatval[i] << endl;
		}
		cout << "SCENARIO\t" << k << "\t _cut_BEN_rmatind \t" << inst->_cut_BEN_rmatind[inst->n_meta_items] << "\t \t"<< inst->_cut_BEN_rmatval[inst->n_meta_items] << endl;
		cout << "RHS\t" << inst->_cut_BEN_RHS << endl;
		cin.get();
#endif

				status=CPXcutcallbackadd (env,cbdata,wherefrom,nzcnt,inst->_cut_BEN_RHS,'G',inst->_cut_BEN_rmatind,inst->_cut_BEN_rmatval,0);
				if(status!=0){
					printf("CPXcutcallbackadd\n");
					exit(-1);
				}

//				cout << "MOD_LOWER ADDED CUT\n\n";
				//				cin.get();

				inst->n_cuts_MOD_LOWER++;

				(*useraction_p)=CPX_CALLBACK_SET;

			}

		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(denominator>inst->TOLL_DERIVATIVE)
		{

			double magic_value= (compute_val_funct(inst,punto) - compute_val_diff(inst,punto) * punto) ;

			for (int llll = 0; llll < 2; llll++)
			{
				if(inst->option==1 && llll==1){continue;}
				if(inst->option==2 && llll==0){continue;}

				if(llll ==0)
				{
					comb_solve_model_BEN_AUX_1(inst,inst->I_TILDE_BEN,inst->AUX_SOL_BEN,k);
				}

				if(llll ==1)
				{
					comb_solve_model_BEN_AUX_2(inst,inst->I_TILDE_BEN,inst->AUX_SOL_BEN,k);
				}

				int nzcnt=inst->n_meta_items+1;
				double first_part=0;

				inst->_cut_BEN_RHS=0.0;

				for ( int j = 0; j < inst->n_items; j++)
				{
					inst->_cut_BEN_RHS +=	inst->AUX_SOL_BEN[inst->n_items+j];
				}

				for (int i = 0; i <inst->n_meta_items; i++)
				{

					inst->_cut_BEN_rmatval[i]=0.0;

					for (int  k = inst->NFS[i]; k < inst->NFS[i+1]; k++ )
					{
						inst->_cut_BEN_rmatval[i] -= (inst->AUX_SOL_BEN[inst->AFS[k]]*denominator);
					}

					first_part -= (inst->_cut_BEN_rmatval[i] *inst->_cut_BEN_Y[i]);

					inst->_cut_BEN_rmatind[i]=i;
				}


				inst->_cut_BEN_rmatval[inst->n_meta_items]=1.0;

				inst->_cut_BEN_rmatind[inst->n_meta_items]=inst->n_meta_items+k;

				inst->_cut_BEN_RHS*=denominator;

				inst->_cut_BEN_RHS+=magic_value;

#ifdef 	DEBUG_CUTS
		cout << "***BEN_UPPER***\n";
		for (int i = 0; i < inst->n_meta_items; i++)
		{
			cout << "META-ITEM\t" << i << "\t _cut_BEN_rmatind \t" << inst->_cut_BEN_rmatind[i] << "\t \t"<< inst->_cut_BEN_rmatval[i] << endl;
		}
		cout << "SCENARIO\t" << k << "\t _cut_BEN_rmatind \t" << inst->_cut_BEN_rmatind[inst->n_meta_items] << "\t \t"<< inst->_cut_BEN_rmatval[inst->n_meta_items] << endl;
		cout << "RHS\t" << inst->_cut_BEN_RHS << endl;
		cin.get();
#endif

				if(inst->_cut_BEN_Y[inst->n_meta_items+k] > (inst->_cut_BEN_RHS + first_part) + inst->TOLL_VIOL)
				{

					status=CPXcutcallbackadd (env,cbdata,wherefrom,nzcnt,inst->_cut_BEN_RHS,'L',inst->_cut_BEN_rmatind,inst->_cut_BEN_rmatval,0);
					if(status!=0)
					{
						printf("CPXcutcallbackadd\n");
						exit(-1);
					}

					if(llll ==0)
					{
						inst->n_cuts_BEN_1++;
					}

					if(llll ==1)
					{
						inst->n_cuts_BEN_2++;
					}


//					cout << "BEN ADDED CUT\n\n";
//					cin.get();

					(*useraction_p)=CPX_CALLBACK_SET;


				}
			}
		}

	}

	return 0;
}


/*****************************************************************/
int CPXPUBLIC myusercutcallback_BEN(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p)
/*****************************************************************/
{



	(*useraction_p)=CPX_CALLBACK_DEFAULT;

	instance *inst=(instance *) cbhandle;

	int num_variables=inst->n_meta_items+inst->m_scenarios;
	int status;


	status=CPXgetcallbacknodex(env,cbdata,wherefrom,inst->_cut_BEN_Y,0,num_variables-1);
	if(status!=0){
		printf("cannot get the x\n");
		exit(-1);
	}

	load_I_TILDE(inst,false,inst->I_TILDE_BEN,inst->_cut_BEN_Y);

	for (int k = 0; k < inst->m_scenarios; k++)
	{

		double  punto=compute_subset_coverage_utility_scenario_fract(inst,inst->_cut_BEN_Y,k);

		double denominator=compute_val_diff(inst,punto);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(denominator<-inst->TOLL_DERIVATIVE)
		{
			return 0;
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(denominator>inst->TOLL_DERIVATIVE)
		{

			double magic_value= (compute_val_funct(inst,punto) - compute_val_diff(inst,punto) * punto) ;


			for (int llll = 0; llll < 2; llll++)
			{

				if(inst->option==1 && llll==1){continue;}
				if(inst->option==2 && llll==0){continue;}

				if(llll ==0)
				{
					comb_solve_model_BEN_AUX_1(inst,inst->I_TILDE_BEN,inst->AUX_SOL_BEN,k);
				}

				if(llll ==1)
				{
					comb_solve_model_BEN_AUX_2(inst,inst->I_TILDE_BEN,inst->AUX_SOL_BEN,k);
				}

				int nzcnt=inst->n_meta_items+1;
				double first_part=0;

				inst->_cut_BEN_RHS=0.0;

				for ( int j = 0; j < inst->n_items; j++)
				{
					inst->_cut_BEN_RHS +=	inst->AUX_SOL_BEN[inst->n_items+j];
				}

				for (int i = 0; i <inst->n_meta_items; i++)
				{

					inst->_cut_BEN_rmatval[i]=0.0;

					for (int  k = inst->NFS[i]; k < inst->NFS[i+1]; k++ )
					{
						inst->_cut_BEN_rmatval[i] -= (inst->AUX_SOL_BEN[inst->AFS[k]]*denominator);
					}

					first_part -= (inst->_cut_BEN_rmatval[i] *inst->_cut_BEN_Y[i]);

					inst->_cut_BEN_rmatind[i]=i;
				}

				inst->_cut_BEN_rmatval[inst->n_meta_items]=1.0;

				inst->_cut_BEN_rmatind[inst->n_meta_items]=inst->n_meta_items+k;

				inst->_cut_BEN_RHS*=denominator;

				inst->_cut_BEN_RHS+=magic_value;

				if(inst->_cut_BEN_Y[inst->n_meta_items+k] > (inst->_cut_BEN_RHS + first_part) + inst->TOLL_VIOL_FRAC_BEN)
				{

					status=CPXcutcallbackadd (env,cbdata,wherefrom,nzcnt,inst->_cut_BEN_RHS,'L',inst->_cut_BEN_rmatind,inst->_cut_BEN_rmatval,0);
					if(status!=0)
					{
						printf("CPXcutcallbackadd\n");
						exit(-1);
					}

					if(llll ==0)
					{
						inst->n_cuts_BEN_FRAC_1++;
					}

					if(llll ==1)
					{
						inst->n_cuts_BEN_FRAC_2++;
					}


					(*useraction_p)=CPX_CALLBACK_SET;


				}
			}
		}

	}

	return 0;
}




/*****************************************************************/
void build_model_BEN(instance *inst)
/*****************************************************************/
{

	///////////////////////////////////////////////////////////////////////////////

	inst->n_cuts_BEN_1=0;
	inst->n_cuts_BEN_FRAC_1=0;
	inst->n_cuts_BEN_2=0;
	inst->n_cuts_BEN_FRAC_2=0;
	inst->n_cuts_MOD_LOWER=0;

	inst->_cut_BEN_rmatval=new double[inst->n_meta_items+1];
	inst->_cut_BEN_Y=new double[inst->n_items+inst->m_scenarios];
	inst->_cut_BEN_rmatind=new int[inst->n_meta_items+1];
	inst->_cut_BEN_local_rho=new double[inst->n_meta_items+1];


	inst->_cut_BEN_super_rho=new double*[inst->m_scenarios];
	inst->_cut_BEN_single_rho=new double*[inst->m_scenarios];

	for ( int k = 0; k < inst->m_scenarios; k++ )
	{
		inst->_cut_BEN_super_rho[k]=new double[inst->n_meta_items];
		inst->_cut_BEN_single_rho[k]=new double[inst->n_meta_items];
	}
	///////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * setting the CPLEX environment

	//opening the environment
	inst->env_BEN=CPXopenCPLEX(&(inst->status));
	if(inst->status!=0)
	{
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	//opening the pointer to the problem
	inst->lp_BEN=CPXcreateprob(inst->env_BEN,&(inst->status),"BEN");
	if(inst->status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * creating the variables *
	inst->ccnt=inst->n_meta_items+inst->m_scenarios;
	inst->obj=(double*) calloc(inst->ccnt,sizeof(double));
	inst->lb=(double*) calloc(inst->ccnt,sizeof(double));
	inst->ub=(double*) calloc(inst->ccnt,sizeof(double));
	inst->c_type=(char*) calloc(inst->ccnt,sizeof(char));


	inst->colname=(char**) calloc(inst->ccnt,sizeof(char*));
	for(int i=0;i<inst->ccnt;i++){inst->colname[i]=(char*) calloc(2000,sizeof(char));}

	int counter=0;
	for ( int j = 0; j < inst->n_meta_items; j++){

		inst->obj[counter]=0.0;
		inst->lb[counter]=0.0;
		inst->ub[counter]=1.0;
		inst->c_type[counter]='B';
		sprintf(inst->colname[counter], "y%d",j);
		counter++;

	}

	for ( int i = 0; i < inst->m_scenarios; i++)
	{
		inst->obj[counter]=1.0;
		inst->lb[counter]=0.0;

		switch (inst->type_of_zed_function)
		{

		case 1 :
			inst->ub[counter]=1.0;
			break;

		case 2 :
			inst->ub[counter]=CPX_INFBOUND;
			break;

		case 3 :
			inst->ub[counter]=1.0;
			break;

		default :
			cout << "wrong value for type_of_zed_function";
			exit(-1);

		}


		inst->c_type[counter]='C';
		sprintf(inst->colname[counter], "nu%d",i);
		counter++;
	}


	inst->status=CPXnewcols(inst->env_BEN,inst->lp_BEN,inst->ccnt,inst->obj,inst->lb,inst->ub,inst->c_type,inst->colname);
	if(inst->status!=0)
	{
		printf("error in CPXnewcols\n");
		exit(-1);
	}

	free(inst->obj);
	free(inst->lb);
	free(inst->ub);
	free(inst->c_type);

	for(int i=0;i<inst->ccnt;i++){free(inst->colname[i]);}
	free(inst->colname);


	// * setting the objective function in the minimization form
	CPXchgobjsen(inst->env_BEN,inst->lp_BEN,CPX_MAX);



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(inst->FLAG_INSTANCE_MP==true)
	{
		inst->rcnt=1;
		inst->nzcnt=inst->n_meta_items;
		inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
		inst->sense=(char*) calloc(inst->rcnt,sizeof(double));
		inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
		inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
		inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

		inst->rhs[0]=1.0;
		inst->sense[0]='L';

		for ( int i = 0; i < inst->n_meta_items; i++ )
		{
			inst->rmatval[i]=inst->weight_item_MP[i];
			inst->rmatind[i]=i;
		}

		inst->rmatbeg[0]=0;

		inst->status=CPXaddrows(inst->env_BEN,inst->lp_BEN,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
		if(inst->status!=0)
		{
			printf("error in CPXaddrows\n");
			exit(-1);
		}

		free(inst->rmatbeg);
		free(inst->rmatval);
		free(inst->rmatind);
		free(inst->rhs);
		free(inst->sense);

	}



	if(inst->FLAG_INSTANCE_LOCATION==true || inst->FLAG_INSTANCE_WORST_CASE==true)
	{

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(inst->KP_constraint>0)
		{

			inst->rcnt=1;
			inst->nzcnt=inst->n_meta_items;
			inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
			inst->sense=(char*) calloc(inst->rcnt,sizeof(double));
			inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
			inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
			inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

			inst->rhs[0]=inst->KP_constraint_CAPACITY;
			inst->sense[0]='L';

			for ( int i = 0; i < inst->n_meta_items; i++ )
			{
				inst->rmatval[i]=inst->KP_constraint_weights[i];
				inst->rmatind[i]=i;
			}

			inst->rmatbeg[0]=0;

			inst->status=CPXaddrows(inst->env_BEN,inst->lp_BEN,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
			if(inst->status!=0)
			{
				printf("error in CPXaddrows\n");
				exit(-1);
			}

			free(inst->rmatbeg);
			free(inst->rmatval);
			free(inst->rmatind);
			free(inst->rhs);
			free(inst->sense);
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(inst->cardinality!=0 && inst->partition_constraints==0)
		{

			inst->rcnt=1;
			inst->nzcnt=inst->n_meta_items;
			inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
			inst->sense=(char*) calloc(inst->rcnt,sizeof(double));
			inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
			inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
			inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

			inst->rhs[0]=inst->cardinality;
			inst->sense[0]='L';

			for ( int i = 0; i < inst->n_meta_items; i++ )
			{
				inst->rmatval[i]=1.0;
				inst->rmatind[i]=i;
			}


			inst->rmatbeg[0]=0;

			inst->status=CPXaddrows(inst->env_BEN,inst->lp_BEN,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
			if(inst->status!=0)
			{
				printf("error in CPXaddrows\n");
				exit(-1);
			}

			free(inst->rmatbeg);
			free(inst->rmatval);
			free(inst->rmatind);
			free(inst->rhs);
			free(inst->sense);

		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(inst->partition_constraints!=0)
		{

			for ( int p = 0; p < inst->n_element_partition; p++ )
			{

				inst->rcnt=1;
				inst->nzcnt=inst->n_meta_items;
				inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
				inst->sense=(char*) calloc(inst->rcnt,sizeof(double));
				inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
				inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
				inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

				inst->rhs[0]=inst->budget_per_element;

				inst->sense[0]='L';

				int counter=0;
				for ( int kk = 0; kk < inst->n_meta_items; kk++ )
				{
					if(inst->meta_item_element_partition[kk] == p)
					{
						inst->rmatval[counter]=1.0;
						inst->rmatind[counter++]=kk;
					}
				}

				if(inst->budget_per_element>=counter)
				{
					cout << "not enough items in element\t" << p << endl;
					exit(-1);
				}

				inst->nzcnt=counter;

				inst->rmatbeg[0]=0;

				inst->status=CPXaddrows(inst->env_BEN,inst->lp_BEN,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
				if(inst->status!=0)
				{
					printf("error in CPXaddrows\n");
					exit(-1);
				}

				free(inst->rmatbeg);
				free(inst->rmatval);
				free(inst->rmatind);
				free(inst->rhs);
				free(inst->sense);
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(int i=0; i<inst->n_meta_items; i++)
	{
		for(int j=i+1; j<inst->n_meta_items; j++)
		{
			if(inst->CONF_MATRIX[i][j]==1)
			{

				//cout << "item conflict\t" << i << "\t" << j << endl;

				inst->rcnt=1;
				inst->nzcnt=2;
				inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
				inst->sense=(char*) calloc(inst->rcnt,sizeof(double));

				inst->rhs[0]=1.0;
				inst->sense[0]='L';


				inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
				inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
				inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

				inst->rmatval[0]=1.0;
				inst->rmatind[0]=i;

				inst->rmatval[1]=1.0;
				inst->rmatind[1]=j;

				inst->rmatbeg[0]=0;

				inst->status=CPXaddrows(inst->env_BEN,inst->lp_BEN,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
				if(inst->status!=0)
				{
					printf("error in CPXaddrows\n");
					exit(-1);
				}

				free(inst->rmatbeg);
				free(inst->rmatval);
				free(inst->rmatind);
				free(inst->rhs);
				free(inst->sense);

			}
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////




	inst->AUX_SOL_BEN=(double*) calloc(2*inst->n_items+1,sizeof(double));
	inst->I_TILDE_BEN=(double*) calloc(inst->n_items,sizeof(double));


	////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "MODEL BUILD (compute_super_rho_scenario and compute_single_rho_scenario)\n";

	for ( int k = 0; k < inst->m_scenarios; k++ )
	{

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

			inst->_cut_BEN_super_rho[k][j]=dummy;
		}


		for ( int j = 0; j < inst->n_meta_items; j++)
		{

			double dummy=0;
			for (int i = inst->NFS[j]; i < inst->NFS[j+1]; i++ )
			{
				dummy+=inst->a[k][inst->AFS[i]];
			}

			inst->_cut_BEN_single_rho[k][j]=dummy;
		}

	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef print_rho
	for ( int k = 0; k < inst->m_scenarios; k++ )
	{
		cout << "SCENARIO\t" << k << endl;
		for ( int j = 0; j < inst->n_meta_items; j++)
		{
			cout << "META-ITEM\t" << j << "\t single_rho \t"<< inst->_cut_BEN_single_rho[k][j] << "\t super_rho \t" << inst->_cut_BEN_super_rho[k][j] << endl;
		}
	}
	cout << "RO COMPUTED\n";
#endif

#ifdef write_prob
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * writing the created ILP model on a file *
	inst->status=CPXwriteprob(inst->env_BEN,inst->lp_BEN,"BEN.lp",NULL);
	if(inst->status!=0) {
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
	cout << "write BEN.lp\n\n";
	exit(-1);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif

}



/*****************************************************************/
void solve_model_BEN(instance *inst)
/*****************************************************************/
{




	CPXsetintparam (inst->env_BEN, CPX_PARAM_SCRIND, CPX_ON);


	//	// * Set relative tolerance *
	//	inst->status = CPXsetdblparam (inst->env_BEN, CPX_PARAM_EPAGAP, 0.0);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPAGAP\n");
	//	}

	// * Set a tolerance *
	inst->status = CPXsetdblparam (inst->env_BEN, CPX_PARAM_EPGAP, inst->TOLL_OPTIMALITY);
	if (inst->status)
	{
		printf ("error for CPX_PARAM_EPGAP\n");
	}

	//	// * Set mip tolerances integrality *
	//	inst->status = CPXsetdblparam (inst->env_BEN, CPX_PARAM_EPINT, 0.0);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPINTP\n");
	//	}
	//


	//this is an important parameter to compute the right objective function value

	//	// * Set Feasibility tolerance *
	inst->status = CPXsetdblparam (inst->env_BEN, CPX_PARAM_EPRHS, 1e-9);
	if (inst->status)
	{
		printf ("error for CPX_PARAM_EPRHS\n");
	}

	// * Set number of CPU*
	inst->status = CPXsetintparam (inst->env_BEN, CPX_PARAM_THREADS, inst->number_of_CPU);
	if (inst->status)
	{
		printf ("error for CPX_PARAM_EPRHS\n");
	}

	// * Set time limit *
	inst->status = CPXsetdblparam (inst->env_BEN, CPX_PARAM_TILIM,inst->timelimit);
	if (inst->status)
	{
		printf ("error for CPX_PARAM_EPRHS\n");
	}


	///////////////////////////////////////////////////////////////////////////////////
	// * solving the MIP model
	clock_t time_start=clock();


	CPXsetintparam(inst->env_BEN, CPX_PARAM_MIPCBREDLP, CPX_OFF);        // let MIP callbacks work on the original model
	CPXsetintparam(inst->env_BEN, CPX_PARAM_PRELINEAR, CPX_OFF);              // assure linear mappings between the presolved and original models
	CPXsetintparam(inst->env_BEN, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);

	//	memory available for working storage: WorkMem in Concert Technology or CPX_PARAM_WORKMEM in the Callable Library
	//
	//	node storage file switch: NodeFileInd in Concert Technology or CPX_PARAM_NODEFILEIND in the Callable Library
	//
	//	tree memory limit: TreLim in Concert Technology or CPX_PARAM_TRELIM in the Callable Library
	//
	//	directory for working files WorkDir in Concert Technology or CPX_PARAM_WORKDIR in the Callable Library

	CPXsetdblparam(inst->env_BEN,CPX_PARAM_TRELIM,6000);


	inst->status = CPXsetlazyconstraintcallbackfunc(inst->env_BEN,mycutcallback_BEN,inst);
	if (inst->status)
	{
		printf ("error for CPXsetlazyconstraintcallbackfunc\n");
	}


	int option_val=inst->option;

	if(inst->option<0)
	{

		inst->option=abs(inst->option);

		cout << "******FRACTIONAL SEPARATION******\n\n";

		inst->status = CPXsetusercutcallbackfunc(inst->env_BEN,myusercutcallback_BEN,inst);
		if (inst->status)
		{
			printf ("error for CPXsetuserconstraintcallbackfunc\n");
		}

	}


	cout << "\nCPXmipopt:\n";
	inst->status=CPXmipopt(inst->env_BEN,inst->lp_BEN);
	if(inst->status!=0)
	{
		printf("error in CPXmipopt\n");
		//exit(-1);
	}

	inst->option=option_val;


	clock_t time_end=clock();
	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;
	///////////////////////////////////////////////////////////////////////////////////


	bool sol_found=true;

	// * getting the solution
	inst->x=(double*) calloc(inst->n_meta_items+inst->m_scenarios,sizeof(double));


	inst->status=CPXgetmipx(inst->env_BEN,inst->lp_BEN,inst->x,0,inst->n_meta_items+inst->m_scenarios-1);
	if(inst->status!=0)
	{
		sol_found=false;
		printf("error in CPXgetmipx\n");
	}

	inst->objval=-1;
	inst->status=CPXgetmipobjval(inst->env_BEN,inst->lp_BEN,&(inst->objval));
	if(inst->status!=0)
	{
		sol_found=false;
		printf("error in CPXgetmipobjval\n");
	}


	cout << fixed << "\n***objval:\t" << inst->objval << endl;


	int open_facilities=-1;
	int satisfied_clients=-1;

	if(sol_found)
	{

		open_facilities=0;
		satisfied_clients=0;

		for ( int j = 0; j < inst->n_meta_items; j++)
		{
			if( (int)(inst->x[j]+0.5) ==1)
			{
				open_facilities++;
			}
		}


#ifdef print_solution
		printf("\n\nSolution\n");
		for ( int j = 0; j < inst->n_meta_items; j++)
		{

			cout << "META ITEMS\t" << j << "\tval:\t" << (int)(inst->x[j]+0.5) << endl;
		}
		printf("\n");

		for ( int j = 0; j < inst->m_scenarios; j++)
		{

			cout << "SCENARIO\t" << j << "\tval:\t" << inst->x[inst->n_meta_items+j] << endl;
		}
		printf("\n");
#endif

	}


	inst->bestobjval=-1;
	inst->status=CPXgetbestobjval(inst->env_BEN,inst->lp_BEN,&(inst->bestobjval));
	if(inst->status!=0)
	{
		printf("error in CPXgetbestobjval\n");
	}

	inst->lpstat=CPXgetstat(inst->env_BEN,inst->lp_BEN);
	inst->nodecount = CPXgetnodecnt(inst->env_BEN, inst->lp_BEN);

	cout << "\nlpstat\t" << inst->lpstat << endl;
	cout << "\nnodecount\t" << inst->nodecount << endl;

	bool *reached_clients=new bool[inst->n_items];

	satisfied_clients=compute_subset_coverage_number(inst,inst->x,reached_clients);

	///////////////////////////////////////////////////////////////////////////////////////////
#ifdef PRINT_FIGURE_TEX
	if(sol_found)
	{
		if(inst->FLAG_INSTANCE_LOCATION==true)
		{
			int *zz=new int[inst->n_items];
			int *yy=new int[inst->n_meta_items];
			for ( int j = 0; j < inst->n_meta_items; j++)
			{
				yy[j]= (int)(inst->x[j]+0.5);
			}
			for ( int i = 0; i < inst->n_items; i++ )
			{
				if(reached_clients[i])
				{
					zz[i]=1;
				}
				else
				{
					zz[i]=0;
				}
			}

			draw_grid_sol_market(inst,yy,zz);

			delete[]zz;
			delete[]yy;
		}
	}
#endif
	///////////////////////////////////////////////////////////////////////////////////////////

	delete[] reached_clients;

	double recomputed_val=compute_subset_coverage_utility(inst,inst->x);

	cout << fixed << "\n***RECOMPUTED SOL:\t" << recomputed_val << endl;

	cout << "\n***open_facilities\t" << open_facilities << endl;
	cout << "***satisfied_clients\t" << satisfied_clients << endl;


	///////////////////////////////////////////////////////////////////////////////////

	int cur_numcols=CPXgetnumcols(inst->env_BEN,inst->lp_BEN);
	int cur_numrows=CPXgetnumrows(inst->env_BEN,inst->lp_BEN);

	cout <<"n_cuts_BEN_1\t" <<  inst->n_cuts_BEN_1 << "\n";
	cout <<"n_cuts_BEN_FRAC_1\t"<<  inst->n_cuts_BEN_FRAC_1 << "\n";
	cout <<"n_cuts_BEN_2\t" <<  inst->n_cuts_BEN_2 << "\n";
	cout <<"n_cuts_BEN_FRAC_2\t"<<  inst->n_cuts_BEN_FRAC_2 << "\n";

	cout <<"n_cuts_MOD_LOWER\t"<<  inst->n_cuts_MOD_LOWER << "\n";


	//////////////////////////////////////////////////////////////////////////////////
	ofstream compact_file;
	compact_file.open("info_EXTENSIVE.txt", ios::app);
	compact_file << fixed

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
			<<  inst->conflict_perc << "\t"

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

			<<  inst->objval<< "\t"
			<<  inst->bestobjval<< "\t"
			<<  inst->lpstat<< "\t"
			<<   solution_time << "\t"

			<<  recomputed_val << "\t"
			<<  inst->nodecount<< "\t"

			<<  cur_numcols << "\t"
			<<  cur_numrows << "\t"

			<<  open_facilities<< "\t"
			<<  satisfied_clients<< "\t"

			<<  inst->item_not_covered<< "\t"
			<<  inst->item_single_covered<< "\t"

			<<  inst->n_cuts_BEN_1 << "\t"
			<<  inst->n_cuts_BEN_FRAC_1 << "\t"
			<<  inst->n_cuts_BEN_2 << "\t"
			<<  inst->n_cuts_BEN_FRAC_2 << "\t"

			<<  inst->n_cuts_MOD_LOWER << "\t"

			<<  inst->n_conflicts << "\t"

			<<   inst->alpha << "\t"

			<<    inst->AVERAGE_DEMAND


			<< endl;
	compact_file.close();
	//////////////////////////////////////////////////////////////////////////////////

#ifdef PRINT_SINGLE_FILE_OUTPUT
	//////////////////////////////////////////////////////////////////////////////////
	char DUMMY_NAME[1000];

	sprintf(DUMMY_NAME,"SOL/OUTER_BENDERS/TEST_ID_%d.sol",inst->TEST_ID);

	ofstream compact_file_single;
	compact_file_single.open(DUMMY_NAME);
	compact_file_single << fixed

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
			<<  inst->conflict_perc << "\t"

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

			<<  inst->objval<< "\t"
			<<  inst->bestobjval<< "\t"
			<<  inst->lpstat<< "\t"
			<<   solution_time << "\t"

			<<  recomputed_val << "\t"
			<<  inst->nodecount<< "\t"

			<<  cur_numcols << "\t"
			<<  cur_numrows << "\t"

			<<  open_facilities<< "\t"
			<<  satisfied_clients<< "\t"

			<<  inst->item_not_covered<< "\t"
			<<  inst->item_single_covered<< "\t"

			<<  inst->n_cuts_BEN_1 << "\t"
			<<  inst->n_cuts_BEN_FRAC_1 << "\t"
			<<  inst->n_cuts_BEN_2 << "\t"
			<<  inst->n_cuts_BEN_FRAC_2 << "\t"

			<<  inst->n_cuts_MOD_LOWER << "\t"

			<<  inst->n_conflicts << "\t"

			<<   inst->alpha << "\t"

			<<    inst->AVERAGE_DEMAND

			<< endl;
	compact_file_single.close();
	//////////////////////////////////////////////////////////////////////////////////
#endif




	free(inst->x);


}


/*****************************************************************/
void clean_model_BEN(instance *inst)
/*****************************************************************/
{

	inst->status=CPXfreeprob(inst->env_BEN,&(inst->lp_BEN));
	if(inst->status!=0) {printf("error in CPXfreeprob\n");exit(-1);}

	inst->status=CPXcloseCPLEX(&(inst->env_BEN));
	if(inst->status!=0) {printf("error in CPXcloseCPLEX\n");exit(-1);}

	free(inst->AUX_SOL_BEN);
	free(inst->I_TILDE_BEN);

	delete []inst->_cut_BEN_rmatval;
	delete []inst->_cut_BEN_Y;
	delete []inst->_cut_BEN_rmatind;
	delete []inst->_cut_BEN_local_rho;

	for ( int k = 0; k < inst->m_scenarios; k++ )
	{
		delete[]inst->_cut_BEN_super_rho[k];
		delete[]inst->_cut_BEN_single_rho[k];
	}
	delete[]inst->_cut_BEN_super_rho;
	delete[]inst->_cut_BEN_single_rho;

}
