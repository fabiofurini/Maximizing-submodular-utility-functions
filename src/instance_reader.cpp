

#include "instance_reader.h"

/***************************************************************************/
void read_file_MP(instance *inst)
/***************************************************************************/
{

	ifstream in(inst->input_file);
	if(!in)
	{
		cout << "File could not be opened. " << endl;
		exit(1);
	}

	//cout << "READING MP FILE\n";

	int n_meta_items,m_scenarios;

	in >> m_scenarios  ;

	in >> n_meta_items  ;

	inst->n_items=n_meta_items;
	inst->n_meta_items=n_meta_items;
	inst->m_scenarios=m_scenarios;

	cout << "n_items\t" << inst->n_items << endl;
	cout << "n_meta_items\t" << inst->n_meta_items << endl;
	cout << "m_scenarios\t" << inst->m_scenarios << endl;

	inst->weight_item_MP=new double[inst->n_items];

	for ( int j = 0; j < inst->n_items; j++)
	{

		in >> inst->weight_item_MP[j];
	}

	inst->a=new double*[inst->m_scenarios];

	for ( int i = 0; i < inst->m_scenarios; i++)
	{
		inst->a[i]=new double[inst->n_items];
	}

	for ( int i = 0; i < inst->m_scenarios; i++)
	{

		for(int j=0;j<inst->n_items;j++)
		{
			in >> inst->a[i][j];
		}
	}

	cout << "BUILDING neighbourhoods\n";
	vector < vector < int > > dummy;

	for ( int j = 0; j < inst->n_items; j++)
	{
		vector < int > local_dummy;


		local_dummy.push_back(j);

		dummy.push_back(local_dummy);

	}

#ifdef print_neighbourhoods
	for ( int i = 0; i < inst->n_items; i++ ){

		cout << "INVESTIMENT\t" << i << "\t number of SECTORS \t " << dummy[i].size() << "\t\t";

		for ( int j = 0; j <  dummy[i].size(); j++){
			cout << dummy[i][j] << " ";
		}

		cout << endl;
	}
	cout << endl;
#endif


	inst->item_OK=new bool[inst->n_items];

	//////////////////////////////////////////////
	for ( int j = 0; j < inst->n_items; j++ )
	{
		inst->item_OK[j]=true;
		if((int)(dummy[j].size())==0)
		{
			inst->item_OK[j]=false;
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

	inst->NBS[0]=0;
	for ( int j = 0; j < inst->n_items; j++ ){

		inst->DM[j]=dummy[j].size();

		inst->NBS[j+1]=inst->NBS[j]+dummy[j].size();

		for ( int jjj = 0; jjj <  (int)dummy[j].size(); jjj++)
		{

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

		for ( int j = 0; j < inst->n_items; j++ ){
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
			cout << "sector\t" << inst->ABS[k]  << endl;
		}
		cin.get();
	}

	for ( int k = 0; k < inst->n_meta_items; k++){
		cout << "Forward star of meta item \t" << k << "\t size of his neighbourhood\t" << inst->DP[k] << endl;
		for (int  kk = inst->NFS[k]; kk < inst->NFS[k+1]; kk++ )
		{
			cout << "item\t" << inst->AFS[kk] << endl;
		}
		cin.get();
	}
#endif



#ifdef print_item_weights

	cout << "\n\nCAPITAL REQUIREMENTS --->>> ITEM\n";
	for ( int j = 0; j < inst->n_items; j++)
	{
		cout << inst->weight_item_MP[j] << "\t";
	}
	cout << endl;

#endif



#ifdef print_a_coefficients

	for ( int i = 0; i < inst->m_scenarios; i++)
	{
		for(int j=0;j<inst->n_items;j++)
		{
			cout << "a " << i << " " << j << " = " << inst->a[i][j] << endl;
		}
	}

#endif


	in.close();



}
