

#include "benders_functions.h"

#define TOLERANCE_SET_TO_ONE 1e-5

/*****************************************************************/
void load_I_TILDE(instance *inst, bool rounding,double *I_TILDE_DFL,double *DFL_BEN_2_Y)
/*****************************************************************/
{


	for ( int j = 0; j < inst->n_items; j++)
	{

		I_TILDE_DFL[j]=0.0;

		for (int  k = inst->NBS[j]; k < inst->NBS[j+1]; k++ )
		{

			I_TILDE_DFL[j]=I_TILDE_DFL[j]+DFL_BEN_2_Y[inst->ABS[k]];
		}

		if(rounding)
		{
			I_TILDE_DFL[j]=(int)(I_TILDE_DFL[j]+0.5);
		}

		if( (1 - TOLERANCE_SET_TO_ONE) < I_TILDE_DFL[j] && I_TILDE_DFL[j] < (1 + TOLERANCE_SET_TO_ONE))
		{
			I_TILDE_DFL[j]=1.0;
		}
	}

}


/*****************************************************************/
void comb_solve_model_BEN_AUX_1(instance *inst,double *I_TILDE,double *AUX_SOL, int k)
/*****************************************************************/
{

	for ( int i = 0; i < inst->n_items; i++ )
	{

		if(I_TILDE[i] < 1 - 0.0001 ||  inst->DM[i]==1)
		{

			//pi
			AUX_SOL[i]=inst->a[k][i];
			//sigma
			AUX_SOL[i+inst->n_items]=0.0;
		}
		else
		{
			//pi
			AUX_SOL[i]=0.0;
			//sigma
			AUX_SOL[i+inst->n_items]=inst->a[k][i];
		}
	}

}


/*****************************************************************/
void comb_solve_model_BEN_AUX_2(instance *inst,double *I_TILDE,double *AUX_SOL,int k)
/*****************************************************************/
{



	for ( int i = 0; i < inst->n_items; i++ )
	{

		if(I_TILDE[i] <= 1 + 0.0001)
		{

			//pi
			AUX_SOL[i]=inst->a[k][i];
			//sigma
			AUX_SOL[i+inst->n_items]=0.0;
		}
		else
		{
			//pi
			AUX_SOL[i]=0.0;
			//sigma
			AUX_SOL[i+inst->n_items]=inst->a[k][i];
		}
	}
}
