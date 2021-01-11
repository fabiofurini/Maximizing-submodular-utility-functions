#ifndef BENDERS_FUNCTIONS_HEADER
#define BENDERS_FUNCTIONS_HEADER

#include "global_variables.h"

#include "global_functions.h"


/*****************************************************************/
void load_I_TILDE(instance *inst, bool rounding,double *I_TILDE_DFL,double *DFL_BEN_2_Y);
/*****************************************************************/

/*****************************************************************/
void comb_solve_model_BEN_AUX_1(instance *inst,double *I_TILDE,double *AUX_SOL, int k);
/*****************************************************************/

/*****************************************************************/
void comb_solve_model_BEN_AUX_2(instance *inst,double *I_TILDE,double *AUX_SOL,int k);
/*****************************************************************/


#endif
