#ifndef DFL_MOD_OUTER_HEADER
#define DFL_MOD_OUTER_HEADER


#include "global_variables.h"
#include "global_functions.h"
#include "coverage_functions.h"
#include "benders_functions.h"


using namespace std;


/*****************************************************************/
void build_model_MOD_OUTER(instance *inst);
/*****************************************************************/

/*****************************************************************/
void solve_model_MOD_OUTER(instance *inst);
/*****************************************************************/

/*****************************************************************/
void clean_model_MOD_OUTER(instance *inst);
/*****************************************************************/



#endif
