#ifndef DFL_MOD_LOCAL_HEADER
#define DFL_MOD_LOCAL_HEADER


#include "global_variables.h"
#include "global_functions.h"
#include "coverage_functions.h"
#include "benders_functions.h"


using namespace std;


/*****************************************************************/
void build_model_MOD_LOCAL(instance *inst);
/*****************************************************************/

/*****************************************************************/
void solve_model_MOD_LOCAL(instance *inst);
/*****************************************************************/

/*****************************************************************/
void clean_model_MOD_LOCAL(instance *inst);
/*****************************************************************/



#endif
