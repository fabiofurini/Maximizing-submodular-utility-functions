#ifndef DFL_BEN_2_HEADER
#define DFL_BEN_2_HEADER


#include "global_variables.h"
#include "global_functions.h"
#include "coverage_functions.h"
#include "benders_functions.h"


using namespace std;


/*****************************************************************/
void build_model_BEN(instance *inst);
/*****************************************************************/

/*****************************************************************/
void solve_model_BEN(instance *inst);
/*****************************************************************/

/*****************************************************************/
void solve_model_BEN_LP(instance *inst);
/*****************************************************************/

/*****************************************************************/
void clean_model_BEN(instance *inst);
/*****************************************************************/



#endif
