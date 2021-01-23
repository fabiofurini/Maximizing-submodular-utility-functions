LIST OF INPUT PARAMETERS:


///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 1: ''input_file''

it is the instance file, we have two modes. The code can read an instance or we can generate one.

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 2: ''param_file''

it is a parameter file in which we have the following parameters:

number_of_CPU 1 
TOLL_DERIVATIVE  -> it determines the tolerance in which the derivative is considered 0 (default 1e-12) 
TOLL_VIOL  -> it determines if a cut (on integer point) is considered violated (default 10e-10) 
TOLL_VIOL_FRAC_BEN  -> this is the violation of the fractional benders cuts (default 10e-6)
FLAG_GREEDY_SOL -> if 1 it load the greedy solution to initialize the models

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 3: ''algorithm''

it is the algorithm, we have 4 algorithms,

11 -> submodular cuts on the entire objective function 
12 -> submodular cuts per scenario 
13 -> outer approximation and then submodular cuts 
14 -> outer approximation and then benders cuts

///////////////////////////////////////////////////////////////////////////////////////////////////////////
*** Parameter 4: ''option''

it determines the cuts to be separated, 1 for the first cut only, 2 for the second cut only and 3 for both.

For the algorithm 14, the code has options -1 -2 and -3 to separate the fractional counterparts of the benders cuts respectively.

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 5: ''timelimit''

the time limit

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 6: ''n_items''

the number of items

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 7: ''m_scenarios''

the number of scenarios

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 8: ''n_meta_items''

the number of meta items.

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 9: ''lambda''

 the value of lambda 

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 10: ''delta''

the covering radius.

The items and the meta items are randomly sampled in a square of value 30. Then the covering radius determines the items that are covered by the meta items, i.e. if the euclidean distance is smaller than the covering radius.

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 11: ''seed''

 the seed for the random number generation

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***parameter 12: ''TOLL_OPTIMALITY''

it determines the optimality gap (for cplex).

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 13: ''ProbScenario''

it is the probability of an item to be in a scenario

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 14: ''cardinality''

it is the RHS value of the cardinality constraint (0 for the other type of constraints)

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 15: ''conflict_perc''

it determines the probability of having a conflict between pairs of meta items (for each conflict the code adds to the models a conflict constraints)

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameters 16 and 17: ''value_a'' and ''distribute_a''

it determines the value of the a_{ij} (or its upper bound). It works combined with parameter 17. In case distribute_a is set to 0 then it works as before, i.e., it sets a_{ij} to value_a in all the scenarios in which the item is present according to the parameter on the probability of an item being in a scenario. Instead if distribute_a is set to 1, a_{ij} is randomly chosen in the interval [1,value_a] and integer. Clearly the probability of an item being in a scenario also set to zero some of the a_{ij}.


///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameters 18: ''KP_constraint''

it decides if a KP constraint is used instead of the cardinality constraints. 

Its value determines the way the weights of the meta-items are determined, we have 5 possibilities (following the classical KP generator):

// 1 -> Uncorrelated: w j u.r. in [1, R]. 
// 2 -> Weakly correlated: w j u.r. in [max{1, p j − R/10}, p j + R/10]. 
// 3 -> Strongly correlated: w j = p j + R/10. 
// 4 -> Almost strongly correlated: w j u.r. in [p j + R/10 − R/500, p j + R/10 + R/500]. 
// 5 -> Subset-sum: w j = p j.

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***parameter 19: ''KP_constraint_R_VALUE''

it is the value of R used above.

As profits of the meta-items, the code uses the average total demand covered (averaged over the scenarios).

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***parameter 20: ''KP_constraint_perc_cap''

it is used to determine the KP constraints as the sum of all item weight multiplied by  KP_constraint_perc_cap/100


///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 21: ''USE_WORST_CASE_INSTANCE''

 if it is set to 1 generate and use the worst case instances, the only input parameter of these instances is the cardinality

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 22: ''type_of_zed_function''

it can be set to 1 to use the lambda based utility function, if set to 2 the code use the identity function, set to  to use the non monotone utility function

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 22,23 and 24: ''partition_constraints'', ''meta_item_per_element'' and ''budget_per_element''


Parameters 22,23 and 24 determine the partition constraints, parameter 22 must be set to 1 to use them instead of the simple cardinality constraint. Then parameter 23 determines the number of meta items per element of the partition (the items are taken in order). 
Finally parameter 24 is the budget, i.e., the rhs of each constraint. So the number of elements in the partition is determined by the number of meta items divided by meta_item_per_element. If the number of meta item is not a multiple of meta_item_per_element the code exits

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 26: ''TEST_ID''

it  is the test ID

///////////////////////////////////////////////////////////////////////////////////////////////////////////
***Parameter 27: ''scale_factor_alpha''

it is the scale factor for alpha (used when type_of_zed_function =3)
