#RUNLIST GENERATOR

"""
Prints the runlist from global variables and lists (defined after the def. of this function)
"""
def printRunList(TEST_ID):
    for algo in _ALGO:
        for option in _OPTION:
            if algo != 14 and option == -3:
                continue
            tl=_TL
            for n_items in _N_ITEMS:
                for m_scenarios in _M_SCENARIOS:
                    for n_meta_items in _N_META_ITEMS:
                        for _lambda in _LAMBDA:
                            for delta in _DELTA:
                                for seed in _SEED:
                                    TOLL=_TOLL=0
                                    for ProbScenario in _PROBSCENARIO:
                                        for cardinality in _CARDINALITY:
                                            for conflict_perc in _CONFLICT_PERC:
                                                value_a=_VALUE_A
                                                distribute_a=_DISTRIBUTE_A
                                                for KP_constraint in _KP_CONSTRAINT:
                                                    for KP_constraint_R_VALUE in _KP_CONSTRAINT_R_VALUE:
                                                        for KP_constraint_perc_cap in _KP_CONSTRAINT_PERC_CAP:
                                                            USE_WORST_CASE_INSTANCE=_USE_WORST_CASE_INSTANCE
                                                            type_of_zed_function=_TYPE_OF_ZED_FUNCTION
                                                            partition_constraints=_PARTITION_CONSTRAINTS
                                                            for meta_item_per_element in _META_ITEM_PER_ELEMENT:
                                                                for budget_per_element in _BUDGET_PER_ELEMENT:
                                                                    #TEST_ID
                                                                    for scale_factor_alpha in _SCALE_FACTOR_ALPHA:
                                                                        print(exe, instanze, paramFile, algo, option, tl, n_items, m_scenarios, n_meta_items, _lambda, delta, seed, TOLL, ProbScenario, cardinality, conflict_perc, value_a, distribute_a, KP_constraint, KP_constraint_R_VALUE, KP_constraint_perc_cap, USE_WORST_CASE_INSTANCE, type_of_zed_function, partition_constraints, meta_item_per_element, budget_per_element, TEST_ID, scale_factor_alpha, sep='\t')

                                                                        TEST_ID = TEST_ID + 1
    return TEST_ID-1
#end of printRunList(TEST_ID) function

#GENERAL VALUES
exe="./INFLU"
instanze="COVERING"
paramFile="param.txt"
TEST_ID = 0
#######################################################################################################
# CARDINALITY INSTANCES
#######################################################################################################

# _ALGO=[14, 13, 12]
# _OPTION=[-3, 3]
_TL=600
_N_ITEMS=[5000, 10000, 20000]
_M_SCENARIOS=[50, 100]
_N_META_ITEMS=[40, 50, 60]
_LAMBDA=[0.25, 0.5, 0.75, 1, 2]
_DELTA=[2, 4, 5, 6]
_SEED=[1, 2, 3, 4, 5]
_TOLL=0
_PROBSCENARIO=[0.25, 0.5]
_CARDINALITY = [10, 15]
_CONFLICT_PERC=[0, 0.01]
_VALUE_A=10
_DISTRIBUTE_A=1
_KP_CONSTRAINT=[0]
_KP_CONSTRAINT_R_VALUE=[0]
_KP_CONSTRAINT_PERC_CAP=[0]
_USE_WORST_CASE_INSTANCE=0
_TYPE_OF_ZED_FUNCTION=1
_PARTITION_CONSTRAINTS=0
_META_ITEM_PER_ELEMENT=[0]
_BUDGET_PER_ELEMENT=[0]
#TEST_ID
_SCALE_FACTOR_ALPHA=[0]

_ALGO=[14]
_OPTION=[-3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

_ALGO=[14]
_OPTION=[3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

_ALGO=[13]
_OPTION=[3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

_ALGO=[12]
_OPTION=[3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

#######################################################################################################
# KP INSTANCES
#######################################################################################################

# _ALGO=[14, 13, 12]
# _OPTION=[-3, 3]
_TL=600
_N_ITEMS=[5000, 10000, 20000]
_M_SCENARIOS=[50, 100]
_N_META_ITEMS=[40, 50, 60]
_LAMBDA=[0.25, 0.5, 0.75, 1, 2]
_DELTA=[2, 4, 5, 6]
_SEED=[1, 2, 3, 4, 5]
_TOLL=0
_PROBSCENARIO=[0.25, 0.5]
_CARDINALITY = [0]
_CONFLICT_PERC=[0]
_VALUE_A=10
_DISTRIBUTE_A=1
_KP_CONSTRAINT=[1,3]
_KP_CONSTRAINT_R_VALUE=[1000]
_KP_CONSTRAINT_PERC_CAP=[50]
_USE_WORST_CASE_INSTANCE=0
_TYPE_OF_ZED_FUNCTION=1
_PARTITION_CONSTRAINTS=0
_META_ITEM_PER_ELEMENT=[0]
_BUDGET_PER_ELEMENT=[0]
#TEST_ID
_SCALE_FACTOR_ALPHA=[0]

_ALGO=[14]
_OPTION=[-3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

_ALGO=[14]
_OPTION=[3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

_ALGO=[13]
_OPTION=[3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

_ALGO=[12]
_OPTION=[3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

#######################################################################################################
# PARTITION INSTANCES
#######################################################################################################

# _ALGO=[14, 13, 12]
# _OPTION=[-3, 3]
_TL=600
_N_ITEMS=[5000, 10000, 20000]
_M_SCENARIOS=[50, 100]
_N_META_ITEMS=[60]
_LAMBDA=[0.25, 0.5, 0.75, 1, 2]
_DELTA=[2, 4, 5, 6]
_SEED=[1, 2, 3, 4, 5]
_TOLL=0
_PROBSCENARIO=[0.25, 0.5]
_CARDINALITY = [0]
_CONFLICT_PERC=[0, 0.01]
_VALUE_A=10
_DISTRIBUTE_A=1
_KP_CONSTRAINT=[0]
_KP_CONSTRAINT_R_VALUE=[0]
_KP_CONSTRAINT_PERC_CAP=[0]
_USE_WORST_CASE_INSTANCE=0
_TYPE_OF_ZED_FUNCTION=1
_PARTITION_CONSTRAINTS=1
_META_ITEM_PER_ELEMENT=[20]
_BUDGET_PER_ELEMENT=[5]
#TEST_ID
_SCALE_FACTOR_ALPHA=[0]

_ALGO=[14]
_OPTION=[-3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

_ALGO=[14]
_OPTION=[3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

_ALGO=[13]
_OPTION=[3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

_ALGO=[12]
_OPTION=[3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

#######################################################################################################
# CARDINALTY DECREASING INSTANCES
#######################################################################################################

# _ALGO=[14, 13, 12]
# _OPTION=[-3, 3]
_TL=600
_N_ITEMS=[5000, 10000]
_M_SCENARIOS=[50, 100]
_N_META_ITEMS=[40]
_LAMBDA=[0.25, 0.5, 0.75, 1, 2]
_DELTA=[2, 4, 5, 6]
_SEED=[1, 2, 3, 4, 5]
_TOLL=0
_PROBSCENARIO=[0.25, 0.5]
_CARDINALITY = [10]
_CONFLICT_PERC=[0]
_VALUE_A=10
_DISTRIBUTE_A=1
_KP_CONSTRAINT=[0]
_KP_CONSTRAINT_R_VALUE=[0]
_KP_CONSTRAINT_PERC_CAP=[0]
_USE_WORST_CASE_INSTANCE=0
_TYPE_OF_ZED_FUNCTION=3
_PARTITION_CONSTRAINTS=0
_META_ITEM_PER_ELEMENT=[0]
_BUDGET_PER_ELEMENT=[0]
#TEST_ID
_SCALE_FACTOR_ALPHA=[0.5]

_ALGO=[14]
_OPTION=[-3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

_ALGO=[14]
_OPTION=[3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

_ALGO=[13]
_OPTION=[3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)

_ALGO=[12]
_OPTION=[3]
TEST_ID=printRunList(TEST_ID+1)
print("TEST_ID=", TEST_ID)
