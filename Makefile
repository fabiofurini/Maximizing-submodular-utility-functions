# COMPILER

#IBM cluster compiler and edu215 compiler
CXX= g++ 
#Cygwin compilers
#CXX= /usr/bin/g++ (this was also used on edu215)
#CXX= c:/MinGW/bin/g++
#CXX= c:/MinGW/bin/c++
#CXX= c:/MinGW/bin/gcc
#CXX= c:/MinGW/bin/mingw32-g++
#CXX= c:/MinGW/bin/mingw32-c++
#CXX= c:/MinGW/bin/mingw32-gcc

# FILES

OBJS =	src/Main.o \
		src/MOD_GLOBAL.o \
		src/MOD_OUTER.o \
		src/MOD_LOCAL.o \
		src/BEN_OUTER.o \
		src/GREEDY_ALG.o \
		src/instance_generator.o \
		src/instance_reader.o \
		src/coverage_functions.o \
		src/benders_functions.o \
		src/global_functions.o


# CPLEX VERSION (LIBS and INCLUDE files)

#cplex 12.7 HERE---FABIO
 CPLEXLIBDIR =  /home/fabio/ILOG/CPLEX_Studio_AcademicResearch129/cplex/lib/x86-64_linux/static_pic
 LP_INCLUDE= /home/fabio/ILOG/CPLEX_Studio_AcademicResearch129/cplex/include/ilcplex

#cplex 12.7 HERE---STEFANO
#CPLEXLIBDIR = /opt/ibm/ILOG/CPLEX_Studio1271/cplex/lib/x86-64_linux/static_pic
#LP_INCLUDE= /opt/ibm/ILOG/CPLEX_Studio1271/cplex/include/ilcplex

#CPLEXLIBDIR =  /Users/fabiofurini/Applications/IBM/ILOG/CPLEX_Studio127/cplex/lib/x86-64_osx/static_pic
#LP_INCLUDE= /Users/fabiofurini/Applications/IBM/ILOG/CPLEX_Studio127/cplex/include/ilcplex

#CPLEXLIBDIR =  /media/646bcabe-12bb-400d-87cb-475f895d4667/fabio/ILOG/CPLEX_Studio_AcademicResearch126/cplex/lib/x86-64_linux/static_pic
#LP_INCLUDE= /media/646bcabe-12bb-400d-87cb-475f895d4667/fabio/ILOG/CPLEX_Studio_AcademicResearch126/cplex/include/ilcplex



# Nothing should be changed

#LP_LIBS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread 
LP_LIBS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -lm -lpthread -ldl
#LP_LIBS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -lm -lpthread -liberty -lbfd


#LP_LIBS    = c:/ilog/cplex100/lib/x86_.net2005_8.0/stat_mda/cplex100.lib
#LP_LIBS    = /usr/ilog/cplex100/lib/x86-64_RHEL3.0_3.2/static_pic/libcplex.a
#LP_LIBS    = /usr/ilog/cplex100/lib/x86-64_sles9.0_3.3/static_pic/libcplex.a
#LP_LIBS    = /usr/ilog/cplex100/lib/x86-64_rhel4.0_3.4/static_pic/libcplex.a

#LP_INCLUDE= c:/ilog/cplex100/include/ilcplex
#LP_INCLUDE= /usr/ilog/cplex100/include/ilcplex

DBG= -O3
#DBG= -g


#DEFS = $(OS_VERSION) $(COMPILER) $(LP_SOLVERS)

#INCDIR = -I. -I$(LP_INCLUDE) -I$(CONCORDE_INCLUDE)
INCDIR = -I. -I$(LP_INCLUDE)

#COMPILER FLAGS

#CXXFLAGS =  $(DBG) $(DEFS) $(INCDIR)
#IBM cluster compiler flags
#CXXFLAGS =  $(DBG) $(INCDIR) -mcpu=powerpc64 -maix64 
#edu215 compiler flags
CXXFLAGS =  $(DBG) $(INCDIR)  -std=c++11

.c.o:
	gcc -c $(CXXFLAGS) $< -o $@

#LDLIBS = $(CONCORDE_LIBS) $(MY_LIBS) $(LP_LIBS) 
#LDLIBS = $(CONCORDE_LIBS) $(LP_LIBS) 
LDLIBS = $(LP_LIBS)

all:INFLU

INFLU: $(OBJS)
		$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDLIBS)

$(OBJS): Makefile

clean:
	rm -f $(OBJS) rm INFLU
