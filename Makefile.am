AUTOMAKE_OPTIONS=foreign


WARN_FLAGS= --pedantic -Wall #  -Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimplicit  -Wimplicit-function-declaration  -Wimplicit-int -Wimport  -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wbad-function-cast -Wnested-externs  -Wold-style-definition -Wpointer-sign -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized -Wundef -Wbad-function-cast #   -Wstrict-prototypes   -Wdeclaration-after-statement  -Wmissing-prototypes  
LANG_FLAGS=-std=gnu99 # because of asm 


# TODO they cannot be serious with that...this is the only way of getting rid of the default -O2 -g  
CFLAGS= 
CXXFLAGS=


FEATURE_FLAGS=-D_USE_ALLREDUCE -D_HYBRID -D_MEASURE_TIME


AM_CFLAGS=-D__SIM_SSE3 -D_GNU_SOURCE  -D_OPTIMIZED_FUNCTIONS  $(SIMD_FLAGS)

if WITH_DEBUG
AM_CFLAGS+=-g3 -O0
AM_CFLAGS+= -D_DEBUG
else
AM_CFLAGS+=-O2 -fomit-frame-pointer -funroll-loops -march=native
endif

if PRODUCTIVE
 AM_CFLAGS+= -D_HELLO
else 
 AM_CFLAGS+= -D_NOT_PRODUCTIVE
endif


AM_CFLAGS+=$(FEATURE_FLAGS) $(WARN_FLAGS) $(LANG_FLAGS)
AM_CXXFLAGS=$(AM_CFLAGS)
AM_LDFLAGS = -lm
CC=mpicc

bin_PROGRAMS=exaML
exaML_SOURCES= examl/axml.c  examl/evaluatePartialGenericSpecial.c examl/genericData.c examl/models.c examl/optimizeModel.c examl/searchAlgo.c examl/topologies.c examl/treeIO.c  examl/evaluateGenericSpecial.c  examl/faultTolerance.c examl/makenewzGenericSpecial.c examl/newviewGenericSpecial.c examl/restartHashTable.c examl/thread.c examl/tree.c examl/bipartitionList.c


if WITH_AVX
 exaML_SOURCES +=  examl/avxLikelihood.c
endif
