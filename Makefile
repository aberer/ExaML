VERSIONFILE=examl/version.h
PROGNAME=examl

include system/Makefile.standard

CC:=$(Q) $(CCACHE) mpicc

SRCDIR:=./examl/


# BEGIN customizable stuff 
mode:=
FEATURES :=  -D_USE_ALLREDUCE  -D_OPTIMIZED_FUNCTIONS -D_NOT_PRODUCTIVE #    -D_USE_RTS
WARNING_OFF := -Wno-declaration-after-statement -std=c99  -Wno-sign-compare
# END


DEPS=$(SRCDIR)/axml.h $(SRCDIR)/globalVariables.h
LFLAGS = -lm
OPT :=  -O2 -fomit-frame-pointer -funroll-loops -march=native
WARN := -Wall -pedantic -Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimplicit  -Wimplicit-function-declaration  -Wimplicit-int -Wimport  -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wbad-function-cast  -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs  -Wold-style-definition -Wstrict-prototypes  -Wdeclaration-after-statement -Wpointer-sign -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized -Wundef -Wbad-function-cast  $(WARNING_OFF)
LANG =  -D_GNU_SOURCE $(WARN) 


ifeq ($(mode),hybrid)
FEATURES:=$(FEATURES) -D_USE_PTHREADS
LFLAGS:=$(LFLAGS) -lpthread
endif

SRC_EXCLUDE:=%_flymake.c %_flymake.h %_flymake_master.c
ifeq ($(HAVE_AVX),0)
SRC_EXCLUDE+=%avxLikelihood.c
endif

CFLAGS:= $(VECT_FLAG) $(OPT) $(LANG) $(FEATURES)
DEBUG_CFLAG := $(VECT_FLAG) $(LANG) $(FEATURES)


firstTarget : release 


standardTargets : vupdate cmpMessage depend 

.PHONY : clean vupdate cmpMessage

vupdate : 
	@echo "[VERSION]"
	@echo "#define VERSION \"$(NEW_VERSION)\"" > $(VERSIONFILE)

clean:
	@echo "[CLEAN]"
	$(RM) $(OBJDIR) $(TEST_OBJ_DIR) $(DEBUG_OBJ_DIR) $(PROFILE_OBJ_DIR)  $(ALL_TARGETS)  *~ \#* callgrind* cachegrind*  gmon.out 

 
include system/Makefile.build
include system/Makefile.depend
