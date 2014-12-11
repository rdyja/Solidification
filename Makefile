#############################################################################
FEMLIBPATH = $(HOME)/Programming/Libraries/TalyFEMLib-6.9.10
PETSC_DIR  = $(HOME)/Programming/Libraries/petsc-3.4.5
PETSC_ARCH = linux-hypre-gnu-dbg
##############################################################################
include ${PETSC_DIR}/conf/variables

CC	=	mpicxx
CFLAGS	=	-std=c++0x $(CC_FLAGS)
LINK	=	$(CC)
LFLAGS	=	
INCFLAG = 	-DPETSC_USE_BOPT_O

TALYFEMLIBPATH = -L${FEMLIBPATH} -ltalyfem

INCPATH	= -I${FEMLIBPATH} ${PETSC_CC_INCLUDES}
LIBPATH = $(TALYFEMLIBPATH) $(PETSC_LIB)
OBJS = main.o solid_input_data.o solid_equation.o solid_model.o enthalpy_model.o\
 solid_material.o solid_grid_field.o

LIBS	=  $(LIBPATH)
VPATH=src

####### Targets.

TARGETS	= solidtaly

##### Build rules.
solidtaly: $(OBJS)
	$(CC) -O -o $@ $(OBJS) $(LIBPATH)	

all: solidtaly
	rm *.o

clean:
	rm -f *.o $(TARGETS) core
	
remres:
	rm *.plt

##### Implicit rules.

.SUFFIXES: .cpp .cxx .cc .C .c

.cpp.o:
	$(CC) -c $(CFLAGS) $(INCPATH) $(INCFLAG)  $<	

.cxx.o:
	$(CC) -c $(CFLAGS) $(INCPATH) $(INCFLAG)  $<	

.cc.o:
	$(CC) -c $(CFLAGS) $(INCPATH) $(INCFLAG)  $<	

.C.o:
	$(CC) -c $(CFLAGS) $(INCPATH) $(INCFLAG)  $<	

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) $(INCFLAG)  $<	






depend:
	$(CXX) -MM $(CXXFLAGS) *.cpp > depend.mk



-include depend.mk
