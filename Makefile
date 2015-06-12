include makefile.config
include ${PETSC_DIR}/conf/variables

CC	=	mpicxx
CFLAGS	=	-std=c++0x $(CC_FLAGS)
LINK	=	$(CC)
LFLAGS	=	
INCFLAG = 	-DPETSC_USE_BOPT_O

TALYFEMLIBPATH = -L${FEMLIBPATH} -ltalyfem
JAZ_INCLUDES = ${FEMLIBPATH}/ExternalLibs/jaz


INCPATH	= -I${FEMLIBPATH} ${PETSC_CC_INCLUDES} -I${JAZ_INCLUDES}
LIBPATH = $(TALYFEMLIBPATH) $(PETSC_LIB) -lhdf5
OBJS = main.o solid_input_data.o solid_equation.o solid_model.o enthalpy_model.o\
 solid_material.o solid_grid_field.o extended_input.o contact_bounds.o newton_bc.o

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
