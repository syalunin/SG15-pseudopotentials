# Run the following commands to initialize oneAPI environment
# . /opt/intel/oneapi/setvars.sh
# export PATH=$PATH:/Users/yalunin/install/mpich/4.0.2/bin

TOPDIR   = $(HOME)/documents/qe-7.3

include $(TOPDIR)/make.inc

MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)$(TOPDIR)/PW/src \
           $(MOD_FLAG)$(TOPDIR)/LR_Modules

PWOBJS   = $(TOPDIR)/PW/src/libpw.a

QEMODS   = $(TOPDIR)/Modules/libqemod.a \
           $(TOPDIR)/upflib/libupf.a \
           $(TOPDIR)/KS_Solvers/libks_solvers.a \
           $(TOPDIR)/FFTXlib/src/libqefft.a \
           $(TOPDIR)/LAXlib/libqela.a \
           $(TOPDIR)/UtilXlib/libutil.a \
           $(TOPDIR)/dft-d3/libdftd3qe.a \
           $(TOPDIR)/PP/src/libpp.a

MODULES  = $(BASEMODS) $(PWOBJS) $(QEMODS)

OBJ     = math.o
TARGET  = main.x 
TEST    = 

all:	fresh $(OBJ) $(TEST) $(TARGET) clean

fresh:
	@rm -f *.x *.o *.mod

$(OBJ):
	$(F90) -c $(basename $@).f90 $(MODFLAGS)

$(TARGET):
	$(F90) -c $(basename $@).f90 -o $(basename $@).o $(MODFLAGS)
	$(LD)  -o $@ $(basename $@).o $(OBJ) $(LDFLAGS) $(MODULES) $(QELIBS)

$(TEST):
	$(F90) -o $@ $(basename $@).f90 math.o $(LAPACK_LIBS)

clean :
	@rm -f *.o *.mod
