# Makefile for ONCVPSP
#
# Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
# University
#

F90        = ifort
FLINKER    = $(F90)
CC         = icc
FCCPP      = iсс

MKLROOT  = $(ONEAPI_ROOT)/mkl/latest

LIBS = -L$(MKLROOT)/lib/ \
$(MKLROOT)/lib/libmkl_core.a \
$(MKLROOT)/lib/libmkl_intel_lp64.a \
$(MKLROOT)/lib/libmkl_sequential.a \
-lpthread -lm -ldl

OBJS_LIBXC = exc_libxc_stub.o

OBJS_ALL =	aeo.o check_data.o const_basis.o dp3int.o eresid.o exchdl.o \
		exc_off_pbe.o excpzca.o excwig.o ldiracfb.o ldiracfs.o \
		lschfb.o lschfs.o lschkb.o lschpb.o lschpse.o lschps.o \
		lschvkbbe.o lschvkbb.o lschvkbs.o modcore.o optimize.o \
		pspot.o qroots.o run_optimize.o sbf8.o sbf_basis.o \
		sbf_rc_der.o tfapot.o vkboutwf.o vkbphsft.o vout.o vpinteg.o \
		vploc.o wf_rc_der.o

OBJS_NRSR =	fphsft.o fpovlp.o gnu_script.o linout.o psatom.o \
		relatom.o run_config.o run_diag.o run_phsft.o run_plot.o \
		run_vkb.o sratom.o upfout.o wellstate.o

OBJS_SR =	oncvpsp.o

OBJS_NR =	oncvpsp_nr.o

OBJS_R =	oncvpsp_r.o fphsft_r.o fpovlp_r.o gnu_script_r.o linout_r.o \
		psatom_r.o relatom.o renorm_r.o run_config_r.o run_diag_r.o \
		run_diag_sr_so_r.o run_phsft_r.o run_plot_r.o run_vkb_r.o \
		sr_so_r.o upfout_r.o wellstate_r.o

OBJS_OPSP = 	$(OBJS_SR) $(OBJS_NRSR) $(OBJS_ALL) $(OBJS_LIBXC)

OBJS_OPSPNR = 	$(OBJS_SR) $(OBJS_NRSR) $(OBJS_ALL) $(OBJS_LIBXC)

OBJS_OPSPR = 	$(OBJS_R) $(OBJS_ALL) $(OBJS_LIBXC)

all:	oncvpsp oncvpspnr oncvpspr clean

oncvpsp:	$(OBJS_OPSP)
		$(FLINKER) -o oncvpsp.x $(OBJS_OPSP) $(LIBS)

oncvpspnr:	$(OBJS_OPSPNR)
		$(FLINKER) -o oncvpspnr.x $(OBJS_OPSPNR) $(LIBS)

oncvpspr:	$(OBJS_OPSPR)
		$(FLINKER) -o oncvpspr.x $(OBJS_OPSPR) $(LIBS)

.SUFFIXES: .c .o .f .f90 .F90

.F90.o:
	$(FCCPP) $(FCCPPFLAGS) $*.F90  > $*_cpp.f90
	$(F90) $(FFLAGS) $(INC) -o $*.o -c $*_cpp.f90

.f90.o:
	$(F90) $(FFLAGS) $(INC) -c $*.f90

.f.o:
	$(F90) $(FFLAGS) $(INC) -c $*.f

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $*.c

clean:
	 /bin/rm -f *.o  *.d  *.mod *.log *_cpp.f90

