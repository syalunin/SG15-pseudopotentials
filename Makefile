# Makefile for ONCVPSP
#
# Copyright (c) 1989-2017 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
# University
#
# . /opt/intel/oneapi/setvars.sh
# export PATH=$PATH:/Users/yalunin/install/mpich/4.0.2/bin

MAKE = make

include make.inc

all:
	cd src ; $(MAKE) all
	./set_path
	cd tests/data ; ./TEST.sh

test:
	cd tests/data ; ./TEST.sh

clean:
	cd src ; $(MAKE) clean ; rm *.x
	cd tests/data ; /bin/rm -f *.out *.diff TEST.report

