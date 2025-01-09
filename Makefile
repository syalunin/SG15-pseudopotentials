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
