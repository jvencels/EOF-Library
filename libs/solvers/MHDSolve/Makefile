SOLVER  = MHDSolve
SRC     = MHDSolve.F90
LIBDIR  = /usr/local/share/elmersolver/lib

all:
	elmerf90 -shared $(SRC) -o lib$(SOLVER).so
	sudo cp lib$(SOLVER).so $(LIBDIR)/$(SOLVER).so

clean:
	rm -f lib$(SOLVER).so
