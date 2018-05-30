CC	      = g++
CC_FLAGS        = -g3 -O3 -Wall -I  /usr/local/include/gsl/
LD_FLAGS	= -L/usr/local/lib -lgsl -lgslcblas -lcblas -lm
BAS		= basicmodel.o models.o io.o utilities.o
DC		= dist_check.o io_check.o utilities.o models_check.o
GD		= makegrid.o models.o io.o utilities.o

basic: $(BAS)
	$(CC) $(CC_FLAGS) $(BAS) -o run_basic  $(LD_FLAGS)
basicmodel.o: basicmodel.cpp
	$(CC) $(CC_FLAGS) -c basicmodel.cpp
models.o: models.cpp
	$(CC) $(CC_FLAGS) -c models.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp

grid: $(GD)
	$(CC) $(CC_FLAGS) $(GD) -o run_grid $(LD_FLAGS)
makegrid.o: makegrid.cpp
	$(CC) $(CC_FLAGS) -c makegrid.cpp
