.SUFFIXES: .o .c
CC = gcc 
DEBUG = -g -Wall
#DEBUG = -Wall
FLAGS = $(DEBUG) -O2 -static -I"C:/Program Files/GSL-1.13/include"
LIBS = -L"c:/Program Files/GSL-1.13/lib" -lgsl -lgslcblas -lm -lstdc++

CORE = hcc_dclinde.o parse_option.o dtsv.o globalMIT.o globalMIT_fixedOrder.o mt19937ar_lib.o

CORE2 = grn_cmp_hcc.o parse_option.o

EXEC = hcc_dclinde

EXEC2 = grn_cmp_hcc

all: $(EXEC) $(EXEC2)

globalMIT.o: globalmit_source/globalMIT.cpp
	g++ -O2 -Wall -c -o globalMIT.o globalmit_source/globalMIT.cpp

globalMIT_fixedOrder.o: globalmit_source/globalMIT_fixedOrder.cpp
	g++ -O2 -Wall -c -o globalMIT_fixedOrder.o globalmit_source/globalMIT_fixedOrder.cpp

.c.o: 
	$(CC) $(FLAGS) -c $< -o $@

core : $(CORE)

$(EXEC) : $(CORE) $(OLSLIB)
	$(CC) $(FLAGS) $(CORE) -o $(EXEC) $(LIBS)
	@echo 'Made '

$(EXEC2) : $(CORE2)
	$(CC) $(FLAGS) $(CORE2) -o $(EXEC2) -lm
	@echo 'Made '

clean:
	rm -f *.o *.bak *.*~
	rm -f $(EXEC)
	@echo 'Made'
