CC = gcc
CFLAGS = -Wall -pedantic -Wshadow -Wextra -g
RM = rm -r
SRC = main
# Assume that SRC (the source code name) is given from shell,
# otherwise assume the name to be "main.c"

$(SRC) : $(SRC).o
	$(CC) $(SRC).o ../../mylib/obj/*.o -lm -fopenmp g.o -o $(SRC)
	rm -r $(SRC).o

$(SRC).o : $(SRC).c
	$(CC) -c $(CFLAGS) -I../../mylib/include/ -fopenmp $(SRC).c

clean :
	$(RM) $(SRC)

