CC = gcc
CFLAGS = -Wall -pedantic -Wshadow -Wextra -g
# Add here the appropriate flag for optimization
RM = rm -r
SRC = main
# Assume that SRC (the source code name) is given from shell,
# otherwise assume the name to be "main.c"

$(SRC) : $(SRC).o
	$(CC) $(SRC).o ../mylib/obj/*.o -lm -lpython3.6m -o $(SRC)
	rm -r $(SRC).o

$(SRC).o : $(SRC).c
	$(CC) -c $(CFLAGS) -I../mylib/include/ -I/usr/include/python3.6 $(SRC).c

clean :
	$(RM) $(SRC)

