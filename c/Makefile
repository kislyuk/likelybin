MANIFEST = *.c *.h Makefile
PROJECT = mcmc-core
HOST = topaz
AUTHOR = kislyuk

CC = gcc
OPTFLAGS = -O3
#DEBUGFLAGS = -DDEBUG -pg
XWARNINGS = -Wextra -Wstrict-prototypes -Wno-unused-variable -Wno-unused-parameter
#XWARNINGS := ${XWARNINGS} -Wc++-compat
#XWARNINGS := ${XWARNINGS} -Wconversion -Wlogical-op
#XWARNINGS := ${XWARNINGS} -Wmissing-declarations
CFLAGS = $(OPTFLAGS) $(DEBUGFLAGS) -std=c99 -pedantic -Wall $(XWARNINGS)
#-march=native
#CFLAGS = $(OPTFLAGS) -Wall -pedantic
# -DDEBUG
CLINKER = gcc
LIBDIR = -L/usr/local/lib
MLIB = -lm
INCDIR = -I./
TRGT = ./${PROJECT}
OBJECTS = *.o
#LINKER_OPTS = -static

.c.o:
	$(CC) $(CFLAGS) $(INCDIR) -c $*.c

#$(TRGT): $(OBJECTS)
$(TRGT):
	$(CC) $(CFLAGS) $(INCDIR) -c *.c
	$(CLINKER) $(OPTFLAGS) $(DEBUGFLAGS) $(OBJECTS) $(LIBDIR) -o $(TRGT) $(MLIB) $(LINKER_OPTS)

clean:
	rm -f *.o ./*~ *~ core $(TRGT)

uninstall:
	rm $(TRGT)

up:
	(cd ..; make up)

test:
	time ./${PROJECT} -o 2 -s 10000 -n 2 < test-input

roll-dist:
	(cd ..; make roll-dist)

profile:
#	OPTFLAGS = -pg
#	make clean
#	make
#	make test
	gprof $(TRGT)
