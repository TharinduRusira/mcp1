# Makefile
# For Q4 in Assignment 1 (CS5270 MCP, 2014 S2)
# MSc in CS, Dept of CSE, UoM
# by Sanath Jayasena, June 2014
# Modified from Prof. Kathy Yelick's http://www.cs.berkeley.edu/~yelick/cs194f07/
#
# Include platform-dependent settings.
#
include Makefile.include
#
LDFLAGS += -Ldcmt0.4/lib -ldcmt

Q4_INCS = black_scholes.h gaussian.h parser.h random.h timer.h util.h
Q4_C_SRCS = black_scholes.c gaussian.c main.c parser.c random.c dcmt0.4/lib/random_seed.c timer.c util.c
Q4_C_OBJS = $(Q4_C_SRCS:.c=.o)
Q4_EXE = Q4.x


all: Q4.x

%.o: %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

Q4.x: $(Q4_C_OBJS) dcmt0.4/lib/libdcmt.a
	$(LINKER) $(CFLAGS) $(Q4_C_OBJS) -o $@ $(LDFLAGS)

dcmt0.4/lib/libdcmt.a: 
	make -C dcmt0.4/lib

black_scholes.o: black_scholes.c black_scholes.h gaussian.h random.h util.h

gaussian.o: gaussian.c gaussian.h util.h

main.o: main.c black_scholes.h parser.h random.h timer.h

parser.o: parser.c parser.h

random.o: random.c random.h

dcmt0.4/lib/random_seed.o: dcmt0.4/lib/random_seed.c

timer.o: timer.c timer.h

util.o: util.c util.h

clean:
	make -C dcmt0.4/lib clean
	rm -f $(Q4_C_OBJS) $(Q4_EXE)
