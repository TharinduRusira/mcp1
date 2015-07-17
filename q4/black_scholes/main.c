// main.c
// For Q4 in Assignment 1 (CS5270 MCP, 2014 S2)
// MSc in CS, Dept of CSE, UoM
// by Sanath Jayasena, June 2014
// Modified from Prof. Kathy Yelick's http://www.cs.berkeley.edu/~yelick/cs194f07/

#include "black_scholes.h"
#include "parser.h"
#include "random.h"
#include "timer.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Usage: ./Q4.x <filename> <nthreads>
 *
 * <filename> (don't include the angle brackets) is the name of 
 * a data file in the current directory containing the parameters
 * for the Black-Scholes simulation.  It has exactly six lines 
 * with no white space.  Put each parameter one to a line, with
 * an endline after it.  Here are the parameters:
 *
 * S
 * E
 * r
 * sigma
 * T
 * M
 *
 * <nthreads> (don't include the angle brackets) is the number of
 * worker threads to use at a time in the benchmark.  The sequential
 * code which we supply to you doesn't use this argument; your code
 * will.
 */
int
main (int argc, char* argv[])
{
  confidence_interval_t interval;
  double S, E, r, sigma, T;
  int M = 0;
  char* filename = NULL;
  int nthreads = 1;
  double t1, t2;
  
  if (argc < 3)
    {
      fprintf (stderr, 
	       "Usage: ./Q4.x <filename> <nthreads>\n\n");
      exit (EXIT_FAILURE);
    }
  filename = argv[1];
  nthreads = to_int (argv[2]);
  parse_parameters (&S, &E, &r, &sigma, &T, &M, filename);

  /* 
   * Make sure init_timer() is only called by one thread,
   * before all the other threads run!
   */
  init_timer ();

  /* Same goes for initializing the PRNG */
  init_prng (random_seed ());

  /*
   * Run the benchmark and time it.
   */
  t1 = get_seconds ();
  black_scholes (&interval, S, E, r, sigma, T, M,nthreads);
  t2 = get_seconds ();

  /*
   * A fun fact about C string literals (i.e., strings enclosed in
   * double quotes) is that the C preprocessor automatically
   * concatenates them if they are separated only by whitespace.
   */
  printf ("Black-Scholes benchmark:\n"
	  "------------------------\n"
	  "S        %g\n"
	  "E        %g\n"
	  "r        %g\n"
	  "sigma    %g\n"
	  "T        %g\n"
	  "M        %d\n",
	  S, E, r, sigma, T, M);
  printf ("Confidence interval: (%g, %g)\n", interval.min, interval.max);
  printf ("Total simulation time: %g seconds\n", t2 - t1);

  return 0;
}



