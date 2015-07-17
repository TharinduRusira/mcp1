/*
 *	Author: Tharindu Rusira, Department of Computer Science and Engineering, University of Moratuwa
 *	Program: Monte Carlo PI calculation serial/parallel implementations
 */

#include <iostream>
#include <math.h>
#include <time.h>
#include <random>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include "omp.h"

#define ITR 1000*1000*1000

#define GET_TIME(x);	if (clock_gettime(CLOCK_MONOTONIC, &(x)) < 0)\
			{ perror("clock_gettime():"); exit(EXIT_FAILURE);} 

// get elapsed time in ms
float elapsed_time_msec(struct timespec *begin, struct timespec *end, long *sec,long *nsec)
{
  if (end->tv_nsec < begin->tv_nsec) 
  {
    *nsec = 1000000000 - (begin->tv_nsec - end->tv_nsec);
    *sec = end->tv_sec - begin->tv_sec -1;
  }
  else 
  {
    *nsec = end->tv_nsec - begin->tv_nsec;
    *sec = end->tv_sec - begin->tv_sec;
  }
  return (float) (*sec) * 1000 + ((float) (*nsec)) / 1000000;
}

const double xy_max = 1;

double pi;
int in_circle;
int n_threads;
pthread_mutex_t mutex;

using namespace std;

bool in_or_out(double x, double y,const double xy_max){
  double center = xy_max/2;
  if(pow((x - center),2)+pow((y - center),2)<= pow(center,2))
    return true;
  return false;
}

float pi_sequential(){
  pi = 0;
  in_circle = 0;
  double x,y;
  
  //time measurement
  struct timespec t1, t2;
  long sec, nsec;
  float time_elapsed; 
  
  srand(time(NULL));
  
  //default_random_engine generator;
  //uniform_real_distribution<double> distribution(0,xy_max);
  //time_t begin = time(NULL);
  
  GET_TIME(t1);
  for(int i=0;i<ITR; i++){
    x = (double)random()/(double)RAND_MAX;//distribution(generator);
    y = (double)random()/(double)RAND_MAX;//distribution(generator);
    if(in_or_out(x,y,xy_max))
      in_circle++;
  }
  pi = 4* ((double)in_circle/(double)ITR);
  GET_TIME(t2);
  
  time_elapsed = elapsed_time_msec(&t1,&t2,&sec,&nsec);
  cout << "Pi (seq) = " << pi << endl;
  cout << "Time (seq) = " << time_elapsed << endl;
  return time_elapsed;
}
///////////////////////////////////////////////////
/*
 * OpenMP parallel version
 * The loop in the above sequential code is massively (embarassingly) parallel. We use this behavor to parallelize the computation
 */

/*
double pi_parallel_omp(const double xy_max){
  double pi = 0 ;
  int in_circle = 0;
  
  default_random_engine generator;
  uniform_real_distribution<double> distribution(0,1);
  
  time_t begin = time(NULL);
  
  omp_set_num_threads(2);
  double x,y;
  #pragma omp parallel for private(x,y) shared(in_circle)
  for(int i=0;i<ITR; i++){
    x = distribution(generator);
    y = distribution(generator);
    if(in_or_out(x,y,xy_max))
      in_circle++; 
  }
  pi = 4* ((double)in_circle/(double)ITR);
  double time_elapsed = time(NULL) - begin;
  cout << "Time (par) = " << time_elapsed << endl;
  return pi;
}
*/
////////////////////////////////////////////////

/*
 * Pthread based implementations
 */
void* threaded_pi(void* arg)
{
 
 double x,y;
 unsigned int* seed = (unsigned int*)&arg; // thread-local seed for the re-entrant random number generator
 int itrs_per_thread = ITR/n_threads;
 
 int my_in_circle=0;
 
 for(int k=0;k<itrs_per_thread;k++)
 {
   x = (double)rand_r(seed)/(double)RAND_MAX;
   y = (double)rand_r(seed)/(double)RAND_MAX;
      
   if(in_or_out(x,y,xy_max))
   {
    my_in_circle++ ;
   }
  } 
  
  pthread_mutex_lock(&mutex);
  in_circle+=my_in_circle;
  pthread_mutex_unlock(&mutex);
  pthread_exit(0);
}

float pi_parallel_pthreads()
{
  pi =0;
  in_circle = 0; // reset shared variable
  pthread_t threads[n_threads];  
  
  //time measurement
  struct timespec t1, t2;
  long sec, nsec;
  float total_pthread_time; 
  
  GET_TIME(t1);
  for(int i=0; i<n_threads; i++)
  {
    // create new thread
    pthread_create(&threads[i],NULL,threaded_pi,(void*)i);
  }
  // collect all results and finalize
  for(int i=0; i<n_threads; i++)
  {
    // create new thread
    pthread_join(threads[i],NULL); 
  }
  pi = 4 * ((double)in_circle/(double)ITR);
  GET_TIME(t2);
  
  total_pthread_time = elapsed_time_msec(&t1,&t2,&sec,&nsec);
  cout << "Pi (pthreads) = " << pi << endl;
  cout << "Pthread time = " << total_pthread_time <<endl;
  // exit the main thread
  //pthread_exit(NULL);
  return total_pthread_time;
}

///////////////////////////////////////////////

int main(int argc, char *args[]) {
  
    int itrs = 30;
     // do argument valildation
    if(argc == 2 && strcmp(args[1],"-s")==0)
    {
      double seq_time =0;
     // sequential mode activated 
     // run 30 iterations and get average time
      for(int i=0;i<itrs;i++)
	seq_time+=pi_sequential();
      
      cout << "Avg. seq. time = " << seq_time/(double)itrs<<endl;
    }
    else if(argc == 3 && strcmp(args[1],"-p")==0)
    {
      double par_time = 0;
      // parallel mode activated. Check for the number of threads now
      // run 30 iterations and get average time
      n_threads = atoi(args[2]);
      for(int i=0;i<itrs;i++)
	par_time += pi_parallel_pthreads();
      
      pthread_mutex_destroy(&mutex);
      cout << "Avg. par. time = " << par_time/(double)itrs<<endl;
    }else 
    {
      cout << "Please provide valid input arguments" << endl;
    }
    return 0;
}
