/*
 *	Author: Tharindu Rusira, Department of Computer Science and Engineering, University of Moratuwa
 *	Program: Vector dot product serial/parallel implementations
 * 	Time measurement code is extracted from prof. Sanath Jayasena's class material.
 */

#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <errno.h>
#include <stdio.h>

#ifndef SIZE
#define SIZE 1000*1000*1000
#endif

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

int n_threads;
pthread_mutex_t mutex;

using namespace std;

// this struct will be passed to each thread
struct thread_arg{
  double* A;
  double* B;
  double* global_sum;
  int thread_id;
};

double* generate_random_vector(int size)
{
  srand(time(NULL));
  double* X = (double*)malloc(size*sizeof(double));
  for(int i=0; i<size;i++)
  {
    X[i] = (double)rand();
  }
  return X;
}

float sequential_vector_dot_product(double* A, double* B)
{
    double dot_product = 0;
    
    struct timespec t1, t2;
    long sec, nsec;
    float time_elapsed; 
    
    GET_TIME(t1);
    for(int i =0;i<SIZE; i++)
    {
      dot_product+= A[i]*B[i];
    }
    GET_TIME(t2);
    time_elapsed = elapsed_time_msec(&t1,&t2,&sec,&nsec);
    cout << "Dot Product (seq) = " << dot_product << endl;
    cout << "Time (seq) = " << time_elapsed << endl;
    
    return time_elapsed;
}

void* threaded_dot_product(void* arg)
{
  // all input sizes are divisible by n_threads (number of threads)
  int elements_per_thread = SIZE/n_threads;
  
  thread_arg my_args = *(thread_arg*)arg;
  
  int start  = my_args.thread_id*elements_per_thread;
  int end = start + elements_per_thread;
  double partial_sum =0;
  
  for(int i = start; i<end;i++)
  {
    partial_sum+=my_args.A[i]*my_args.B[i];
  }
  
  pthread_mutex_lock(&mutex);
  *(my_args.global_sum) += partial_sum;
  pthread_mutex_unlock(&mutex);
  
  pthread_exit(NULL);
 
}

float parallel_vector_dot_product(double* A, double* B)
{
  double parallel_dot_product = 0;
  pthread_t threads[n_threads];  
  
  struct timespec t1, t2;
  long sec, nsec;
  float total_threaded_time; 
  
  GET_TIME(t1);
  
  for(int i=0; i<n_threads; i++)
  {
    thread_arg thread_input;
    thread_input.A = A;
    thread_input.B = B;
    thread_input.global_sum = &parallel_dot_product;
    thread_input.thread_id = i;
    // create new thread
    pthread_create(&threads[i],NULL,threaded_dot_product,(void*)&thread_input);
  }
  // collect all results and finalize
  for(int i=0; i<n_threads; i++)
  {
    // create new thread
    pthread_join(threads[i],NULL); 
  }
  GET_TIME(t2);
  
  total_threaded_time = elapsed_time_msec(&t1,&t2,&sec,&nsec);
  
  cout << "Dot product (par) = "<< parallel_dot_product << endl;
  cout << "Total time (par) = " << total_threaded_time << endl;
  return total_threaded_time;
}

int main(int argc, char* args[]) {
    int itrs = 30;
    double* A = generate_random_vector(SIZE) ;
    double* B = generate_random_vector(SIZE);
     // do argument valildation
    if(argc == 2 && strcmp(args[1],"-s")==0)
    {
     float seq_time =0;
     
     // sequential mode activated 
     // run 30 iterations and get average time
      for(int i=0;i<itrs;i++)
	seq_time+=sequential_vector_dot_product(A,B);
      
      cout << "Avg. seq. time = " << seq_time/(float)itrs<<endl;
    }
    else if(argc == 3 && strcmp(args[1],"-p")==0)
    {
      float par_time = 0;
      // parallel mode activated. Check for the number of threads now
      // run 30 iterations and get average time
      n_threads = atoi(args[2]);
      for(int i=0;i<itrs;i++)
	par_time += parallel_vector_dot_product(A,B);
      
      pthread_mutex_destroy(&mutex);
      
      cout << "Avg. par. time = " << par_time/(float)itrs<<endl;
    }else 
    {
      cout << "Please provide valid input arguments" << endl;
    }
    
    free(A); free(B);
    return 0;
}
