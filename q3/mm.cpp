#include <iostream>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <errno.h>
#include <stdio.h>

#ifndef SIZE
#define SIZE 1500
#endif

#define GET_TIME(x);	if (clock_gettime(CLOCK_MONOTONIC, &(x)) < 0)\
			{ perror("clock_gettime():"); exit(EXIT_FAILURE);}

int n_threads;

struct thread_input{
  double** A;
  double** B;
  double** C;
  int thread_id;
};

using namespace std;

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

double** generate_random_matix()
{
  double** MATX = (double**)malloc(SIZE*sizeof(double*)); // this is bad practice to pass local pointers, but YOLO !
  for(int i=0; i< SIZE; i++)
  {
      MATX[i] = (double*)malloc(SIZE*sizeof(double)); // assign space for each column
  }
  for(int row=0; row < SIZE; row++)
  {
    for(int col=0; col < SIZE; col++)
    {
      MATX[row][col] = ((double)rand()/(double)RAND_MAX) + 1; // generate random numbers between 1 and 2

    }

  }

  return MATX;
}

float sequential_mm(double** MATA, double** MATB, double** MATC)
{
  float total_seq_time = 0;
  
  struct timespec t1, t2;
  long sec, nsec;
  float time_elapsed;

  GET_TIME(t1);
  // matrix multiplication. Sequential Version
  for(int row_c=0; row_c<SIZE;row_c++)
  {
    for(int col_c=0;col_c<SIZE;col_c++)
    {
      MATC[row_c][col_c]= 0;
      for(int k=0;k<SIZE;k++)
      {
	MATC[row_c][col_c]+=MATA[row_c][k]*MATB[k][col_c]; // Transposing B is a possible optimization, not implemented here.
      }
     
    }

  }
  
  GET_TIME(t2);

  total_seq_time = elapsed_time_msec(&t1,&t2,&sec,&nsec);
  cout << "Mat mul time (seq) = " << total_seq_time << endl;
  return total_seq_time;
}

void* threaded_mm(void* arg)
{
  const int lines_per_thread = SIZE/n_threads;
  const int lines_per_thread_last = lines_per_thread + SIZE%n_threads;
  
  // we identify each threads' region in A and B,
  thread_input my_args = *(thread_input*)arg;
  int start = my_args.thread_id* lines_per_thread; 
  int end;
  
  if(my_args.thread_id != (n_threads-1))
  {
    end = start + lines_per_thread;
  }else
  {
    end = start + lines_per_thread_last;  
  }
    
  // calculate matrix multiplication within the calculated boundary
  for(int i = start; i < end; i++)
  {
    for(int j= 0; j < SIZE;j++)
    {
      my_args.C[i][j] = 0;
      for(int k =0; k<SIZE; k++)
      {
	my_args.C[i][j] += my_args.A[i][k]*my_args.B[k][j];
      }
    }
  }
}

float parallel_mm(double** ParMatA, double** ParMatB, double** ParMatC)
{  
  pthread_t threads[n_threads];
  thread_input t_in[n_threads];
  // time calculation
  struct timespec t1, t2;
  long sec, nsec;
  float total_parallel_time;
  
  GET_TIME(t1);
  
  for(int i=0; i<n_threads; i++)
  {
    t_in[i].A = ParMatA;
    t_in[i].B = ParMatB;
    t_in[i].C = ParMatC;
    t_in[i].thread_id = i;
    // create new thread
    pthread_create(&threads[i],NULL,threaded_mm,(void*)&t_in[i]);
  }
  // collect all results and finalize
  for(int i=0; i<n_threads; i++)
  {
    // create new thread
    pthread_join(threads[i],NULL); 
  }
  GET_TIME(t2);
  
  total_parallel_time = elapsed_time_msec(&t1,&t2,&sec,&nsec);
  cout << "Mat mul time (par) = " << total_parallel_time << endl;
  return total_parallel_time;
}


int main(int argc, char* args[]) {
  
  double** MATA;
  double** MATB;
  double** ParMatA;
  double** ParMatB;
  
  int itrs = 35; // remove first two values 
  
  if(argc == 2 && strcmp(args[1],"-s")==0)
    {
     float seq_time =0;
     
    // sequential mode activated 
    MATA = generate_random_matix();
    MATB = generate_random_matix();
    double** MATC = (double**)malloc(SIZE*sizeof(double*));	
    for(int i=0;i<SIZE;i++)
    {
      MATC[i] = (double*)malloc(SIZE*sizeof(double));
    }
     
     // run 30 iterations and get average time
      for(int i=0;i<itrs;i++)
	{
		if(i<5){sequential_mm(MATA,MATB,MATC);continue;} // warm up
		seq_time+=sequential_mm(MATA,MATB, MATC);
	}	
      cout << "Avg. seq. time = " << seq_time/(float)(itrs-5)<<endl;
      free(MATC);
    }
    else if(argc == 3 && strcmp(args[1],"-p")==0)
    {
      float par_time = 0;
      // parallel mode activated. Check for the number of threads now
      ParMatA = generate_random_matix();
      ParMatB = generate_random_matix();
      
      double** ParMatC = (double**)malloc(SIZE*sizeof(double*));					
      for(int i=0;i<SIZE; i++)
      {
	ParMatC[i] = (double*)malloc(SIZE*sizeof(double));
      }
      // run 30 iterations and get average time
      n_threads = atoi(args[2]);
      for(int i=0;i<itrs;i++)
	{
		if(i<5){parallel_mm(ParMatA,ParMatB,ParMatC);continue;}
		par_time += parallel_mm(ParMatA, ParMatB, ParMatC);
	}
      
      cout << "Avg. par. time = " << par_time/(float)(itrs-5)<<endl;
      free(ParMatC);
    }else if(argc == 2 && strcmp(args[1],"-v")==0)
    {
      // verify mode
      
      MATA = generate_random_matix();
      MATB = generate_random_matix();
      ParMatA = MATA;
      ParMatB = MATB;
      
      double** MATC = (double**)malloc(SIZE*sizeof(double*));	
      double** ParMatC = (double**)malloc(SIZE*sizeof(double*));					
      for(int i=0;i<SIZE;i++)
      {
	MATC[i] = (double*)malloc(SIZE*sizeof(double));
	ParMatC[i] = (double*)malloc(SIZE*sizeof(double));
      }
      
      sequential_mm(MATA,MATB,MATC);
      n_threads = 8; // use this default value
      parallel_mm(ParMatA,ParMatB,ParMatC);
      
      // now we find the maximum difference between two resulting matrices
      float maxdiff = 0.0;
      
      for(int i=0;i< SIZE; i++)
      {
	for(int j=0;j< SIZE; j++)
	{
	  //find the maxdiff of any two corresponding elements
	  if(fabs(MATC[i][j]-ParMatC[i][j])>maxdiff)
	  {
	    maxdiff = fabs(MATC[i][j]-ParMatC[i][j]);
	  }
	}
      }
      cout << "Maximum difference between any two corresponding elements = " << maxdiff << endl;
      free(MATC); 
      free(ParMatC);
    }
    else 
    {
      cout << "Please provide valid input arguments" << endl;
    }
    return 0;
}
