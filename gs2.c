#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/*** Skeleton for Lab 1 ***/

/***** Globals ******/
float **a; /* The coefficients */
float *x;  /* The unknowns */
float *b;  /* The constants */
float err; /* The absolute relative error */
int num = 0;  /* number of unknowns */

float *temp;

/****** Function declarations */
void check_matrix(); /* Check whether the matrix will converge */
void get_input();  /* Read input from file */

/********************************/



/* Function definitions: functions are ordered alphabetically ****/
/*****************************************************************/

/* 
   Conditions for convergence (diagonal dominance):
   1. diagonal element >= sum of all other elements of the row
   2. At least one diagonal element > sum of all other elements of the row
 */
void check_matrix()
{
  int bigger = 0; /* Set to 1 if at least one diag element > sum  */
  int i, j;
  float sum = 0;
  float aii = 0;
  
  for(i = 0; i < num; i++)
  {
    sum = 0;
    aii = fabs(a[i][i]);
    
    for(j = 0; j < num; j++)
       if( j != i)
	 sum += fabs(a[i][j]);
       
    if( aii < sum)
    {
      printf("The matrix will not converge.\n");
      exit(1);
    }
    
    if(aii > sum)
      bigger++;
    
  }
  
  if( !bigger )
  {
     printf("The matrix will not converge\n");
     exit(1);
  }
}


/******************************************************/
/* Read input from file */
/* After this function returns:
 * a[][] will be filled with coefficients and you can access them using a[i][j] for element (i,j)
 * x[] will contain the initial values of x
 * b[] will contain the constants (i.e. the right-hand-side of the equations
 * num will have number of variables
 * err will have the absolute error that you need to reach
 */
void get_input(char filename[])
{
  FILE * fp;
  int i,j;  
 
  fp = fopen(filename, "r");
  if(!fp)
  {
    printf("Cannot open file %s\n", filename);
    exit(1);
  }

 fscanf(fp,"%d ",&num);
 fscanf(fp,"%f ",&err);

 /* Now, time to allocate the matrices and vectors */
 a = (float**)malloc(num * sizeof(float*));
 if( !a)
  {
	printf("Cannot allocate a!\n");
	exit(1);
  }

 for(i = 0; i < num; i++) 
  {
    a[i] = (float *)malloc(num * sizeof(float)); 
    if( !a[i])
  	{
		printf("Cannot allocate a[%d]!\n",i);
		exit(1);
  	}
  }
 
 x = (float *) malloc(num * sizeof(float));
 if( !x)
  {
	printf("Cannot allocate x!\n");
	exit(1);
  }


 b = (float *) malloc(num * sizeof(float));
 if( !b)
  {
	printf("Cannot allocate b!\n");
	exit(1);
  }

 /* Now .. Filling the blanks */ 

 /* The initial values of Xs */
 for(i = 0; i < num; i++)
	fscanf(fp,"%f ", &x[i]);
 
 for(i = 0; i < num; i++)
 {
   for(j = 0; j < num; j++)
     fscanf(fp,"%f ",&a[i][j]);
   
   /* reading the b element */
   fscanf(fp,"%f ",&b[i]);
 }
 
 fclose(fp); 

}


/* 
 * Calculates and updates new unknown values for the current process 
 * 
 * Returns true if any
 */
int calc_unknowns(int my_rank, int comm_sz)
{
  int high_err = 0;
  for (int j = my_rank * num/comm_sz; j < (my_rank+1) * num/comm_sz; j++) 
  {
    float new = b[j];
    for (int p = 0; p < num; p++) 
    {
      if (p != j) 
        new -= x[p] * a[j][p];
    }
    new = new / a[j][j];
    if (fabs((new - x[j]) / new) > err)
      high_err = 1;
    temp[j - (my_rank * num/comm_sz)] = new;
  }
}

int fill_in(int rank, int comm_sz) {
  /* Fill in new unknowns */
  int first_i = rank * num/comm_sz
  for (int j = first_i; j < (rank+1) * num/comm_sz; j++)
  {
    if (fabs((temp[j - first_i] - x[j]) / temp[j - first_i]) > err)
    x[j] = temp[j - first_i];
  }

  return high_err;
}


/************************************************************/


int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  int comm_sz;
  int my_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  //float new = 0;
  int i;
  int temp_err; /* Does this process have an error margin that's too high  */ 

	temp = (float *) malloc((num / comm_sz) * sizeof(float));
	if( !temp)
	{
		printf("Cannot allocate temp!\n");
		exit(1);
	}


  /* Check number of arguments */
  if( argc != 2)
  {
    printf("Usage: gsref filename\n");
    exit(1);
  }

  /* Read the input file and fill the global data structure above */ 
  get_input(argv[1]);

  /* Check for convergence condition */
  /* This function will exit the program if the coffeicient will never converge to 
  * the needed absolute error. 
  * This is not expected to happen for this programming assignment.
  */
  check_matrix();


  MPI_Barrier(MPI_COMM_WORLD);

  if (my_rank == 0) /* Master Process */
  {
    int high_err = 1; /* Is the margin of error too high */
    int nit = 0; /* # of iterations */

    while (high_err) 
    {
      nit++;

      /* Broadcast current unknowns */
      MPI_Bcast(x, num, MPI_FLOAT, 0, MPI_COMM_WORLD);

      /* Calculate this processes unknowns */
      high_err = calc_unknowns(my_rank, comm_sz);
      
      /* Receive/update unkowns and check for completion */
      for (i = 1; i < comm_sz; i++) 
      {
        //MPI_Recv(&x[i * num / comm_sz], num / comm_sz, MPI_FLOAT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(temp, num / comm_sz, MPI_FLOAT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //MPI_Recv(&temp_err, 1, MPI_INT, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (!high_err)
          high_err = temp_err;
      }

      /* Let other processes know if we're done */
      MPI_Bcast(&high_err, 1, MPI_INT, 0, MPI_COMM_WORLD);      

    }

    /* Writing to the stdout */
    /* Keep that same format */
    for( i = 0; i < num; i++)
      printf("%f\n", x[i]);

    printf("total number of iterations: %d\n", nit);

  } else /* Other Processes */
  {
    int keep_going = 1;

    while (keep_going) 
    {
      /* Receive current unkowns */
      MPI_Bcast(x, num, MPI_FLOAT, 0, MPI_COMM_WORLD);

      /* Calculate this processes unknowns */
      temp_err = calc_unknowns(my_rank, comm_sz);

      /* Send out new unknown */
      MPI_Send(&x[my_rank * num / comm_sz], num / comm_sz, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
      MPI_Send(&temp_err, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);


      /* See if we're done */      
      MPI_Bcast(&keep_going, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
  }

  // printf("Process %d:", my_rank);
  // for (i = 0; i < test_sz; i++) {
  // 	printf(" %f,", test[i]);
  // }
  // printf("\n");

  // do {
  //   done = 1;

  //   // Get new for assigned unknowns
  //   // my unknowns => (my_rank * num/comm_sz) -> ((my_rank+1) * num/comm_sz - 1)


  //   if (my_rank != 0) {
  //     // Get the updated unknowns

  //     // Send new to master
  //     for (int j = my_rank * num/comm_sz; j < (my_rank+1) * num/comm_sz; j++) {
  //       x[j] = b[j];
  //       for (int p = 0; p < num; p++) {
  //         if (p != j) {
  //           new -= x[p] * a[j][p];
  //         }
  //       }
  //     }
  //     // Recieve whether or not to continue from master

  //   } else {
  //     int nit = 0;  number of iterations 
  //     // Send the updated unknowns
  //     for (i = 0; i < num; i++) {
  //       MPI_Send(x, )
  //     }

  //     x[0] = 25;
  //     printf("%s\n", "hi there");

  //     nit++;
  //   }
  //   // for (int j = 0; j < num; j++) {
  //   //   old = x[j];

  //   //   // Now calculate new
  //   //   new = b[j];
  //   //   for (int p = 0; p < num; p++) {
  //   //     if (p != j) {
  //   //       new -= x[p] * a[j][p];
  //   //     }
  //   //   }
  //   //   new = new / a[j][j];

  //   //   x[j] = new;

  //   //   if (fabs((new - old) / new) > err)
  //   //     done = 0; 
  //   // }
  // } while (done == 0);

  MPI_Finalize();

  exit(0);
}
