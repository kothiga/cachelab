/*
 * CPSC 4210
 *  - High Performance Parallel Computing
 *
 *    Name: Austin Kothig
 *      ID: 001182645
 *     Sem: Spring 2018
 *
 * Purpose: Transpose a Matrix using good 
 *          temporal locality
 *
 */

/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"

/* Prototype for user Made functions */
void transpose_64(int M, int N, int A[N][M], int B[M][N]);

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N]) {

   //-- By default Block Size is set to 16
   int i, j, k, l, b = 16;
   int temp, d;

   //-- Set Block Size for 32x32 is 8
   if (N == 32 && M == 32) { b = 8; }
  
   //-- For a 64x64 transpose use a different approach
   if (N == 64 && M == 64) { transpose_64(M, N, A, B); }

   //-- For 32x32, and 61x67 transpose, this approach is okay
   //-- if 32x32 -- Block Size = 8
   //-- if 61x67 -- Block Size = 16
   else {
      //-- Loop through the array in blocks
      for (i = 0; i < N; i += b) {
	 for (j = 0; j < M; j += b) {

	    //-- In the current block, make sure we're in bounds
	    for (k = 0; k < b && i + k < N; k++) {
	       for (l = 0; l < b && j + l < M; l++) {

		  //-- Make sure we are not looking
		  //-- at a diagonal element
		  if (j+l != i+k) {
		     
		     //-- Put this element in the transposed matrix
		     B[j+l][i+k] = A[i+k][j+l];
		  
		  } else {

		     //-- Store the diagonal element
		     temp = A[i+k][j+l];
		     d = i+k; 
		  }
	       }

	       //-- The diagonal is in mem so store
	       if (i == j) {
		  B[d][d] = temp;
	       }
	    }
	 }
      }
   }
}


/*
 * transpose_64
 *
 * Special way of doing the 64x64 transposition.
 * Break bigger blocks (size 8x8) into smaller blocks
 * (Four 4x4 blocks). Use a temp array of size 8,
 * as a buffer for moving blocks of data
 *
 */
void transpose_64(int M, int N, int A[N][M], int B[M][N]) {

   //-- our most outter block size is of size 8
   //-- inner blocks are of size 4
   int i, j, k, l;
   //int temp[8];
   int temp0, temp1, temp2, temp3;
   int temp4, temp5, temp6, temp7;

   //-- Loop throgh the Big Blocks
   for (i = 0; i < N; i += 8) {
      for (j = 0; j < M; j += 8) {

	 //-- construct top two mini blocks into B
	 for (k = i; k < i+4; k++) {

	    //-- get the next 8 var and
	    //-- store them in temporary memory
	    temp0 = A[k][j+0];
	    temp1 = A[k][j+1];
	    temp2 = A[k][j+2];
	    temp3 = A[k][j+3];
	    temp4 = A[k][j+4];
	    temp5 = A[k][j+5];
	    temp6 = A[k][j+6];
	    temp7 = A[k][j+7];


	    //-- store the first half of
	    //-- local reg into B
	    B[j+0][k] = temp0;
	    B[j+1][k] = temp1;
	    B[j+2][k] = temp2;
	    B[j+3][k] = temp3;

	    
	    //-- store the second half of
	    //-- local reg into B
	    B[j+0][k+4] = temp4;
	    B[j+1][k+4] = temp5;
	    B[j+2][k+4] = temp6;
	    B[j+3][k+4] = temp7;
	 }

	 
	 //-- construct bottom two mini blocks into B
	 for (k = j; k < j+4; k++) {

	    //-- get the next 8 var and
	    //-- store them in local reg
	    temp4 = A[i+4][k];
	    temp5 = A[i+5][k];
	    temp6 = A[i+6][k];
	    temp7 = A[i+7][k];
	    

	    //-- retrieve the temporarily stored
	    //-- chunk from before, load it into cache
	    temp0 = B[k][i+4];
	    temp1 = B[k][i+5];
	    temp2 = B[k][i+6];
	    temp3 = B[k][i+7];
	    
	    
	    //-- now that this block is back in cache
	    //-- fill this block with the values it should have
	    B[k][i+4] = temp4;
	    B[k][i+5] = temp5;
	    B[k][i+6] = temp6;
	    B[k][i+7] = temp7;

	    
	    //-- store the next section
	    B[k+4][i+0] = temp0;
	    B[k+4][i+1] = temp1;
	    B[k+4][i+2] = temp2;
	    B[k+4][i+3] = temp3;

	    
	    //-- store the last block straight from A
	    for (l = 0; l < 4; l++) {
	       B[k+4][i+l+4] = A[i+l+4][k+4];
	    }
	 }
      }
   }
}





/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
   int i, j, tmp;

   for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
	 tmp = A[i][j];
	 B[j][i] = tmp;
      }
   }    

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
   /* Register your solution function */
   registerTransFunction(transpose_submit, transpose_submit_desc); 

   /* Register any additional transpose functions */
   registerTransFunction(trans, trans_desc); 

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
   int i, j;

   for (i = 0; i < N; i++) {
      for (j = 0; j < M; ++j) {
	 if (A[i][j] != B[j][i]) {
	    return 0;
	 }
      }
   }
   return 1;
}

