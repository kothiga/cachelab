/*
 * CPSC 4210
 *  - High Performance Parallel Computing
 *
 *    Name: Austin Kothig
 *      ID: 001182645
 *     Sem: Spring 2018
 *
 * Purpose: Simulate the cache management using
 *          Least Resently Used policy
 *
 */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "unistd.h"
#include "getopt.h"
#include "cachelab.h"

#define debug 0

struct cache_line {
   int valid;
   unsigned long int tag;
   int LRU;   
};

/* Function Prototypes */
void help();

int main(int argc, char** argv) {

   int   v =   0 ; //-- verbose flag. 1 if stated
   int   s =  -1 ; //-- # of set bits (read in from command line)
   int   E =  -1 ; //-- associativity 
   int   b =  -1 ; //-- # of byte offset bits
   char *t = "-1"; //-- Trace File 
   
   unsigned long int set_index_mask;

   //-- loop through arguments
   int opt;
   while ((opt = getopt(argc, argv, "hvs:E:b:t:")) != -1) {
      switch (opt) {
	 case 'h': help(); return 0; break;
	 case 'v': v = 1;            break;
	 case 's': s = atoi(optarg); break;
	 case 'E': E = atoi(optarg); break;
	 case 'b': b = atoi(optarg); break;
	 case 't': t = optarg;       break;
	 default :
	    printf("wrong argument\n");
	    return 0; break;
      }
   }

   //-- check to see if we missed any arguments
   if (s == -1 || E == -1 || b == -1 ||
       (strcmp("-1", t) == 0) ) {
      printf("./csim: Missing required command line argument\n");
      help();
      return 0;
   }

   //-- used for generating tags from addresses
   set_index_mask = pow(2,s)-1;

   //-- debug inputs, print them out to see
   #if debug
   printf("v is %d\n", v);
   printf("s is %d\n", s);
   printf("E is %d\n", E);
   printf("b is %d\n", b);
   printf("t is %s\n", t);
   printf("Set Mask is %lu\n", set_index_mask);
   #endif

   
   //-- allocate the cache
   struct cache_line **cache;
   cache = (struct cache_line **)malloc(sizeof(struct cache_line *)*(set_index_mask+1));
   for (int i = 0; i < set_index_mask+1; i++) {
      cache[i] = (struct cache_line *)malloc(sizeof(struct cache_line) * E);
   }

   //-- make sure all values in the cache
   //-- are initialized to 0
   for (int i = 0; i < set_index_mask+1; i++) {
      for (int j = 0; j < E; j++) {
	 cache[i][j].tag = 0;
	 cache[i][j].valid = 0;
	 cache[i][j].LRU = 0;
      }
   }
      
   //-- initialize reading in
   FILE *pFile;

   //-- open the trace file
   pFile = fopen(t,"r");
   char identifier;
   unsigned long int address;
   int size;

   
   //-- start simulating the cache
   int hit = 0, miss = 0, evict = 0, ts = 0;
   while (fscanf(pFile," %c %lx,%d", &identifier, &address, &size) > 0) {

      //-- sim should ignore all instruction cache accesses
      if (identifier == 'I') {
	 continue;
      }

      //-- extract set index from address by bit shifting
      //-- 'b' bits and 'anding' with our mask, then get
      //-- the tag bit from the passed in address
      unsigned long int index, tag;
      index = ((address >> b) & set_index_mask);
      tag = (address >> (b + s));
      
      char* verbose;
      int found = 0, evic = 0, val = 0;

      //-- go through each cache block
      //-- and look in the extracted index
      for (int i = 0; i < E; i++) {

	 //-- look for a valid entry
	 if (cache[index][i].valid == 1) {

	    //-- compare the tags
	    if (cache[index][i].tag == tag) {
	       
	       //-- say that we have used this
	       cache[index][i].LRU = ts++;
	       found = 1;
	       break;
	    }
	    
	 } else {
	    
	    //-- if we found a non-valid, we can insert here
	    val = 1;
	    cache[index][i].valid = 1;
	    cache[index][i].tag = tag;
	    cache[index][i].LRU = ts++;
	    break;
	 }
      }

      //-- if the data was not found, and we did not insert 
      //-- into an  empty valid bit, look for the smallest LRU 
      //-- and replace it in the cache
      if (found == 0 && val == 0) {
	 int j = 0;
	 for (int i = 1; i < E; i++) {
	    if (cache[index][j].LRU > cache[index][i].LRU) {
	       j = i;
	    }
	 }
	 evic = 1;
	 cache[index][j].valid = 1;
	 cache[index][j].tag = tag;
	 cache[index][j].LRU = ts++;	 
      }

      //-- based on the passed in identifier update the
      //-- hit, miss, and evict variables as necessary.
      //-- Set verbose too incase user specified -v
      switch (identifier) {

	 //-- Load Instruction
	 case 'L':
	    if (found == 1) {       //-- found=1 evic=X
	       verbose = "hit";
	       hit++;
	    } else if (evic == 1) { //-- found=0 evic=1
	       verbose = "miss eviction";
	       miss++; evict++;
	    } else {                //-- found=0 evic=0
	       verbose = "miss";
	       miss++;
	    }
	    break;

         //-- Store Instruction
	 case 'S':
	    if (found == 1) {       //-- found=1 evic=X
	       verbose = "hit";
	       hit++;
	    } else if (evic == 1) { //-- found=0 evic=1
	       verbose = "miss eviction";
	       miss++; evict++;
	    } else {                //-- found=0 evic=0
	       verbose = "miss";
	       miss++;
	    }
	    break;

         //-- Data Modify Instruction
	 case 'M':
	    if (found == 1) {       //-- found=1 evic=X
	       verbose = "hit hit";
	       hit++; hit++;
	    } else if (evic == 1) { //-- found=0 evic=1
	       verbose = "miss eviction hit";
	       miss++; evict++; hit++;
	    } else {                //-- found=0 evic=0
	       verbose = "miss hit";
	       miss++; hit++;
	    }
	    break;
      }
      
      if (v) { //-- print out this line if using verbose
	 printf("%c %lx,%d %s\n", identifier, address, size, verbose);
      }
   }
   
   //-- print the results
   printSummary(hit, miss, evict);

   //-- deallocate file
   fclose(pFile);

   //-- deallocate the array
   for (int i = 0; i < set_index_mask+1; i++) {
      free(cache[i]);
   }
   free(cache);
   
   return 0;
}


void help() {
   printf("Usage: ./csim [-hv] -s <num> -E <num> -b <num> -t <file>\n");
   printf("Options:\n");
   printf("  -h\t\tPrint this help message.\n");
   printf("  -v\t\tOptional verbose flag.\n");
   printf("  -s <num>\tNumber of set index bits.\n");
   printf("  -E <num>\tNumber of lines per set.\n");
   printf("  -b <num>\tNumber of block offset bits.\n");
   printf("  -t <file>\tTrace file.\n\n");
   printf("Examples:\n");
   printf("linux> ./csim -s 4 -E 1 -b 4 -t traces/yi.trace\n");
   printf("linux> ./csim -v -s 8 -E 2 -b 4 -t traces/yi.trace\n");
}
