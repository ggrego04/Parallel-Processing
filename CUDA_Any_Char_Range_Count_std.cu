/********************************************************
 * numberCountCUDA_Unified_solution.c.c
 * Login to machine 103ws??.in.cs.ucy.ac.cy
 * Compile Using:
 *  nvcc  -O3 numberCountCUDA_Unified_solution.cu
 *  ./a.out 100000099 100 3
 *  Number of Arguments 4
 * Array Size: 100000099 SetSize = 100 CountNumber 3
 * Initializing the array with Uniformly Distributed Data...
 * Counting 3 in Uniformly Distributed Data
 * Time: 0.10 Sec. SERIAL: Number if instances found 1000212 (1.000211%).
 * Kernel ONLY:Time: 0.04 Sec.
 * Time: 0.04 Sec. CUDA: Number if instances found 1000212 (1.000211%).
 */
#include <stdio.h>
#include "support.h"
// This Version is configured for GTX750Ti
// Consider one Block per SM
#define ID 999999
#define NUMBER_BLOCKS 10
#define NUMBER_THREADS_PER_BLOCK 32*32
// The size of the Information sent to threads
#define INFO_SIZE 4

void initData(char * vector, int size) {
	int i;
	srand (ID);
	for (i = 0; i < size; i++){
		vector[i] = rand()%26+'a';
		//vector[i] = 'a';
	}
}
void verify(unsigned int  * counters_CPU,unsigned int  * counters_AVX2){
	int j;
	for(j = 0; j < 26; ++j) 
		if (counters_CPU[j]!=counters_AVX2[j]) 
			printf("Does not Verify at value %d\n",j);

}

/* Count the instances of countNumber */
void countNumb_Serial(int size, unsigned int * counters, char * theArray, char from, char to){
  int i;
  for (i = 0; i < size; i++){
	 if (theArray[i]>=from && theArray[i]<=to)                          
	    counters[theArray[i]-'a']++;
  }
}

/***********************************
 * CUDA Implementation
 **********************************/
__global__ void numCountKernel(char *dArray, unsigned int *dCounters, unsigned int * dInfo){
	int i;
	int arraySize = dInfo[0];
	int totalNumberOfThreads =dInfo[1];
    char from = dInfo[2];
    char to = dInfo[3];
	int blockID = blockIdx.x;
	//int threadID = threadIdx.x;
	// Calculate the Unique Thread ID
	int threadUniqueID = threadIdx.x + blockIdx.x * blockDim.x;
    // Each WARP check contigues values in the global memory
	for (i=threadUniqueID;i<arraySize;i+=totalNumberOfThreads)
		if (dArray[i]>=from && dArray[i]<=to) 
            atomicAdd(&dCounters[dArray[i]-'a'+26*blockID],1);
}

void countNumber_CUDA(int arraySize, unsigned int * counters, char * theArray, char from, char to, int blocks, int threads){
	  int totalNumberOfThreads = blocks*threads;
	  int i=0;
	  unsigned int *Counters_CUDA;
      unsigned int *Info;
	  // We will need one Counter per BLOCK
	  int CountersInBytes = blocks*sizeof(unsigned int);
	   // Set the information send to threads
	   // dInfo: [0] the size of the dArray, [1] the countNumber, [1] the Total # of threads
	   // Allocate memory on the CPU for hInfo
	   // *** Using Unified Memory for CUDA ***
	   cudaMallocManaged(&Info, sizeof(unsigned int)*INFO_SIZE);
	   Info[0] = arraySize;
	   Info[1] = totalNumberOfThreads;
       Info[2] = from;
       Info[3] = to;
	
	// Allocate Memory on the CPU for the Counters
	cudaMallocManaged(&Counters_CUDA, sizeof(unsigned int)*CountersInBytes);
	
	// kernel invocation code
	dim3 dimBlock(threads);
	dim3 dimGrid(blocks);
	//Total Time: ?? Sec. Vs Kernel Only Time ?? Sec.
	//Computation ?? Sec Data Transfer ??Sec.
	startTime(1);
	numCountKernel<<<dimGrid, dimBlock>>>(theArray,Counters_CUDA,Info);
	//numCountKernel_shareVar<<<dimGrid, dimBlock>>>(dArray,dCounters,dInfo);
	cudaThreadSynchronize();
	stopTime(1);printf("Kernel ONLY:");elapsedTime(1);printf("\n");
	cudaError_t  code = cudaGetLastError();
	if (code != cudaSuccess){ 
		printf("6) **%s** in %s at line %d\n", cudaGetErrorString(code), __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	int j;
	// Reduction on Host
	for(i=0;i<blocks;i++){
		for(j=0;j<26;j++)
		counters[j] += Counters_CUDA[j+26*i];
	}
    // free the memory allocated on the GPU
    cudaFree(Counters_CUDA);
	cudaFree(Info);
	
}

int main(int argc, char **argv)
{
	int i;
    int blocks=NUMBER_BLOCKS;
    int threads=NUMBER_THREADS_PER_BLOCK;
	int  ArraySize=0, countNumber=0;
    char from, to;

	if (argc < 4){
		printf("Use: numberCount.out <ArraySize> <from> <to>\n");
		printf("NOTE: if numberOfChunks/chankSize > 100 the the number represents chankSize else is numberOfChunks\n");
	}else{
		ArraySize = atoi(argv[1]);
		from = argv[2][0];
		to = argv[3][0];
	}
	//printf("Number of Arguments %d\n",argc);
	//printf("Array Size: %d ArraySize = %d CountNumber %d\n", ArraySize, setSize,countNumber);
     __attribute__ ((aligned (256))) char * theArray = (char *) malloc(sizeof(char) * ArraySize);

    unsigned int  * counters_CPU = (unsigned int *) malloc(sizeof(unsigned int) * 26);
	unsigned int  * counters_CUDA = (unsigned int *) malloc(sizeof(unsigned int) * 26); 
    
	// *** Using Unified Memory for CUDA ***
	cudaMallocManaged(&theArray, sizeof(int)*ArraySize);
  /**********************************************************************/
    initData(theArray, ArraySize);
  /**********************************************************************/
	#ifdef DEBUG
		printf("Uniformly Distributed Data:\n");
		printArray(theArray,ArraySize);
	#endif
	printf("Counting %d in Uniformly Distributed Data\n",countNumber);
    
	startTime(0);
	countNumb_Serial(ArraySize, counters_CPU, theArray, from, to);
	stopTime(0);for (i=0; i<26; i++) printf("%c Found %d times.\n", 'a'+i,counters_CPU[i]);
  	printf("Time CPU:"); elapsedTime(0);
    
	startTime(0);
	countNumber_CUDA(ArraySize, counters_CUDA,theArray, from, to, blocks, threads);
	stopTime(0);for (i=0; i<26; i++) printf("%c Found %d times.\n", 'a'+i,counters_CUDA[i]);
  	printf("Time CUDA:"); elapsedTime(0);

	verify(counters_CPU, counters_CUDA);
  /**********************************************************************/
   cudaFree(theArray);
  return 0;
} 	

