/********************************************************
 * AVX2_Solution.c
 * Login to machine 103ws??.in.cs.ucy.ac.cy 
 * Compile using:
 * gcc -mavx2 -lpthread -O3 AVX2_Any_Char_Range_Count_std.c -o AVX2_Any_Char_Range_Count_std.out
 * ./AVX2_Any_Char_Range_Count_std.out  1000000000 d t threads                 
 *  ==========================================
 ***************************************************************************
 * Author: Petros Panayi, PhD, CS, UCY, Jan. 2020
 ***************************************************************************
 * Reference: https://software.intel.com/sites/landingpage/IntrinsicsGuide/
 ***************************************************************************/

#include <immintrin.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>
#include "support.h"

#define ARRAY_SIZE 256*1024*1024
#define COMPARE_VALUE 'z'		
#define ID 999999
#define CONST_CHUNKSIZE 128*1024*1024

static int chunk_start;

typedef struct{
    int id;  //the number of the particular thread
    int length; //the size of the table
    int numt;   //the total number of threads
    char f;
    char t;
    int* count;  //the counter for the particular character
    char* array; //the contents of the table
    int chunk_size;
}arg_struct ;

void initData(char * vector, int size) {
	int i;
	srand (ID);
	for (i = 0; i < size; i++){
		vector[i] = rand()%26+'a';
		//vector[i] = 'a';
	}
}
void verify(int  * counters_CPU,int  * counters_AVX2){
	int j;
	for(j = 0; j < 26; ++j) 
		if (counters_CPU[j]!=counters_AVX2[j]) 
			printf("Does not Verify at value %d\n",j);

}

void printArray(char * charArray, int size) {
	int i;
	for (i = 0; i < size; i++)
		printf("%c", charArray[i]);
	printf("\n");
}

int getChunkStart(void *arguments){
	arg_struct *args=(arg_struct *)arguments;	//we declare our structs arguments   
	chunk_start+=CONST_CHUNKSIZE;	//We move to the next chunk
	return ((chunk_start<args->length)?chunk_start:-1); 
	//if the chunk is bigger than what it needs to be then we return -1 to end it
}

void count_all_chars_serial(char * vector, int  * counters, int size, char from, char to){
	int i;
	for (i=0; i<size; i++){
		if (vector[i]>=from && vector[i]<=to){
			//printf("%c\n",vector[i]);
			counters[vector[i]-'a']++;
		}
	}
}

  
void count_all_chars_AVX2(char * vector, int * counters, int size, char from, char to){
    char j;    
    for(j=from; j<=to; j++){      //a loop starts between from and two
      int k=0;
      int sum=0;
      char c = j;                 //each j is our character
    	__m256i y= _mm256_set1_epi8(c);     //this __m256i variable contains char c  
      while(k<size/32){
        __m256i s= _mm256_set1_epi8(0);   //declaration and initialization of the __m256i sum
        int i;
      	for(i=0;i<128;i++){
    		  __m256i x=  _mm256_loadu_si256((__m256i*)(vector +32*k));    //loads the vector value in x
    		
    		  __m256i z=_mm256_cmpeq_epi8 (x,y);                          //compares x with the character c
          __m256i check = _mm256_and_si256 (z,  _mm256_set1_epi8(1)); //this checks if the value of z is -1
          //if it is, then it returns 1 in check
    		  s = _mm256_add_epi8(s,check);                    //the sum is increased by the value of check, which is 1
    	      k++;
    	  }
        __m128i part1 = _mm256_extracti128_si256 (s,0);   //we split s and the first half is stored in part1
        __m128i part2 = _mm256_extracti128_si256 (s,1);   //the second one is stored in part2
        __m256i p1 = _mm256_cvtepi8_epi16 (part1);        //we extend part1 from 8 bits in 16
        __m256i p2 = _mm256_cvtepi8_epi16 (part2);        //we extend part2 from 8 bits in 16
        __m256i result = _mm256_hadd_epi16 (p1,p2);       //we horizontally add the parts
         result = _mm256_hadd_epi16 (result,result);      //it needs to be horizontally added 4 times
         result = _mm256_hadd_epi16 (result,result);      //because it adds pairs of 4 each time, which means
         result = _mm256_hadd_epi16 (result,result);      //for two parts of 16 bits each it needs to hadd 4 times
         
        sum += (int)(((char*)&result)[0]+((char*)&result)[16]);     //the sum is increased by the addition of 
        //the two parts in result
      } 
       counters[j-'a'] += sum; 
    }    
         
}



pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER; //mutex created
     
void *any_Char_Range_Count_pthreads_worker(void *arguments){
	int i,j,k;
    unsigned int  * counter = (unsigned int *) malloc(sizeof(unsigned int) * 26);
    for(k=0;k<26;k++) counter[k]=0;
  int sum=0;
	arg_struct * args = (arg_struct *) arguments;  //we assign the contents of arguments in the struct
	for(i=0;i<args->length;i++){
  //from the start until the end of the array 
		if(((args->array[i])>=(args->f))&&((args->array[i])<=(args->t))){
			//if the character is between from and to
            counter[(args->array[i])-'a']++;  //the counter is increased by one
		}        
	}     
    
  pthread_mutex_lock(&mutex);  //lock the mutex
  for(j=0;j<26;j++){
     (args->count[j]) += counter[j]; //we indicate that the args->count in the struct is equal to the sum
  } 
  pthread_mutex_unlock(&mutex);  //unlock the mutex
  pthread_exit((void*)arguments); //we exit the function
}

void any_Char_Range_Count_pthreads(char * vector, int * counters, int size, char from, char to, int numThreads){ 
  int i,j;
  char k;
  pthread_t threads[numThreads];   //we create a table with threads using the number of threads
  pthread_attr_t attr;             //we create th attributes
  int rc;
  int length_per_thread=size/numThreads;  //calculation of the size of each table of the thread
  void *status;
  arg_struct *args=(arg_struct *)malloc(sizeof(arg_struct) * numThreads);
  
  //initialization of the struct
    for(i=0;i<numThreads;i++){
    //we assign each attribute of each struct with the respected value 
        args[i].id = i;
        args[i].length = length_per_thread;
        args[i].array = vector+i*length_per_thread;
        args[i].f=from;
        args[i].t=to;
        args[i].count=counters;
        args[i].numt=numThreads;
        args[i].chunk_size=0;
    }
    for(i=0;i<numThreads;i++){ 
    //we call pthread_create for each struct and the threads are created
        rc=pthread_create(&threads[i],&attr,any_Char_Range_Count_pthreads_worker,(void *)&args[i]);
        if (rc) {
  			   printf("ERROR; return code from pthread_create()is %d\n", rc);
  			   exit(-1);
  		  }
    }
    for(i=0; i<numThreads; i++) {
    //we wait for each thread to finish
      rc = pthread_join(threads[i], &status);
      if (rc) {
  			printf("ERROR; return code from pthread_join()is %d\n", rc);
  			exit(-1);
  		}
    }  
    /* Free attribute and wait for the other threads */
    pthread_attr_destroy(&attr);
}

               

void *any_Char_Range_Count_pthreads_AVX2_worker(void *arguments){	
    int i=0,s,j,y;
	arg_struct *args=(arg_struct *)arguments;   	
	startTime(args->id);			//start thread lifecycle
	for(i=args->f; i<=args->t; i++){  	//for each character
		__m256i ch = _mm256_set1_epi8(i); //this __m256i variable contains char i (32 times) 
		s = 0; 
		j = (args->chunk_size*args->id);	
		while(j < (args->chunk_size*(args->id+1))){ 
			__m256i sum = _mm256_set1_epi8(0); //declaration and initialization of the __m256i sum
			__m256i m1 = _mm256_loadu_si256((__m256i*)(args->array + j)); //loads the char value in m1
			__m256i m2 = _mm256_cmpeq_epi8 (m1, ch); 		 //compares m1 with the character ch
      			__m256i flag = _mm256_and_si256 (m2, _mm256_set1_epi8(1)); //this checks if the value of m2 is -1
			sum = _mm256_add_epi8(sum,flag); 			//the sum is increased by the value of flag
			__m128i part0 = _mm256_extracti128_si256 (sum,0);	//we split s and the first half is stored in part0
            		__m128i part1 = _mm256_extracti128_si256 (sum,1); //and the second half is stored in part1
            		__m256i p0 = _mm256_cvtepi8_epi16 (part0);     //we extend part0 from 8 bits in 16
            		__m256i p1 = _mm256_cvtepi8_epi16 (part1);      //we extend part1 from 8 bits in 16
            		__m256i result = _mm256_hadd_epi16 (p0,p1);       //we horizontally add the parts
            	  	result = _mm256_hadd_epi16 (result,result);//it needs to be horizontally added 4 times
            		result = _mm256_hadd_epi16 (result,result);//because it adds pairs of 4 each time, which means
            		result = _mm256_hadd_epi16 (result,result); //for two parts of 16 bits each it needs to hadd 4 times
     			s += (int)(((char*)&result)[0] +((char*)&result)[16]); //the sum is increased by the addition of 
        //the two parts in result
     			j=j+32;
		}
		args->count[i-'a'] += s;
	}
	stopTime(args->id);	//stop thread lifecycle
	pthread_exit(NULL);	//terminate thread
}
void any_Char_Range_Count_pthreads_AVX2( char * theArray, int *counters, int arraySize, char from, char to, int numThreads){
    int i,j,rc;
	pthread_t threads[numThreads];                      //an array of our threads
    pthread_attr_t attr; 
	arg_struct *args[numThreads];                           //an array of our structs
	int chunk_size = arraySize/numThreads;              //the size of each thread so that it will be equally divided between them
     //initialization of the struct
    for (i=0; i < numThreads; i++){                     
     //we assign each attribute of each struct with the respected value
			args[i] = (arg_struct *)malloc(sizeof(arg_struct));  
			args[i]->id = i;                                                                
			args[i]->length=arraySize;                           
			args[i]->array = theArray;                       
			args[i]->f = from;                               
			args[i]->t = to;                                   
            args[i]->numt=numThreads;
			args[i]->count=counters;              
			args[i]->chunk_size = (i==numThreads-1)?arraySize-i*chunk_size:chunk_size; 
            //thread is created  
            pthread_create(&threads[i], &attr, any_Char_Range_Count_pthreads_AVX2_worker, (void*)args[i]);           
    	}
        
    for (i=0; i < numThreads; i++){            //we wait for all the threads to finish 
        pthread_join(threads[i], NULL);             
		          
	}
    return;
}

void *any_Char_Range_Count_pthreads_dynamic_worker(void *arguments)
{	int i,start;
	arg_struct *args=(arg_struct *)arguments;       //we create the struct            
	startTime(args->id);                    //start thread lifecycle
	while(1){
		pthread_mutex_lock(&mutex);               //locks the mutex
		start=getChunkStart((void*)args);	  //we calculate the next start
		pthread_mutex_unlock(&mutex);             //unlocks the mutex
		int until = start+CONST_CHUNKSIZE;  //we search from start until the next chunk_size
		if((start+CONST_CHUNKSIZE)>(args->length)) //the last chunk size could be smaller
			until = args->length;	  //we search from start until the end of the array
		if (start<0) {		//if start is equal to -1 we stop the thread lifecycle
		stopTime(args->id);	
		pthread_exit(NULL);     //terminate thread
		}
		for(i = start; i < until; i++){
			if (args->array[i]>=args->f && args->array[i]<=args->t)  //if character is between from-to
				args->count[args->array[i]-'a']++;               //the counter is added by 1
		}
	}
	stopTime(args->id);               //stop thread lifecycle
	pthread_exit(NULL);            //terminate thread
}

void any_Char_Range_Count_pthreads_dynamic(char * theArray, int *counters, int arraySize, char from, char to,  int numThreads){
	int i=0, j;
	chunk_start=-CONST_CHUNKSIZE;		//chunk_start is equal to -the chunksize because it will be added in getChunkStart 
	pthread_t threads[numThreads];     //an array of threads
	arg_struct *args[numThreads];          //an array of struct package_t with size equal to the number of threads
    //initialization of the struct
    	for (i=0; i < numThreads; i++){    
         //we assign each attribute of each struct with the respected value              
			args[i] = (arg_struct *)malloc(sizeof(arg_struct));    //we create the struct
			args[i]->id = i;                                                          
			args[i]->length=arraySize;                     
			args[i]->array = theArray;                      
			args[i]->f = from;                               
			args[i]->t = to;                                  
            args[i]->numt=numThreads;
			args[i]->count=counters;
            //thread is created
      pthread_create(&threads[i], NULL, any_Char_Range_Count_pthreads_dynamic_worker, (void*)args[i]);      
    }

     for (i=0; i < numThreads; i++){                  //we wait for all the threads to finish
        pthread_join(threads[i], NULL);                  
		
	}
     
}


void *count_all_chars_pthreads_AVX2_dynamic_worker(void *arguments){	
    int i=0,s,j,y;
	arg_struct *args=(arg_struct *)arguments;
	startTime(args->id);           //start thread lifecycle
	int start=0;
	while(1){
	pthread_mutex_lock(&mutex);        //mutex locked
		start=getChunkStart((void*)args); //find the next chunkstart   
		pthread_mutex_unlock(&mutex);   //mutex unlocked
		 if (start<0) {
			pthread_exit(NULL);
		}
        
		if (start+CONST_CHUNKSIZE > args->length) {
	      for(i=args->f; i<=args->t; i++){     //for each character
      		__m256i ch = _mm256_set1_epi8(i);      //this __m256i variable contains char i (32 times) 		
      		s = 0; 
      		j =start;
      		while(j < args->length){ 
      			    __m256i sum = _mm256_set1_epi8(0);        //declaration and initialization of the __m256i sum
      				__m256i m1 = _mm256_loadu_si256((__m256i*)(args->array + j));  //loads the char value in m1
      				__m256i m2 = _mm256_cmpeq_epi8 (m1, ch);           //compares m1 with the character ch
            	    __m256i flag = _mm256_and_si256 (m2, _mm256_set1_epi8(1));       //this checks if the value of m2 is -1
      			   	sum = _mm256_add_epi8(sum,flag);                        //the sum is increased by the value of flag
      				 	
      			          __m128i part0 = _mm256_extracti128_si256 (sum,0);  //we split s and the first half is stored in part0
                  		__m128i part1 = _mm256_extracti128_si256 (sum,1);   //and the second half is stored in part1
                  		__m256i p0 = _mm256_cvtepi8_epi16 (part0);    //we extend part1 from 8 bits in 16
                  		__m256i p1 = _mm256_cvtepi8_epi16 (part1);   //we extend part1 from 8 bits in 16
                  		__m256i result = _mm256_hadd_epi16 (p0,p1);  //we horizontally add the parts
                  	  result = _mm256_hadd_epi16 (result,result);    //it needs to be horizontally added 4 times
                  		result = _mm256_hadd_epi16 (result,result);  //because it adds pairs of 4 each time, which means
                  		result = _mm256_hadd_epi16 (result,result);  //for two parts of 16 bits each it needs to hadd 4 times
      
           			s += (int)(((char*)&result)[0] +((char*)&result)[16]); //the sum is increased by the addition of 
                     //the two parts in result
           			j=j+32;
      		}
      		args->count[i-'a'] += s;      //the final sum is added
      	}
	}else{
	
      	  for(i=args->f; i<=args->t; i++){
      		__m256i ch = _mm256_set1_epi8(i); 
      		                                           //for now on the same procedure as above is executed
      		s = 0;                                     //the only difference is that the chunk size might be
      		j =start;                                  //smaller because it's the last chunk
      		while(j < start+CONST_CHUNKSIZE){ 
      
      			  __m256i sum = _mm256_set1_epi8(0);   
      				__m256i m1 = _mm256_loadu_si256((__m256i*)(args->array + j)); 
      				__m256i m2 = _mm256_cmpeq_epi8 (m1, ch); 
            	__m256i flag = _mm256_and_si256 (m2, _mm256_set1_epi8(1)); 
      				sum = _mm256_add_epi8(sum,flag); 
      				 	
      			          __m128i part0 = _mm256_extracti128_si256 (sum,0);
                  		__m128i part1 = _mm256_extracti128_si256 (sum,1); 
                  		__m256i p0 = _mm256_cvtepi8_epi16 (part0); 
                  		__m256i p1 = _mm256_cvtepi8_epi16 (part1);
                  		__m256i result = _mm256_hadd_epi16 (p0,p1); 
                  	  result = _mm256_hadd_epi16 (result,result);
                  		result = _mm256_hadd_epi16 (result,result); 
                  		result = _mm256_hadd_epi16 (result,result); 
      
           			s += (int)(((char*)&result)[0] +((char*)&result)[16]); 
           			j=j+32;
      		}
      		args->count[i-'a'] += s;      
      	}
	   }
	}
	stopTime(args->id);    //stop thread lifecycle
	pthread_exit(NULL);       //terminate thread
}


void any_Char_Range_Count_AVX2_pthreads_dynamic(char * theArray , int  * counters,int size,char from, char to, int numThreads) {
	pthread_t threads[numThreads];  //an array of threads        
	arg_struct *args[numThreads];      //an array of our structs, the size of our number of threads
    chunk_start=-CONST_CHUNKSIZE; //chunk_start is equal to -the chunksize because it will be added in getChunkStart 
				
    pthread_attr_t attr; 
	                        
	int i,j;
     //initialization of the struct
	for (i=0; i < numThreads; i++){ 
     //we assign each attribute of each struct with the respected value                               
		args[i] = (arg_struct *)malloc(sizeof(arg_struct));      //we create the struct  
		args[i]->id = i;                                                        
		args[i]->array = theArray;                            
		args[i]->length = size;                               
		args[i]->f = from;                                   
		args[i]->t = to;                                       
		args[i]->count=counters;                   
        args[i]->numt=numThreads;
        //thread is created
		pthread_create(&threads[i], &attr, count_all_chars_pthreads_AVX2_dynamic_worker, (void*)args[i]);      
	}

	 for (i=0; i < numThreads; i++){                  //we wait for the threads to finish 
		pthread_join(threads[i], NULL);                  
		
	}
}  

void any_Char_Range_Count_omp(char * vector, int  * counters, int size, char from, char to, int numThreads){
	int i=0;

   //omp operation using the number of threads
	#pragma omp parallel num_threads(numThreads) private(i)
	{
	//private counter
	unsigned int  * count = (unsigned int *) malloc(sizeof(unsigned int) * 26);
		//an omp for starts
		#pragma omp for schedule(static)
		for (i=0; i<size; i++){
				//checks if character is in the specified range
		      if (vector[i]>=from && vector[i]<=to){
					//private counter increased
			     count[vector[i]-'a']++;
		      }
	    }
		//critical section adds the private counter results to counters
		#pragma omp critical
		for(i=0; i<26;i++){
			counters[i]+=count[i];	
		}
	} 

} 

void any_Char_Range_Count_omp_dynamic(char * vector, int  * counters, int size, char from, char to, int numThreads){
	int i=0;

   //omp operation using the number of threads 
	#pragma omp parallel num_threads(numThreads) private(i)
	{
	//private counter
	unsigned int  * count = (unsigned int *) malloc(sizeof(unsigned int) * 26);
		//an omp for starts dynamically with chunk size equal to size/threads 
		#pragma omp for schedule(dynamic, size/numThreads)
		//checks if character is in the specified range
		for (i=0; i<size; i++){
				//private counter increased
		      if (vector[i]>=from && vector[i]<=to){
					//private counter increased
			     count[vector[i]-'a']++;
		      }
	    }
		//critical section adds the private counter results to counters
		#pragma omp critical
		for(i=0; i<26;i++){
			counters[i]+=count[i];	
		}
	} 

} 


void any_Char_Range_Count_omp_AVX2(char * vector, int * counters, int size, char from, char to, int numThreads){
    char j;  
	 //omp operation using the number of threads  
    #pragma omp parallel num_threads(numThreads) private(j)
    {
	 //private counter
    unsigned int  * count = (unsigned int *) malloc(sizeof(unsigned int) * 26);
	 //an omp for starts
    #pragma omp for schedule(static)
    for(j=from; j<=to; j++){      //a loop starts between from and two
      int k=0;
      int sum=0;
      char c = j;                 //each j is our character
    	__m256i y= _mm256_set1_epi8(c);     //this __m256i variable contains char c  
      while(k<size/32){
        __m256i s= _mm256_set1_epi8(0);   //declaration and initialization of the __m256i sum
        int i;
      	for(i=0;i<128;i++){
    		  __m256i x=  _mm256_loadu_si256((__m256i*)(vector +32*k));    //loads the vector value in x
    		
    		  __m256i z=_mm256_cmpeq_epi8 (x,y);                          //compares x with the character c
          __m256i check = _mm256_and_si256 (z,  _mm256_set1_epi8(1)); //this checks if the value of z is -1
          //if it is, then it returns 1 in check
    		  s = _mm256_add_epi8(s,check);                    //the sum is increased by the value of check, which is 1
    	      k++;
    	  }
        __m128i part1 = _mm256_extracti128_si256 (s,0);   //we split s and the first half is stored in part1
        __m128i part2 = _mm256_extracti128_si256 (s,1);   //the second one is stored in part2
        __m256i p1 = _mm256_cvtepi8_epi16 (part1);        //we extend part1 from 8 bits in 16
        __m256i p2 = _mm256_cvtepi8_epi16 (part2);        //we extend part2 from 8 bits in 16
        __m256i result = _mm256_hadd_epi16 (p1,p2);       //we horizontally add the parts
         result = _mm256_hadd_epi16 (result,result);      //it needs to be horizontally added 4 times
         result = _mm256_hadd_epi16 (result,result);      //because it adds pairs of 4 each time, which means
         result = _mm256_hadd_epi16 (result,result);      //for two parts of 16 bits each it needs to hadd 4 times
         
        sum += (int)(((char*)&result)[0]+((char*)&result)[16]);     //the sum is increased by the addition of 
        //the two parts in result
      } 
		//private counter increased
       count[j-'a'] += sum; 
    }
	//critical section adds the private counter results to counters
    #pragma omp critical
		for(j=0; j<26;j++){
			counters[j]+=count[j];	
		} 
    }   
         
}       
                    
void any_Char_Range_Count_omp_AVX2_dynamic(char * vector, int * counters, int size, char from, char to, int numThreads){
    char j;   
	 //omp operation using the number of threads 
    #pragma omp parallel num_threads(numThreads) private(j)
    {
	 //private counter
    unsigned int  * count = (unsigned int *) malloc(sizeof(unsigned int) * 26);
	 //an omp for starts dynamically with chunk size equal to size/threads
    #pragma omp for schedule(dynamic, size/numThreads)
    for(j=from; j<=to; j++){      //a loop starts between from and two
      int k=0;
      int sum=0;
      char c = j;                 //each j is our character
    	__m256i y= _mm256_set1_epi8(c);     //this __m256i variable contains char c  
      while(k<size/32){
        __m256i s= _mm256_set1_epi8(0);   //declaration and initialization of the __m256i sum
        int i;
      	for(i=0;i<128;i++){
    		  __m256i x=  _mm256_loadu_si256((__m256i*)(vector +32*k));    //loads the vector value in x
    		
    		  __m256i z=_mm256_cmpeq_epi8 (x,y);                          //compares x with the character c
          __m256i check = _mm256_and_si256 (z,  _mm256_set1_epi8(1)); //this checks if the value of z is -1
          //if it is, then it returns 1 in check
    		  s = _mm256_add_epi8(s,check);                    //the sum is increased by the value of check, which is 1
    	      k++;
    	  }
        __m128i part1 = _mm256_extracti128_si256 (s,0);   //we split s and the first half is stored in part1
        __m128i part2 = _mm256_extracti128_si256 (s,1);   //the second one is stored in part2
        __m256i p1 = _mm256_cvtepi8_epi16 (part1);        //we extend part1 from 8 bits in 16
        __m256i p2 = _mm256_cvtepi8_epi16 (part2);        //we extend part2 from 8 bits in 16
        __m256i result = _mm256_hadd_epi16 (p1,p2);       //we horizontally add the parts
         result = _mm256_hadd_epi16 (result,result);      //it needs to be horizontally added 4 times
         result = _mm256_hadd_epi16 (result,result);      //because it adds pairs of 4 each time, which means
         result = _mm256_hadd_epi16 (result,result);      //for two parts of 16 bits each it needs to hadd 4 times
         
        sum += (int)(((char*)&result)[0]+((char*)&result)[16]);     //the sum is increased by the addition of 
        //the two parts in result
      } 
		 //private counter increased
       count[j-'a'] += sum; 
    }
	 //critical section adds the private counter results to counters
    #pragma omp critical
		for(j=0; j<26;j++){
			counters[j]+=count[j];	
		} 
    }   
         
} 

int main(int argc, char **argv) {
	int arraySize, sum, i,method;
    int numThreads=0;
	char from, to;
	if (argc == 6) {
		arraySize = atoi(argv[1]);
		from = argv[2][0];
		to = argv[3][0];
        
		numThreads = atoi(argv[4]);
        method = atoi(argv[5]);
	} else {
		printf("Usage: ./a.out NumberOfCharsInMB charFrom charTo NumberOfThreads Method \n");
		printf("Methods: 0=AVX2, 1=pThreads_serial, 2=pThreads_serial_dynamic, 3=pThreads_AVX2, 4=pThreads_AVX2_dynamic, 5=OMP_serial, 6=OMP_serial_dynamic, 7=OMP_AVX2, 8=OMP_AVX2_dynamic. \n");
        return 0;
	}

	printf("%d %c %c\t",arraySize,from,to);
	 __attribute__ ((aligned (256))) char * charArray = (char *) malloc(sizeof(char) * arraySize);
	unsigned int  * counters_CPU = (unsigned int *) malloc(sizeof(unsigned int) * 26);
	unsigned int  * counters_AVX2 = (unsigned int *) malloc(sizeof(unsigned int) * 26);
    unsigned int  * counters_pthreads_CPU = (unsigned int *) malloc(sizeof(unsigned int) * 26);
    unsigned int  * counters_pthreads_CPU_dyn = (unsigned int *) malloc(sizeof(unsigned int) * 26);
    unsigned int  * counters_pthreads_AVX2 = (unsigned int *) malloc(sizeof(unsigned int) * 26);
    unsigned int  * counters_pthreads_AVX2_dyn = (unsigned int *) malloc(sizeof(unsigned int) * 26);
    unsigned int  * counters_omp_CPU = (unsigned int *) malloc(sizeof(unsigned int) * 26);
    unsigned int  * counters_omp_CPU_dyn = (unsigned int *) malloc(sizeof(unsigned int) * 26);
    unsigned int  * counters_omp_AVX2 = (unsigned int *) malloc(sizeof(unsigned int) * 26);
    unsigned int  * counters_omp_AVX2_dyn = (unsigned int *) malloc(sizeof(unsigned int) * 26);
	initData(charArray, arraySize);
	//printArray(charArray, arraySize);
	for(i=0;i<26;i++) counters_CPU[i]=0;
	for(i=0;i<26;i++) counters_AVX2[i]=0;
    for(i=0;i<26;i++) counters_pthreads_CPU[i]=0;
    for(i=0;i<26;i++) counters_pthreads_CPU_dyn[i]=0;
    for(i=0;i<26;i++) counters_pthreads_AVX2[i]=0;
    for(i=0;i<26;i++) counters_pthreads_AVX2_dyn[i]=0;
    for(i=0;i<26;i++) counters_omp_CPU[i]=0;	
    for(i=0;i<26;i++) counters_omp_CPU_dyn[i]=0;
    for(i=0;i<26;i++) counters_omp_AVX2[i]=0;
    for(i=0;i<26;i++) counters_omp_AVX2_dyn[i]=0;
    
    
  	startTime(0);
  	count_all_chars_serial(charArray, counters_CPU, arraySize, from , to);
  	stopTime(0); 
  	for (i=0; i<26; i++) printf("%c Found %d times.\n", 'a'+i,counters_CPU[i]);
  	printf("Time CPU:");elapsedTime(0);	
    
    if(method==0){
        startTime(0);
    	count_all_chars_AVX2(charArray, counters_AVX2, arraySize, from, to);
    	stopTime(0); 
    	for (i=0; i<26; i++) printf("%c Found %d times.\n", 'a'+i,counters_AVX2[i]);
    	printf("\tTime AVX2:");elapsedTime(0);printf("\n");
        
        verify(counters_CPU,counters_AVX2); 
          
    }else if(method==1){
        startTime(0);
    	any_Char_Range_Count_pthreads(charArray, counters_pthreads_CPU, arraySize, from , to, numThreads);
    	stopTime(0); 
    	for (i=0; i<26; i++) printf("%c Found %d times.\n", 'a'+i,counters_pthreads_CPU[i]);
    	printf("Time pthreads:");elapsedTime(0);printf("\n");
        
        verify(counters_CPU,counters_pthreads_CPU);
     
    }else if(method==2){
        startTime(0);
    	any_Char_Range_Count_pthreads_dynamic(charArray, counters_pthreads_CPU_dyn, arraySize, from , to, numThreads);
    	stopTime(0); 
    	for (i=0; i<26; i++) printf("%c Found %d times.\n", 'a'+i,counters_pthreads_CPU_dyn[i]);
    	printf("Time pthreads dynamic:");elapsedTime(0);printf("\n");
      
        verify(counters_CPU,counters_pthreads_CPU_dyn); 
    }else if(method==3){
    
        startTime(0);
    	any_Char_Range_Count_pthreads_AVX2(charArray, counters_pthreads_AVX2, arraySize, from , to, numThreads);
    	stopTime(0); 
    	for (i=0; i<26; i++) printf("%c Found %d times.\n", 'a'+i,counters_pthreads_AVX2[i]);
    	printf("Time pthreads AVX2:");elapsedTime(0);printf("\n");
        
    	verify(counters_CPU,counters_pthreads_AVX2);   
    } else if(method==4){
        startTime(0);
    	any_Char_Range_Count_AVX2_pthreads_dynamic(charArray, counters_pthreads_AVX2_dyn, arraySize, from , to, numThreads);
    	stopTime(0); 
    	for (i=0; i<26; i++) printf("%c Found %d times.\n", 'a'+i,counters_pthreads_AVX2_dyn[i]);
    	printf("Time pthreads AVX2 dynamic:");elapsedTime(0); printf("\n");
        
    	verify(counters_CPU,counters_pthreads_AVX2_dyn);   
    }else if(method==5){
        startTime(0);
    	any_Char_Range_Count_omp(charArray, counters_omp_CPU, arraySize, from , to, numThreads);
    	stopTime(0); 
    	for (i=0; i<26; i++) printf("%c Found %d times.\n", 'a'+i,counters_omp_CPU[i]);
    	printf("Time omp CPU:");elapsedTime(0); printf("\n");
        
    	verify(counters_CPU,counters_omp_CPU);   
    }else if(method==6){
        startTime(0);
    	any_Char_Range_Count_omp_dynamic(charArray, counters_omp_CPU_dyn, arraySize, from , to, numThreads);
    	stopTime(0); 
    	for (i=0; i<26; i++) printf("%c Found %d times.\n", 'a'+i,counters_omp_CPU_dyn[i]);
    	printf("Time omp CPU dynamic:");elapsedTime(0); printf("\n");
        
    	verify(counters_CPU,counters_omp_CPU_dyn);   
    }else if(method==7){
        startTime(0);
      any_Char_Range_Count_omp_AVX2(charArray, counters_omp_AVX2, arraySize, from , to, numThreads);
    	stopTime(0); 
    	for (i=0; i<26; i++) printf("%c Found %d times.\n", 'a'+i,counters_omp_AVX2[i]);
    	printf("Time omp AVX2:");elapsedTime(0); printf("\n");
        
    	verify(counters_CPU,counters_omp_AVX2);   
    }else if(method==8){
        startTime(0);
    	any_Char_Range_Count_omp_AVX2_dynamic(charArray, counters_omp_AVX2_dyn, arraySize, from , to, numThreads);
    	stopTime(0); 
    	for (i=0; i<26; i++) printf("%c Found %d times.\n", 'a'+i,counters_omp_AVX2_dyn[i]);
    	printf("Time omp AVX2 dynamic:");elapsedTime(0); printf("\n");
        
    	verify(counters_CPU,counters_omp_AVX2_dyn);   
    } else {
		printf("Usage: ./a.out NumberOfCharsInMB charFrom charTo NumberOfThreads Method \n");
		printf("Methods: 0=AVX2, 1=pThreads_serial, 2=pThreads_serial_dynamic, 3=pThreads_AVX2, 4=pThreads_AVX2_dynamic, 5=OMP_serial, 6=OMP_serial_dynamic, 7=OMP_AVX2, 8=OMP_AVX2_dynamic. \n");
        return 0;
	}
    
  return 0;
}