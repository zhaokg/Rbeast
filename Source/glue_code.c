#include "abc_000_warning.h"

#define IMPORT_NUMPY
#include "abc_001_config.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
 
#ifndef ARM64_OS
	#include <immintrin.h> //https://stackoverflow.com/questions/56049110/including-the-correct-intrinsic-header
    #include "abc_math_avx.h"
#endif

#include "abc_datatype.h"
#include "abc_blas_lapack_lib.h"
#include "abc_common.h"
#include "abc_ide_util.h"
#include "abc_pthread.h"
#include "abc_timer.h"
#include "abc_vec.h"
#include "abc_date.h"

#include "beastv2_io.h"

#if  ( !defined(R_RELEASE) && !defined(M_RELEASE) && !defined(P_RELEASE)  ) || defined( PLAY_MODE)
#include "mrbeast_header.h"
#include "mrbeast_io.h" 
#include "sbmfast.h"
#include "sbmfast_io.h"
#endif

#include "globalvars.h"

#if defined(WIN64_OS) 
//extern void DllExport WinMainDemoST(BEAST_OPTIONS_PTR  option);
//extern void DllExport WinMainDemoTrend(BEAST_OPTIONS_PTR  option);
#endif


#include "abc_cpu.h"
#include "abc_timer.h"
#include "abc_rand.h" 

#define __IS_STRING_EQUAL(a, b)  (strcicmp(a, #b) == 0)
 

static void  GetArg_IsQuiteMode(VOIDPTR prhs[], int nrhs) {

 	if (nrhs >= 6 && IsStruct( prhs[5L] ) ) {
 			VOIDPTR tmp;
			GLOBAL_QUIET_MODE = (tmp = GetField123Check(prhs[5L], "quiet", 3)) ? GetScalar(tmp) : 0L;
			return;	 
	} // if (nrhs >= 5)
	GLOBAL_QUIET_MODE = 0;
	return;
}



void * mainFunction(void *prhs[], int nrhs) {

	if (nrhs >= 7) {
		int avxOption = GetScalar(prhs[nrhs - 1]); 
		SetupRoutines_UserChoice(avxOption);
	}	else {
		GetArg_IsQuiteMode(prhs, nrhs);      // Set GLOBAL_QUITE_MODE
		SetupRoutines_AutoByCPU(1L);          // prnitlevel=1: only print hhe libray choice result
	}
	//SetupVectorFunction_Generic();
	//SetupPCG_GENERIC();	 	
	//print_funcs();

	if (nrhs == 0 ) 	{
		r_error("ERROR: Essential input paramaters are missing!\n");
		//r_printf("\033[1m\033[32m" "\033[4m" "%cERROR: Essential input paramaters are missing!\n",149);
		return IDE_NULL;
	}
	if ( !IsChar(prhs[0]) )	{
		r_error("ERROR: The very first parameter must be a string specifying the algorithm name!\n");
		return IDE_NULL;
	}

	//mexCallMATLAB(0, NULL, 2, prhs, "fprintf");

	#define __STRING_LEN__ 20
	char  algorithm[__STRING_LEN__ +1];
	GetCharArray(prhs[0], algorithm, __STRING_LEN__);

	// Now 'algorithm' should be eitehr the default algorith or the value from Opt.

	void * ANS  = NULL;
	int    nptr = 0;
	if      (__IS_STRING_EQUAL(algorithm, beastv4Demo)) 	{

		#ifdef WIN64_OS
			#if P_INTERFACE ==1
					// Covert the second arg into a Numpy Array. the pointer returend
					// is a new ref that MUST be dec-refed at the end
					prhs[1] = CvtToPyArray_NewRef(prhs[1]);
			#endif

			//BEAST2_OPTIONS    option = {0,}; 
			BEAST2_OPTIONS      option = { {{0,},}, }; //Warning from MacOS: suggest braces around initialization of subobject [-Wmissing-braces]
			option.io.out.result	 = NULL;		// result to be allocated in OUput_allocMEM
			GLOBAL_OPTIONS           = (BEAST2_OPTIONS_PTR)&option;
				
			if ( BEAST2_GetArgs(prhs,nrhs,&option) ==0 ) {
				BEAST2_DeallocateTimeSeriesIO(&(option.io));
				return IDE_NULL;
			}
		
	
			option.extra.computeCredible = TRUE;
			ANS = PROTECT(BEAST2_Output_AllocMEM(&option)); nptr++;
			
			void DllExport BEAST2_WinMain(VOID_PTR  option);
			BEAST2_WinMain((BEAST2_OPTIONS_PTR)GLOBAL_OPTIONS);
			BEAST2_DeallocateTimeSeriesIO(&(option.io));

			#if P_INTERFACE ==1
					Py_XDECREF(prhs[1]);
			#endif
		#else
			r_printf("WARNING: The GUI interface is supporte only on the Windows 64 operating system.\n");
		#endif	

	}
	else if (__IS_STRING_EQUAL(algorithm, beastv4))
	{
		#if P_INTERFACE ==1
		// Covert the second arg into a Numpy Array. the pointer returend
		// is a new ref that MUST be dec-refed at the end
			prhs[1] = CvtToPyArray_NewRef(prhs[1]);
		#endif
		/*
		// Initialize mutex and condition variable objects
		pthread_mutex_init(&MUTEX_WRITE, NULL);
		pthread_cond_init(&CONDITION_WRITE, NULL);
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		DATA_AVAILABLE_WRITE = false;
		pthread_create(&THREADID_WRITE, &attr, WRITE, (VOID_PTR )&threadParWrite);
		pthread_attr_destroy(&attr); 	
		*/

		/**********************************/
	    //stackoverflow.com/questions/10828294/c-and-c-partial-initialization-of-automatic-structure
	    //stackoverflow.com/questions/11152160/initializing-a-struct-to-0
		//BEAST2_OPTIONS      option = {0,};		
		BEAST2_OPTIONS      option = { {{0,},}, }; //Warning from MacOS: suggest braces around initialization of subobject [-Wmissing-braces]
	
		if (BEAST2_GetArgs(prhs, nrhs, &option) == 0) {
			BEAST2_DeallocateTimeSeriesIO(&(option.io));
		    #if P_INTERFACE ==1
				Py_XDECREF(prhs[1]);
		    #endif
			return IDE_NULL;
		}		
		
		option.io.out.result		 = NULL;		 	// result to be allocated in OUput_allocMEM

		if (option.io.q == 1) {
			ANS = PROTECT(BEAST2_Output_AllocMEM(&option)); nptr++;	
		} else {			
			//memset(&option.extra, 0, sizeof(option.extra));
			option.extra.computeSeasonAmp = 0;
			option.extra.computeTrendSlope = 0;
			option.extra.tallyIncDecTrendJump= 0;
			option.extra.tallyPosNegTrendJump = 0;
			option.extra.tallyPosNegOutliers = 0;
			option.extra.tallyPosNegSeasonJump = 0;
 
			option.extra.computeTrendChngpt = 1;
			option.extra.computeSeasonChngpt = 1;
			//option.extra.computeOutlierChngpt = 1;
			BEAST2_print_options(&option); 
			ANS = PROTECT(BEAST2_Output_AllocMEM(&option)); nptr++;
		}
		/**********************************/
	
		GLOBAL_OPTIONS = (BEAST2_OPTIONS_PTR)&option;
		if (option.io.numOfPixels ==1) {
			beast2_main_corev4();
			BEAST2_DeallocateTimeSeriesIO(&(option.io));
			r_printf("\n");
		} else {
			/**********************************/
			I32 NUM_THREADS			= option.extra.numParThreads;  // I32 NUM_CORES_TOUSE= option.extra.numCPUCoresToUse;;			
			I32 NUM_THREADS_PER_CPU = option.extra.numThreadsPerCPU;
			I32 NUM_CORES           = GetNumCores();

			NUM_CORES	        = max(NUM_CORES, 1L);			
			NUM_THREADS_PER_CPU = max(NUM_THREADS_PER_CPU, 1L);
			
			NUM_THREADS = (NUM_THREADS <= 0) ? NUM_CORES * NUM_THREADS_PER_CPU : NUM_THREADS;
			NUM_THREADS = min(NUM_THREADS, option.io.numOfPixels);

			/**********************************/
			NUM_OF_PROCESSED_PIXELS		 = 0;
			NUM_OF_PROCESSED_GOOD_PIXELS = 0;
			NEXT_PIXEL_INDEX			 = 1;
						
			//Initialize mutex and condition variable objects
			// threadID,  mutex, and condvaar are  global variables.
			pthread_mutex_init(&mutex, NULL);
			pthread_cond_init(&condVar, NULL);
						
			// For portability, explicitly create threads in a joinable state 
			pthread_attr_t attr;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

			thread_id = malloc(sizeof(pthread_t) * NUM_THREADS);

  		   int8_t *thread_stat = malloc(sizeof(int8_t) * NUM_THREADS);
			for (I32 i = 0; i < NUM_THREADS; i++) {
             #if defined(LINUX_OS) || defined (WIN32_OS) || defined (WIN64_OS) 
				cpu_set_t cpuset;
				CPU_ZERO(&cpuset);
				CPU_SET( i%NUM_CORES, &cpuset);
				//https://stackoverflow.com/questions/25472441/pthread-affinity-before-create-threads
				pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset);
				thread_stat[i]=pthread_create( &thread_id[i], &attr, beast2_main_corev4_mthrd, (void*)NULL);
				//https://stackoverflow.com/questions/1407786/how-to-set-cpu-affinity-of-a-particular-pthread
				//sched_setaffinity(thread_id[i], sizeof(cpuset), &cpuset)
			 #elif defined(MAC_OS)
				cpu_set_t cpuset;
				CPU_ZERO(&cpuset);
				CPU_SET(i % NUM_CORES, &cpuset);
				//https://stackoverflow.com/questions/25472441/pthread-affinity-before-create-threads			
				thread_stat[i] = pthread_create(&thread_id[i], &attr, beast2_main_corev4_mthrd, (void*)NULL);
				pthread_setaffinity_np(thread_id[i], sizeof(cpu_set_t), &cpuset);
		     #elif defined (SOLARIS_OS)
				cpu_set_t cpuset;
				CPU_ZERO(&cpuset);
				CPU_SET(i % NUM_CORES, &cpuset);
				//https://stackoverflow.com/questions/25472441/pthread-affinity-before-create-threads			
				thread_stat[i] = pthread_create(&thread_id[i], &attr, beast2_main_corev4_mthrd, (void*)NULL);
				sched_getaffinity(thread_id[i], sizeof(cpu_set_t), &cpuset);
			 #else
				thread_stat[i] = pthread_create(&thread_id[i], &attr, beast2_main_corev4_mthrd, (void*)NULL);
			 #endif
				if ( 0 == thread_stat[i]) {
					r_printf("Parallel computing: thread#%-02d generated ... \n", i + 1);
			    } else {
					r_printf("Parallel computing: thread#%-02d failed to generate ... \n", i + 1);
				}
				
			}
			r_printf("Rbeast: Waiting on %d threads...\n", NUM_THREADS);
			pthread_attr_destroy(&attr);


			IDE_USER_INTERRUPT = 0;
			#if R_INTERACE==1
				r_printf("Press and hold the ESCAPE key or the STOP button to interrupt and quit while running.\n" );
			#elif M_INTERFACE==1
				r_printf("Press and hold CTR+C to interrupt and quit while running.\n");
			#endif

			if (option.extra.printProgressBar) {
				// Print a blank line to be backspaced by the follow
				void* StrBuf      = malloc(option.extra.consoleWidth * 3);
				PERCENT_COMPLETED = 0;
				REMAINING_TIME    = 10000;
				printProgress2(PERCENT_COMPLETED, REMAINING_TIME, option.extra.consoleWidth, StrBuf, 1);

				// https://www.mathworks.com/matlabcentral/answers/101658-is-it-possible-to-start-new-threads-from-a-c-mex-file
				// https://stackoverflow.com/questions/28227313/multithreaded-pthreads-matlab-mex-function-causes-matlab-to-crash-after-exitin
				// https://stackoverflow.com/questions/54010898/how-do-you-print-to-console-in-a-multi-threaded-mex-function
				while (PERCENT_COMPLETED < 1.f && NEXT_PIXEL_INDEX < option.io.numOfPixels && IDE_USER_INTERRUPT==0) {
										
					printProgress2(PERCENT_COMPLETED, REMAINING_TIME, option.extra.consoleWidth, StrBuf, 0);
					if (CheckInterrupt()) {
						ConsumeInterruptSignal();// only needed for Matlab
						IDE_USER_INTERRUPT = 1;
						r_printf("Quitting due to unexpected user interruption...\n");
					}					
					Sleep_ms(2 * 1000);
				}

				if (IDE_USER_INTERRUPT == 0) {
					printProgress2(1.0, 0, option.extra.consoleWidth, StrBuf, 0);
				}

				free(StrBuf);
			}  

			// Wait for all threads to complete
			r_printf("\nFinalizing ... \n");
			for (I32 i = 0; i < NUM_THREADS; i++) {
				//pthread_join(thread_id[i], NULL);
				if (thread_stat[0] == 0) {
				// threads that were gnerated successfully
					I64 ret = 0;
					pthread_join(thread_id[i], &ret);
					//r_printf("\nstack size %d.\n", ret/1024/1024);
					r_printf("Parallel computing : Thread # % -02d finished ... \n", i);
				}				
			}
			if (IDE_USER_INTERRUPT==0)
				r_printf("\nRbeast: Waited on %d threads. Done.\n", NUM_THREADS);
			else
				r_printf("\nQuitted unexpectedly upon the user's interruption.\n");

			// Clean up and exit		
			pthread_mutex_destroy(&mutex);
			pthread_cond_destroy(&condVar);
			// pthread_exit(NULL);	

			free(thread_id);
			free(thread_stat);

			BEAST2_DeallocateTimeSeriesIO(&(option.io));
		}
		
	     #if P_INTERFACE ==1
				Py_XDECREF(prhs[1]);
		 #endif
 
	}
	#if !defined(R_RELEASE) &&  !defined(M_RELEASE) &&  !defined(P_RELEASE)
	else if  (__IS_STRING_EQUAL(algorithm, "mrbeast"))
	{
		MV_OPTIONS         option;
		MV_IO              io;
		MV_RESULT          result;
		memset(&io,     0, sizeof(MV_IO));
		memset(&result, 0, sizeof(MV_RESULT));
	 
		option.io               = &io;
		option.io->out.result   = &result;
		option.io->isRegularOrdered  = 1;
		
		if (!MV_Get1stArg_Data(prhs, nrhs, &option)     ||
			!MV_Get2ndArg_MetaData(prhs, nrhs, &option) ||
			!MV_Get3rdArg_Prior(prhs, nrhs, &option) ||
			!MV_Get4thArg_MCMC(prhs, nrhs, &option) ||
			!MV_Get5thArg_FLAGS(prhs, nrhs, &option)) {
			return IDE_NULL;
		}

		MV_print_options(&option);
		ANS = PROTECT(MV_AllocateOutput(&option)); nptr++;

		GLOBAL_OPTIONS = (MV_OPTIONS_PTR)&option;
		mrbeast_main_core();
		MV_DeallocateTimeSeriesIO(option.io);
	} 
	#endif
	else if (__IS_STRING_EQUAL(algorithm, tsextract)) {
		extern void* BEAST2_TsExtract(void* o, void* pindex);

		//prhs[o] is "tsextract"
		ANS = PROTECT(BEAST2_TsExtract(prhs[1], prhs[2]));
		nptr++;
	}
	else if (__IS_STRING_EQUAL(algorithm, svd)) {
		typedef struct SVDBasisMEM {
			int N, P, Ncycle;
			I32PTR nPtsPerTime; //P
			F32PTR mean;       //P for both XsumPerTime and Mean
			I32PTR goodTimeIndices; //P
			F32PTR Ytrue;// Ncycley*P;
			F32PTR Ycur;//Ncyle(*P
			I32PTR Igood;// Ncycley*P;
			I32PTR NgoodPerPeriod;// Ncycley;
			F32PTR A, B; //P


			F32PTR M; // PXP
			F32PTR Mcopy; // PXP
			F32PTR XtX; // P*P;
			F32PTR Bcoeff; // P;

		} SVDBasisMEM;

	 
		void compute_seasonal_svdbasis(F32PTR y, F32PTR Yout, int  Kmax, SVDBasisMEM* mem);
		float Y[] = {2.2859,2.0691,1.8224,1.5187,1.5328,2.0847,2.8942,3.2902,3.728,5.8975,6.5847,10.2307,10.9694,17.5811,14.1326,11.5202,10.2331,6.6177,5.0504,4.2145,1.8737,3.197,2.8468,2.5349,1.6543,1.3169,1.8825,2.2085,1.8276,2.5948,3.6769,4.822,4.8609,6.4938,6.9191,10.2873,11.2253,13.0343,14.1471,15.8413,11.9918,7.2993,5.5117,6.621,5.5967,2.8521,3.9991,3.2917,4.342,1.6677,1.9921,2.329,1.494,1.697,2.6342,2.5178,2.8399,4.4247,6.831,10.3966,13.3897,14.5965,14.5742,12.6073,11.3888,10.045,9.0508,5.148,4.2168,2.952,2.7271,2.3879,2.0889,1.5052,2.004,1.7263,1.4218,1.8645,3.7781,3.2067,2.5755,3.0108,5.7034,8.6385,9.9321,11.321,13.1164,13.2138,10.3973,8.117,6.3119,5.269,6.7344,3.2315,2.6521,2.6256,2.0983,1.844,1.9795,1.4476,1.8262,1.7339,2.4015,2.1856,3.1992,3.5943,5.1744,7.0701,10.6768,12.4415,15.6008,14.9724,10.063,8.6803,9.1394,6.0462,4.1622,2.6468,3.2652,2.0314,2.2326,1.8504,3.217,1.9163,1.1103,1.8004,2.5058,2.0672,3.3731,3.3932,6.4282,6.4809,10.5954,12.5961,15.711,14.2856,12.7219,12.4061,7.4394,9.5144,5.696,3.1871,3.2476,3.844,2.523,2.0841,2.1426,1.7278,2.0496,1.7432,1.8544,2.4466,3.3078,2.5855,4.2597,6.4083,7.4898,11.4539,12.9476,13.7927,15.2775,13.0988,8.255,7.3098,6.1615,3.3724,1.9019,2.6075,4.4399,2.2925,2.268,2.0642,1.9616,2.0934,2.0652,1.974,2.6556,3.2962,3.8939,4.9837,7.8963,11.7756,11.005,13.3557,11.9917,11.4185,9.9252,8.1113,5.6223,4.4981,3.5695,2.3075,2.5229,2.3264,1.9057,1.4906,1.5844,2.3672,2.0908,2.3647,3.0101,2.4715,4.5199,6.8071,7.7233,9.065,12.7859,14.8723,14.2513,16.931,12.1172,8.4426,6.9039,6.276,5.3299,3.2293};
		
		
		for (int i = 0; i < 9*24; i++) {
			Y[i] = getNaN();
		}
		Y[10] = 1;
		Y[20] = 6;
		//Y[111] = getNaN();
		MemPointers MEM = { .init = mem_init };
		MEM.init(&MEM);
		SVDBasisMEM s = {.N=24*9, .P=24};

		I64 Get_Alloc_SVDBasisMEM(int N, int P, SVDBasisMEM* s, VOID_PTR bufBase); // from beastv2_svdbasis

		I32 totalSize=Get_Alloc_SVDBasisMEM(24*9, 24, &s, NULL );
		VOIDPTR buf = malloc(totalSize);
		Get_Alloc_SVDBasisMEM(24 * 9, 24, &s, buf);

		F32PTR Yout;
		ANS = CreateF64NumMatrix(24*9, 24, &Yout);

		compute_seasonal_svdbasis(Y,  Yout, 9, &s);
		f32_to_f64_inplace(Yout, 24 *9*24);
		free(buf);
		MEM.free_all(&MEM);
	}
	else if (__IS_STRING_EQUAL(algorithm, print)) {
		extern void* BEAST2_PrintResult(void* o, void* pindex);
		//prhs[o] is "tsextract"
		BEAST2_PrintResult(prhs[1], prhs[2]); 
	}
	else if (__IS_STRING_EQUAL(algorithm, disp)) { 
	//prhs[o] is "disp"
		IDEPrintObject(prhs[1]); 
	}
	else if (__IS_STRING_EQUAL(algorithm, datenum)) {
		int y = GetScalar(prhs[1]);
		int m = GetScalar(prhs[2]);
		int d = GetScalar(prhs[3]);
		int yorg = GetScalar(prhs[4]);
		int morg = GetScalar(prhs[5]);
		int dorg = GetScalar(prhs[6]);

		int org = JulianDayNum_from_civil_ag1(yorg, morg, dorg);
		int JDN = JulianDayNum_from_civil_ag1(y, m, d);
		int dnum = JDN - org;
		r_printf("datenum: %d", dnum);		
	}
	else if (__IS_STRING_EQUAL(algorithm, cpu)) {
       #if defined(WIN64_OS) || defined(WIN32_OS)
		  int GetCPUInfo();
		  GetCPUInfo();
        #endif
	}
	/*
		// Initialize mutex and condition variable objects
		pthread_mutex_init(&MUTEX_WRITE, NULL);
		pthread_cond_init(&CONDITION_WRITE, NULL);
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		DATA_AVAILABLE_WRITE = false;
		pthread_create(&THREADID_WRITE, &attr, WRITE, (VOID_PTR )&threadParWrite);
		pthread_attr_destroy(&attr); 	
	*/
	/*
	else if (__IS_STRING_EQUAL(algorithm, beastMT))
	{
	

		BEAST_OPTIONS         option;
		BEAST_IO              io;
		BEAST_RESULT          result;
		memset(&io,     0, sizeof(BEAST_IO));
		memset(&result, 0, sizeof(BEAST_RESULT));
		option.io = &io;
		option.io->out.result  = &result;
		option.io->isRegularOrdered = 1;


		if (!BEAST_Get1stArg_Data(prhs, nrhs, &option) ||
			!BEAST_Get2ndArg_MetaData(prhs, nrhs, &option) ||
			!BEAST_Get3rdArg_Prior(prhs, nrhs, &option) ||
			!BEAST_Get4thArg_MCMC(prhs, nrhs, &option) ||
			!BEAST_Get5thArg_FLAGS(prhs, nrhs, &option)) {
			return NULL;
		}

		BEAST_print_options(&option);
		ANS = PROTECT(BEAST_AllocateOutput(&option)); nptr++;
		GLOBAL_OPTIONS = (BEAST_OPTIONS_PTR)&option;
	 
		
 
		pthread_mutex_init(&mutex, NULL); //Initialize mutex and condition variable objects
		pthread_cond_init(&condVar, NULL);

		// For portability, explicitly create threads in a joinable state 
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
 
		cpu_set_t cpus;

		
		I32 NUM_THREADS;

		I32 nCores   = GetNumCores();
		nCores       = max(nCores, 1);

		NUM_THREADS  = option.extra.numCPUCoresToUse;
		if (NUM_THREADS == 0) {
			NUM_THREADS = nCores - 1;
		} else if (NUM_THREADS < 0)	{
			NUM_THREADS = nCores + NUM_THREADS;
		}
		NUM_THREADS = min(nCores- 1, NUM_THREADS);
		NUM_THREADS = max(NUM_THREADS, 1);
		NUM_THREADS = min(NUM_THREADS, option.io->numOfPixels);	

		thread_id    = malloc(sizeof(pthread_t) * NUM_THREADS);

		NEXT_PIXEL_INDEX = 1;
		for (I32 i = 0; i < NUM_THREADS; i++) 	{	
			 CPU_ZERO(&cpus);
			 CPU_SET(i, &cpus);
			 pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);

			 extern 	int beast_main_mthrd(void *);
			 pthread_create(&thread_id[i], &attr, beast_main_mthrd, (void *)NULL);
		}
		pthread_attr_destroy(&attr);

		if (option.extra.printProgressBar) {
			// Print a blank line to be backspaced by the follow
			void* BUF;
			BUF = malloc(option.extra.consoleWidth * 3);
			printProgress2(0, 0, option.extra.consoleWidth, BUF, 1);

			PERCENT_COMPLETED = 0;
			REMAINING_TIME    = 10000;
			// https://www.mathworks.com/matlabcentral/answers/101658-is-it-possible-to-start-new-threads-from-a-c-mex-file
			// https://stackoverflow.com/questions/28227313/multithreaded-pthreads-matlab-mex-function-causes-matlab-to-crash-after-exitin
			// https://stackoverflow.com/questions/54010898/how-do-you-print-to-console-in-a-multi-threaded-mex-function
			while (PERCENT_COMPLETED < 1.f && NEXT_PIXEL_INDEX < option.io->numOfPixels) {			
				printProgress2(PERCENT_COMPLETED, REMAINING_TIME, option.extra.consoleWidth, BUF, 0);
				Sleep_ms(2 * 1000);
			}
			printProgress2(1.0, 0, option.extra.consoleWidth, BUF, 0);
			free(BUF);
		} else	{
			r_printf("\nRbeast: Waiting on %d threads...\n", NUM_THREADS);
		}


		// Wait for all threads to complete
		for (I32 i = 0; i<NUM_THREADS; i++) {
			pthread_join(thread_id[i], NULL);
		}
		r_printf("\nRbeast: Waited on %d threads. Done.\n", NUM_THREADS);

		//Clean up and exit		
		pthread_mutex_destroy(&mutex);
		pthread_cond_destroy(&condVar);
		//pthread_exit(NULL);	
		free(thread_id);
	}
	 */

 
	//else if (IS_STRING_EQUAL(algorithm, "SBM_ST"))			//SBM_ST(nlhs, plhs, nrhs, prhs);	
 	//else if (IS_STRING_EQUAL(algorithm, "SBM_ST_BIC"))  	//SBM_ST_BIC(nlhs, plhs, nrhs, prhs);
	
#ifdef PLAY_MODE   
	else if (__IS_STRING_EQUAL(algorithm, SLIDE))
	{
		SOptions beastOption = { 0, };
		GLOBAL_OPTIONS = &beastOption;
		
		MemPointers MEM = (MemPointers){ .init = mem_init };
		MEM.init(&MEM);
		if (!GetArgs(prhs,nrhs,& beastOption, &MEM)) {
			MEM.free_all(&MEM);
			return IDE_NULL;
		}

		void SBM_print_options(SOptions * opt);
		SBM_print_options(&beastOption);

		DATAL0 data0;   GLOBAL_DATA0 = &data0;
		AllocInitDataL0(&data0, &beastOption, &MEM);


		//print_options(&beastOption);
		////plhs[0] = SF_allocate_output(&result, &beastOption);
		ANS = SBM_Output_AllocMEM(&beastOption);
		sbmfast_slide();

		free(beastOption.io);
		MEM.free_all(&MEM);
		
	} 
 
	else if (__IS_STRING_EQUAL(algorithm, QR))
	{
		SOptions beastOption = { 0, };
		GLOBAL_OPTIONS = &beastOption;
 

		MemPointers MEM = (MemPointers){ .init = mem_init };
		MEM.init(&MEM);
		if (!GetArgs(prhs, nrhs, &beastOption, &MEM)) {
			MEM.free_all(&MEM);
			return IDE_NULL;
		}

		void SBM_print_options(SOptions * opt);
		SBM_print_options(&beastOption);

		DATAL0 data0;   GLOBAL_DATA0 = &data0;
		AllocInitDataL0(&data0, &beastOption, &MEM);

		//print_options(&beastOption);
		////plhs[0] = SF_allocate_output(&result, &beastOption);
 
		ANS = SBM_Output_AllocMEM(&beastOption);
		sbmfast_slide_robust();

		free(beastOption.io);
		MEM.free_all(&MEM);
	}
	else if (__IS_STRING_EQUAL(algorithm, SWEEP))
	{
		SOptions beastOption = {0,};  GLOBAL_OPTIONS = &beastOption;
		
		MemPointers MEM = (MemPointers){ .init = mem_init };
		MEM.init(&MEM);
		if (!GetArgs(prhs,nrhs,& beastOption, &MEM)) {
			MEM.free_all(&MEM);
			return IDE_NULL;
		}

		void SBM_print_options(SOptions * opt);
		SBM_print_options(&beastOption);

		DATAL0 data0; 
		GLOBAL_DATA0 = &data0;
		AllocInitDataL0(&data0, &beastOption, &MEM);
		

		//print_options(&beastOption);
		////plhs[0] = SF_allocate_output(&result, &beastOption);
		ANS = SBM_Output_AllocMEM(&beastOption);
		sbmfast_sweep();

		free(beastOption.io);
		MEM.free_all(&MEM);
	}
#endif
	/*
	else if (strcicmp( algName, "sbm") == 0)
	{
	OPTIONS beastOption;
	RESULT  result;
	GLOBAL_OPTIONS = (OPTIONS_PTR)&beastOption;
	GLOBAL_RESULT = (RESULT_PTR)&result;

	if (!M_sbm_check_input(nrhs, prhs) ||
	!M_sbm_read_input(&beastOption, nrhs, prhs,  missing) ||
	!check_options(&beastOption,  missing)
	)
	{
	return;
	}

	print_options(&beastOption);
	plhs[0] = M_sbm_allocate_output(&result, &beastOption);
	sbm();
	}
	*/
	
 
	UNPROTECT(nptr);

	return ANS == NULL ? IDE_NULL : ANS;	
}
 

#if R_INTERFACE==1
#include <R_ext/libextern.h>
#include "Rembedded.h"

//	R_FlushConsole(): a R functin to flush the print
#if defined(MSVC_COMPILER)
SEXP DllExport rexFunction1(SEXP rList, SEXP dummy) {
#else
SEXP DllExport rexFunction(SEXP rList, SEXP dummy) {
#endif
	if (!isNewList(rList)) 	return R_NilValue; 
	// stat.ethz.ch/pipermail/r-help/2001-September/015081.html
	// SEXP   prhs = VECTOR_DATA(inList);
	SEXP   prhs[10];
	int    nrhs = length(rList);
	nrhs = min(10L, nrhs);
	for (int i = 0; i < nrhs; i++) 	
		prhs[i] = VECTOR_ELT(rList, i);

	SEXP ans;
	PROTECT(ans=mainFunction(prhs, nrhs));
	UNPROTECT(1);
	return ans == NULL ?  R_NilValue: ans;
}

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}
#if (defined(WIN64_OS) || defined(WIN32_OS)) 
	SEXP TetrisSetTimer(SEXP action, SEXP seconds, SEXP envior);
	static const R_CallMethodDef CallEntries[] = {
		#if defined(MSVC_COMPILER)
			CALLDEF(rexFunction1,    2),
		#else
			CALLDEF(rexFunction,    2),
		#endif
		CALLDEF(TetrisSetTimer, 3),
		{ NULL, NULL, 0 }
	};
#else
	static const R_CallMethodDef CallEntries[] = {
					CALLDEF(rexFunction, 2),
					{ NULL, NULL, 0 }
			};
#endif


void  R_init_Rbeast(DllInfo *dll) {
	  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
	  R_useDynamicSymbols(dll, FALSE);
	  //R_forceSymbols(dll, TRUE);
}

#elif M_INTERFACE == 1

/*******************************************************************************************************
* function: mexFunction() - Matlab interface function.         
* INPUTS:                                                      
*   nlhs - The number of output variables to be assigned by  the mex function.       
*   plhs[] - Empty array of MATLAB matricies, of size nlhs. This is to be filled in by this application.
*   nrhs - Number of input arguments to this mex funciton. 
*   prhs[] - Array of MATLAB matricies from which input data  is taken.
*******************************************************************************************************/

#include "abc_date.h"
// void DllExport mexFunction(int nlhs, mxArray * _restrict plhs[],   int nrhs, const mxArray * _restrict prhs[]) {
/* 
   Restirct is really pf no use here. More important, the decleration of mexFunction in mex.h has no restrict keyword;
   MVSC is OK with the difference; gcc, however, failed with a complaint of conflicting types for ‘mexFunction’. 
   What flag can circuvume this? The best thing I could find is
   https://stackoverflow.com/questions/66951324/conflicting-types-compiling-a-ld-preload-wrapper
*/
void DllExport mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

	
	//https://stackoverflow.com/questions/19813718/mex-files-how-to-return-an-already-allocated-matlab-array
	//https://www.mathworks.com/matlabcentral/answers/77048-return-large-unchange-mxarray-from-mex
	//https://stackoverflow.com/questions/18847833/is-it-possible-return-cell-array-that-contains-one-instance-in-several-cells/18849127#18849127
	//https://undocumentedmatlab.com/articles/matlabs-internal-memory-representation/
	// NOT ALWAYS WORK
	//plhs[0] = prhs[0];
	//mxCreateSharedDataCopy
	
	mxArray * ans = mainFunction(prhs, nrhs);
	plhs[0] = ans;
	return;


	//F32PTR T = mxGetData(prhs[0]);
	//F32PTR Y = mxGetData(prhs[1]);

	//I32  N        = max(mxGetM(prhs[0]), mxGetN(prhs[0]));

	//chol_full(A, U, N, N);
	//pack_chol_update(x, U1, N);
	//inplace_chol_addCol(A,  N, 1, N);
	//pack_chol(B, U1, N);	
	//pack_chol(B, U, N);
	//pack_solve_L(U, x, N);
	//pack_solve_U(U, x, N);

	//memcpy(data.L2.Y, data.Y, sizeof(F32)*data.L1.n);
	//memcpy(data.Y0.Y, data.Y, sizeof(F32)*data.L1.n);

	//plhs[0] = mxCreateNumericMatrix(data.S.KMAX, N, mxSINGLE_CLASS, mxREAL);
	//F32PTR O1 = mxGetData(plhs[0]);

	//U: col (1, N-w) row(1+w, N)

	/*
	data.L4.sKnot[0] = 101, data.L4.sKnot[1] = 201, data.L4.sKnot[2] = 301;
	data.L4.tKnot[0] = 101, data.L4.tKnot[1] = 201, data.L4.tKnot[2] = 301;

	data.L4.sOrder[0] = 2, data.L4.sOrder[1] = 1, data.L4.sOrder[2] = 3;
	data.L4.tOrder[0] = 2, data.L4.tOrder[1] = 2, data.L4.tOrder[2] = 2;
	data.L4.sNum = 2;
	data.L4.tNum = 2;
	*/

	//get_U_s_e(&data, 1, data.n, 4, 1);
	//memcpy(O1, data.L2.Y, sizeof(F32)*N);

	//F32  bestRSS;

	//memcpy(O1, data.L2.Y, sizeof(F32)* 300);	
	//memcpy(O1, data.RSS_segSplitting, sizeof(F32)* 300);	
	//memcpy(O1, knot.XtY, sizeof(F32)* 17);
	//MEM.free_all(&MEM);
	//return;


}

#elif P_INTERFACE == 1

DllExport PyObject* pexFunction(PyObject* self, PyObject* args, PyObject* kwds) {

	int nargs = PyTuple_Size(args);
	int nkwds = PyDict_Size(kwds);

	if (nargs == 0) 	return Py_None;


	VOIDPTR   prhs[10] = {NULL,};
	int       nrhs = nargs; 
	for (int i = 0; i < nargs; i++) {
		prhs[i] = PyTuple_GetItem(args, i);
	}

	VOID_PTR ans= mainFunction(prhs, nrhs);

	return ans != NULL ? ans : IDE_NULL;
}


static PyMethodDef methods[] = {
  { "Rbeast",          &pexFunction, METH_VARARGS | METH_KEYWORDS, "Hello world function" },
  { "setClassObjects", &setClassObjects, METH_VARARGS , "Hello world function" },
 // { "hello", &hello, METH_VARARGS, "Hello world function" },
  { NULL, NULL, 0, NULL }
};

static struct PyModuleDef module_def = {
  PyModuleDef_HEAD_INIT, // always required
  "cbeast",               // module name
  "Testing module",      // description
  -1,                    // module size (-1 indicates we don't use this feature)
  methods,               // method table
};

PyMODINIT_FUNC PyInit_Rbeast() {
  
  PyObject *m= PyModule_Create(&module_def);

  /***********************************************************/
  BarObject_Type.tp_richcompare = PyBaseObject_Type.tp_richcompare;
  BarObject_Type.tp_hash         = PyBaseObject_Type.tp_hash;
  if (PyType_Ready(&BarObject_Type) < 0)
      return NULL;
   
  Py_INCREF(&BarObject_Type);
  if (PyModule_AddObject(m, "pyobject", (PyObject*)&BarObject_Type)) {
      Py_DECREF(&BarObject_Type);
      Py_DECREF(m);
      return NULL;
  }
  
//  printf(" tp_call %#x \n", ((PyObject*)&BarObject_Type)->ob_type->tp_call);

  /*
  if (PyType_Ready(&CustomType) < 0)
      return NULL; 

  Py_INCREF(&CustomType);
  if (PyModule_AddObject(m, "Custom", (PyObject*)&CustomType) < 0) {
      Py_DECREF(&CustomType);
      Py_DECREF(m);
      return NULL;
  }
  */
   
  //r_printf("Initialization done!\n");

  import_array();  // Load NumPy

  currentModule= m;
  return m;
}
#endif



#if R_INTERFACE==11111111111111

SEXP DllExport sbm2(SEXP Y, SEXP opt)
{
	if (!R_sbm_check_input(Y, opt))		 			return R_NilValue;

	char	missing[31];
	BEAST_OPTIONS beastOption;

	R_sbm_read_input(&beastOption, Y, opt, missing);

	if (!sbm_check_options(&beastOption, missing))	return R_NilValue;
	print_options(&beastOption);

	BEAST_RESULT	result;
	SEXP	ANS;
	PROTECT(ANS = R_allocate_output(&result, &beastOption));

	GLOBAL_OPTIONS = (BEAST_OPTIONS_PTR)&beastOption;
	GLOBAL_RESULT = (BEAST_RESULT_PTR)&result;
	sbm();
		UNPROTECT(1);
	return ANS;
}
#endif
 

#include "abc_000_warning.h"