#include "abc_000_warning.h"


// Need to be placed before time.h
#if ( defined(COMPILER_CLANG)||defined(COMPILER_GCC)||defined(COMPILER_SOLARIS) ) && !(defined(__APPLE__)||defined(__MACH__))
 
	#ifndef  _GNU_SOURCE
		#define _GNU_SOURCE
	#endif

    // stackoverflow.com/questions/40515557/compilation-error-on-clock-gettime-and-clock-monotonic
    // This turns on  "clock_getres(), clock_gettime(), clock_settime()".
	#ifndef  _POSIX_C_SOURCE
       #define _POSIX_C_SOURCE 199309L
	#endif

#endif

#include "abc_timer.h"

/*********************************************************************************************/
//
// Cross-platform sleep function 
//
/*********************************************************************************************/

// Google AI:
// sleep()	    Seconds	         <unistd.h>	POSIX.1-2001 (standard)
// usleep()	    Microseconds	 <unistd.h>	Obsolete in POSIX.1-2008
// nanosleep()	Nanoseconds	     <time.h>	Recommended POSIX (high-resolution)
// Nowadays, usleep() is obsolete: stackoverflow.com/questions/42861913/why-use-usleep-and-not-sleep

// stackoverflow.com/questions/264350/is-there-an-alternative-for-sleep-in-c/264378#264378
// stackoverflow.com/questions/10918206/cross-platform-sleep-function-for-c
#ifdef _WIN32
   #define WIN32_LEAN_AND_MEAN
    #include <windows.h>   // for Sleep,timeBeginPeriod,SetWaitableTimer, ...
#elif _POSIX_C_SOURCE >= 199309L || defined (__MUSL__)  // MUSL library
    #include <time.h>      // for  timespec, nanosleep
#else
    #include <unistd.h>    // for usleep
    #include <sys/time.h>  // for timeval
#endif


void Sleep_ms(int milliseconds) {
	
    #ifdef WIN32
	
	   // Goolge AI: Sleep has low precision (around 10-15 ms resolution
	  if   (milliseconds > 100){
		   Sleep(milliseconds);
	  } 
	  /* 
	  // Need to undefine "WIN32_LEAN_AND_MEAN" to make timeBeginPeriod visible
	  // also need to specify the "Winmm.lib" library for timeBeginPeriod in the linking 
	  else if (milliseconds > 10){
	  // https://learn.microsoft.com/en-us/windows/win32/api/timeapi/nf-timeapi-timebeginperiod	  
	  //   1-ms resolution
		    timeBeginPeriod(1);
			Sleep(milliseconds);		
			timeEndPeriod(1);
   	  }
	  */
	  else {
       // https://stackoverflow.com/questions/62273909/sleeping-using-windows-api-createwaitabletimer-executes-differently-on-xeon-comp
	   // https://learn.microsoft.com/en-us/windows/win32/sync/using-waitable-timer-objects
	      LARGE_INTEGER ft;		  	  
		  I64 usec = milliseconds*1000;		  
  		  ft.QuadPart = -(10*usec); // Convert to 100 nanosecond interval, negative value indicates relative time

		  HANDLE timer = CreateWaitableTimer(NULL, TRUE, NULL);
		  if (timer) {			
			SetWaitableTimer(timer, &ft, 0, NULL, NULL, 0);
			WaitForSingleObject(timer, INFINITE);
			CloseHandle(timer);		  		    
		  }
	  }
	  
    #elif _POSIX_C_SOURCE >= 199309L || defined (__MUSL__)  // MUSL library
	
	// https://gist.github.com/scivision/dbbbf33c2faf5a16f31fd6d144adc314
	// The POSIX nanosleep() function is available in Linux, macOS, BSD, Unix, etc. MinGW / MSYS2 implements 
	// nanosleep for certain platforms -- x86_64, but not ARM currently. Cygwin implements nanosleep() for x86_64 
	// and ARM. MSVC Visual Studio currently does not implement nanosleep().
        struct timespec ts;
        ts.tv_sec  =  milliseconds / 1000;
        ts.tv_nsec = (milliseconds % 1000) * 1000000;
        nanosleep( &ts, NULL );
		
    #else	
        /*
        if (milliseconds >= 1000)
            sleep(milliseconds / 1000);
        usleep((milliseconds % 1000) * 1000);
        */
        struct timeval tv; 
        tv.tv_sec  = milliseconds / 1000;
        tv.tv_usec = milliseconds % 1000 * 1000; 
        select(0, NULL, NULL, NULL, &tv);
    #endif
}



/*********************************************************************************************/
//
// Cross-platform rdtsc function 
//
/*********************************************************************************************/


// stackoverflow.com/questions/9887839/how-to-count-clock-cycles-with-rdtsc-in-gcc-x86
#ifdef COMPILER_MSVC

    #include <intrin.h>    // __rdtsc
	
#elif defined(COMPILER_SOLARIS)

    //https://docs.oracle.com/cd/E24457_01/html/E21991/gliwk.html
	// for some reason, mmintrinc.h deson't contain rdtsc
    //#include <mmintrin.h> 
    
   //https://github.com/BackupTheBerlios/iometer-svn/blob/cd9ca025b44a0e80015a040b838d6dc7d486642a/tags/initial/iometer/src/rdtsc.c
   #include <sys/time.h>
   static INLINE  unsigned long long __rdtsc(void) {        
        return  gethrtime(); //https://community.oracle.com/tech/developers/discussion/2490026/about-gethrtime-function
   }
 
    /*
    //https://github.com/michealbrownm/phantom/blob/670d80162ac2e05ebcaffa7a61365068d84a3672/3rd/ed25519-donna/test-ticks.h
    #elif defined(COMPILER_GCC)
        uint32_t lo, hi;
        __asm__ __volatile__("rdtsc" : "=a" (lo), "=d" (hi));
        return ((uint64_t)lo | ((uint64_t)hi << 32));
    #else
    #elif defined(CPU_SPARC) && !defined(OS_OPENBSD)
		uint64_t t;
		__asm__ __volatile__("rd %%tick, %0" : "=r" (t));
		return t;
    */
	
#elif defined(cpu_ARM64)
    /*
    https://stackoverflow.com/questions/44349757/rdtsc-equivalent-for-ios
    There appears to be a clang builtin for this, __builtin_readcyclecounter(). 
    On x86, it compiles to rdtsc; on ARM64, it compiles to mrs x0, PMCCNTR_EL0. 
    On ARM it always returns 0 though. ARM64 is good enough for me though, since
    32-bit support was dropped in iOS 11.
    */
	
  static INLINE  U64 __rdtsc(void)   {
	  return  __builtin_readcyclecounter(); 
  }

  //https://stackoverflow.com/questions/12631856/difference-between-rdtscp-rdtsc-memory-and-cpuid-rdtsc
  //https://stackoverflow.com/questions/40454157/is-there-an-equivalent-instruction-to-rdtsc-in-arm
   static INLINE U64 rdtsc(void) {
        U64 val;
        /*
         * According to ARM DDI 0487F.c, from Armv8.0 to Armv8.5 inclusive, the
         * system counter is at least 56 bits wide; from Armv8.6, the counter
         * must be 64 bits wide.  So the system counter could be less than 64
         * bits wide and it is attributed with the flag 'cap_user_time_short'
         * is true.
         */
        asm volatile("mrs %0, cntvct_el0" : "=r" (val));
        return val;
    }

#else
	
    #include <x86intrin.h> // __rdtsc
	
#endif

static INLINE unsigned long long readTSC(void) {
    // _mm_lfence();  // optionally wait for earlier insns to retire before reading the clock
    return __rdtsc();
    // _mm_lfence();  // optionally block later instructions until rdtsc retires
}

static unsigned long long ReadTSC_InternalTime_TIC;
unsigned long long  tic(void) { return ReadTSC_InternalTime_TIC = readTSC();}
unsigned long long  toc(void) {	return readTSC()- ReadTSC_InternalTime_TIC;}

/*********************************************************************************************/
//
// Define a set of timer functions
//
/*********************************************************************************************/

#ifdef COMPILER_MSVC
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>        // QueryPerformanceFrequency, QueryPerformanceCounter
#elif ( defined(COMPILER_CLANG)||defined(COMPILER_GCC)|| defined(COMPILER_SOLARIS) ) && !(defined(__APPLE__)||defined(__MACH__))
	#include "time.h"           //  define CLOCK_REAL_TIME	
#elif defined( __MACH__)
    #include <mach/mach_time.h>  // mach_timebase_info, mach_absolute_time
#endif


#if defined(COMPILER_MSVC)	
	static LARGE_INTEGER T0;
#elif ( defined(COMPILER_CLANG)||defined(COMPILER_GCC)|| defined(COMPILER_SOLARIS) ) && !(defined(__APPLE__)||defined(__MACH__))
	static struct timespec T0;
#elif defined(__MACH__)
	static U64 T0;
#endif

static F64 fctConversionToSecond   = 0.0;
static F64 elapsedTimeAtBreakpoint;

 // Get conversion factors from counts to seceonds
static F64 __Get_fctConversionToSecond (void) {
	
	F64 _fctConversionToSecond;
#if defined(COMPILER_MSVC)
	LARGE_INTEGER Frequency;
	QueryPerformanceFrequency(&Frequency); 
	_fctConversionToSecond = 1. / (double)Frequency.QuadPart;
#elif defined(__MACH__)
	mach_timebase_info_data_t timebase;
	mach_timebase_info(&timebase);
	F64 conversionTonanosecs;
	conversionTonanosecs   = timebase.numer / timebase.denom;
	_fctConversionToSecond  = conversionTonanosecs/1e9;
#else
	_fctConversionToSecond=1.0;
#endif
   return _fctConversionToSecond;
}


void Timer_Start() {
	
	if ( fctConversionToSecond ==0.01){
		fctConversionToSecond = __Get_fctConversionToSecond();
	}
	
#if defined(COMPILER_MSVC)  
    //http: //forums.codeguru.com/showthread.php?261134-QueryPerformanceCounter-to-milliseconds
	QueryPerformanceCounter(&T0);	
#elif ( defined(COMPILER_CLANG)||defined(COMPILER_GCC)||defined(COMPILER_SOLARIS) ) && !(defined(__APPLE__)||defined(__MACH__))
    //http: //linuxmogeb.blogspot.com/2013/10/how-does-clockgettime-work.html
	//CLOCK_MONOTONIC vs CLOCK_REALTIME 
	clock_gettime(CLOCK_MONOTONIC, &T0);	
#elif defined(__MACH__) 
    // https: //stackoverflow.com/questions/5167269/clock-gettime-alternative-in-mac-os-x
	T0 = mach_absolute_time();	
#endif

}

F64 Timer_ElapsedSecond() {	
	F64 elapsedTime;
#if defined(COMPILER_MSVC)  
//http: //forums.codeguru.com/showthread.php?261134-QueryPerformanceCounter-to-milliseconds
	LARGE_INTEGER T;
	QueryPerformanceCounter(&T);	
	elapsedTime = (double)(T.QuadPart-T0.QuadPart) * fctConversionToSecond;
#elif ( defined(COMPILER_CLANG)||defined(COMPILER_GCC)||defined(COMPILER_SOLARIS) ) && !(defined(__APPLE__)||defined(__MACH__))
//CLOCK_MONOTONIC vs CLOCK_REALTIME //http: //linuxmogeb.blogspot.com/2013/10/how-does-clockgettime-work.html
	struct timespec T;
	clock_gettime(CLOCK_MONOTONIC, &T);
	elapsedTime = (double) (T.tv_sec - T0.tv_sec) + (double)(T.tv_nsec-T0.tv_nsec)/1000000000LL;
#elif defined(__MACH__) 
// https: //stackoverflow.com/questions/5167269/clock-gettime-alternative-in-mac-os-x
	U64 T       =  mach_absolute_time();
	elapsedTime =  (T - T0)*fctConversionToSecond;
#endif
	return elapsedTime;
}

void Timer_SetBreakPt() {
	elapsedTimeAtBreakpoint = Timer_ElapsedSecond();
}

F64 Timer_ElapsedSinceBeakPt() {
	return Timer_ElapsedSecond() - elapsedTimeAtBreakpoint;
}

U64 TimerGetTickCount() {
	U64 tick;
#if	defined(COMPILER_MSVC) 
	tick = GetTickCount64();
#elif ( defined(COMPILER_CLANG)||defined(COMPILER_GCC)||defined(COMPILER_SOLARIS) ) && !(defined(__APPLE__)||defined(__MACH__))
	// https:// stackoverflow.com/questions/17432730/precedence-of-over
	struct timespec T;
	clock_gettime(CLOCK_REALTIME, &T);
	tick =T.tv_sec * 1000000000LL + T.tv_nsec;
#elif defined(__MACH__) //https: //stackoverflow.com/questions/5167269/clock-gettime-alternative-in-mac-os-x
	tick= mach_absolute_time();
#endif
	return tick;
}


#include "abc_000_warning.h"