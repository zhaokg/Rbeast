#pragma once

#include "abc_000_macro.h"
#include "abc_datatype.h"

#ifdef __cplusplus
extern "C" {
#endif



#ifdef COMPILER_MSVC
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>          //QueryPerformanceFrequency QueryPerformanceCounter
#elif defined( __MACH__)
    #include <mach/mach_time.h>
#endif


//stackoverflow.com/questions/264350/is-there-an-alternative-for-sleep-in-c/264378#264378
#ifdef _WIN32
    #include <windows.h>
#elif _POSIX_C_SOURCE >= 199309L || defined (__MUSL__)  // MUSL library
    #include <time.h>  // for  timespec
    #include <time.h>    // for nanosleep
#else
    #include <unistd.h> // for usleep
    #include <sys/time.h>  // for timeval
#endif
 
//https://stackoverflow.com/questions/10918206/cross-platform-sleep-function-for-c
static INLINE void Sleep_ms(int milliseconds) {
    // cross-platform sleep function
    #ifdef WIN32
        Sleep(milliseconds);
    #elif _POSIX_C_SOURCE >= 199309L || defined (__MUSL__)  // MUSL library
        struct timespec ts;
        ts.tv_sec  = milliseconds / 1000;
        ts.tv_nsec = (milliseconds % 1000) * 1000000;
        nanosleep(&ts, NULL);
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


//warning: function declaration isn't a prototype [-Wstrict-prototypes]
//https://stackoverflow.com/questions/42125/warning-error-function-declaration-isnt-a-prototype
extern void InitTimerFunc(void);
extern void StartTimer(void);
extern F64  GetElapsedSecondsSinceStart(void);
extern void SetBreakPointForStartedTimer(void);
extern F64  GetElaspedTimeFromBreakPoint(void);
extern U64  TimerGetTickCount(void);


// stackoverflow.com/questions/9887839/how-to-count-clock-cycles-with-rdtsc-in-gcc-x86
#ifdef COMPILER_MSVC
    #include <intrin.h>    // __rdtsc
#elif defined(COMPILER_SOLARIS)
    //https://docs.oracle.com/cd/E24457_01/html/E21991/gliwk.html
    //#include <mmintrin.h> // __rdtsc for some reason, mmintrinc.h deson't contain rdtsc
    
   //https://github.com/BackupTheBerlios/iometer-svn/blob/cd9ca025b44a0e80015a040b838d6dc7d486642a/tags/initial/iometer/src/rdtsc.c
   #include <sys/time.h>
    static INLINE   unsigned long long __rdtsc(void) {
        //https://community.oracle.com/tech/developers/discussion/2490026/about-gethrtime-function
        return  gethrtime();
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
  static INLINE  U64 __rdtsc(void)   {  return  __builtin_readcyclecounter();  }

  //https://stackoverflow.com/questions/12631856/difference-between-rdtscp-rdtsc-memory-and-cpuid-rdtsc
  //https://stackoverflow.com/questions/40454157/is-there-an-equivalent-instruction-to-rdtsc-in-arm
    static INLINE U64 rdtsc(void)     {
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

extern void               tic(void);
extern unsigned long long toc(void);



#ifdef __cplusplus
}
#endif