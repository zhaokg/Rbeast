#include "abc_000_warning.h"
#include "abc_000_macro.h"

//https://stackoverflow.com/questions/150355/programmatically-find-the-number-of-cores-on-a-machine
#if defined(_WIN32) || defined(WIN64_OS)
    //#define WIN32_LEAN_AND_MEAN
	//#include <windows.h>
    #include "abc_pthread.h"
#elif defined(MAC_OS)
	#include <sys/param.h>
	#include <sys/sysctl.h>
#elif defined(LINUX_OS) || defined(SOLARIS_OS)
//https://docs.oracle.com/cd/E36784_01/html/E36874/sysconf-3c.html
	#include <unistd.h>
#endif

#include <inttypes.h>
#include "abc_ide_util.h" // for printf only

int CountSetBits32(uint32_t x)  {
    //https://en.wikipedia.org/wiki/Hamming_weight
    //https://stackoverflow.com/questions/109023/how-to-count-the-number-of-set-bits-in-a-32-bit-integer
    x = x - ((x >> 1) & 0x55555555);
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
    x = (x + (x >> 4)) & 0x0F0F0F0F;
    x = x + (x >> 8);
    x = x + (x >> 16);
    return x & 0x0000003F;
}

int CountSetBits64(uint64_t x) {
    //uint64_t is an unsigned 64-bit integer variable type (defined in C99 version of C language)
    static const uint64_t m1 = 0x5555555555555555; //binary: 0101...
    static const uint64_t m2 = 0x3333333333333333; //binary: 00110011..
    static const uint64_t m4 = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
    static const uint64_t m8 = 0x00ff00ff00ff00ff; //binary:  8 zeros,  8 ones ...
    static const uint64_t m16 = 0x0000ffff0000ffff; //binary: 16 zeros, 16 ones ...
    static const uint64_t m32 = 0x00000000ffffffff; //binary: 32 zeros, 32 ones
    static const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
    x += x >> 8;  //put count of each 16 bits into their lowest 8 bits
    x += x >> 16;  //put count of each 32 bits into their lowest 8 bits
    x += x >> 32;  //put count of each 64 bits into their lowest 8 bits
    return x & 0x7f;
}

int GetNumCores()  {
    //https://stackoverflow.com/questions/150355/programmatically-find-the-number-of-cores-on-a-machine
#if defined(_WIN32) || defined(WIN64_OS)    
    uint32_t count;
	count = GetCPUInfo(); 
#elif defined(__MACH__)
    int nm[2];
    size_t   len = 4;
    uint32_t count;

    nm[0] = CTL_HW; 
    nm[1] = HW_AVAILCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);

    if(count < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &count, &len, NULL, 0);
        if(count < 1) { count = 1; }
    }
#else
    uint32_t count;
    count = sysconf(_SC_NPROCESSORS_ONLN);
#endif
	return count >= 1 ? count : 1L;
}


#if defined(_WIN32) || defined(WIN64_OS)

typedef struct __CPUINFO {
    //docs.microsoft.com/en-us/windows/win32/procthread/what-s-new-in-processes-and-threads
    DWORD (WINAPI* GetActiveProcessorGroupCount)     ();
    DWORD (WINAPI* GetActiveProcessorCount)          (WORD GroupNumber);
    BOOL  (WINAPI* GetLogicalProcessorInformationEx) (
            LOGICAL_PROCESSOR_RELATIONSHIP           RelationshipType,        
            PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX Buffer,
            PDWORD                                   ReturnedLength );
    HANDLE (WINAPI*  CreateRemoteThreadEx) (
            HANDLE                       hProcess,
            LPSECURITY_ATTRIBUTES        lpThreadAttributes,
            SIZE_T                       dwStackSize,
            LPTHREAD_START_ROUTINE       lpStartAddress,
            LPVOID                       lpParameter,
            DWORD                        dwCreationFlags,
            LPPROC_THREAD_ATTRIBUTE_LIST lpAttributeList,
            LPDWORD                      lpThreadId   );
    BOOL (WINAPI* GetThreadGroupAffinity) (
            HANDLE          hThread, 
            PGROUP_AFFINITY GroupAffinity );
    void  (WINAPI*  GetCurrentProcessorNumberEx) ( PPROCESSOR_NUMBER ProcNumber );
    DWORD (WINAPI* GetCurrentProcessorNumber) (VOID ); //since Windows vista
    BOOL  (WINAPI*  GetNumaHighestNodeNumber)(  PULONG HighestNodeNumber   );
    BOOL  (WINAPI*  GetNumaProcessorNodeEx)  (
                PPROCESSOR_NUMBER Processor,
                PUSHORT           NodeNumber   );
 
    BOOL(WINAPI* GetNumaAvailableMemoryNode) (
                UCHAR      Node,
                PULONGLONG AvailableBytes   ); //available since Win XP

    char isInitilized;
} CPUFUCINFO;

typedef struct __CPUINFO1 {
    DWORD logicalProcessorCount;
    DWORD numaNodeCount;    
    DWORD processorCoreCount;
    DWORD processorL1CacheCount;
    DWORD processorL2CacheCount;
    DWORD processorL3CacheCount ;
    DWORD processorPackageCount;
    DWORD processorGroupCount;
    DWORD coreCountPerGrp[10];
    DWORD currentGroup;
    DWORD currentCoreNumber;
    uint64_t currentThreadAffinity;
    char cpuGroup[256];
    char numCPUCoresToUseber[256];
    //uint64_t currentGroupAffinity;

} CPUINFO;

static CPUFUCINFO cpuFunc = {0,};
static CPUINFO    cpuInfo = {0,};

static void InitCPUFuncs( ) {
    if (cpuFunc.isInitilized)  {
        return;
    }
    HMODULE kerHandle  = GetModuleHandle(TEXT("kernel32"));
    cpuFunc.GetActiveProcessorGroupCount     =  GetProcAddress(kerHandle, "GetActiveProcessorGroupCount");
    cpuFunc.GetActiveProcessorCount          =  GetProcAddress(kerHandle, "GetActiveProcessorCount");
    cpuFunc.GetLogicalProcessorInformationEx = GetProcAddress(kerHandle, "GetLogicalProcessorInformationEx");
    cpuFunc.CreateRemoteThreadEx             = GetProcAddress(kerHandle, "CreateRemoteThreadEx");
    cpuFunc.GetThreadGroupAffinity           = GetProcAddress(kerHandle, "GetThreadGroupAffinity");
    cpuFunc.GetCurrentProcessorNumberEx      = GetProcAddress(kerHandle, "GetCurrentProcessorNumberEx");
    cpuFunc.GetCurrentProcessorNumber        = GetProcAddress(kerHandle, "GetCurrentProcessorNumber");
    cpuFunc.GetNumaHighestNodeNumber         = GetProcAddress(kerHandle, "GetNumaHighestNodeNumber");
    cpuFunc.GetNumaProcessorNodeEx           = GetProcAddress(kerHandle, "GetNumaProcessorNodeEx");
    cpuFunc.GetNumaAvailableMemoryNode       = GetProcAddress(kerHandle, "GetNumaAvailableMemoryNode");
    cpuFunc.isInitilized = 1;
}

static int  GetCoreNumbers_WIN32() {

    cpuInfo = (CPUINFO){ 0, };

    SYSTEM_INFO   sysinfo;
    GetSystemInfo(&sysinfo);   
    cpuInfo.logicalProcessorCount = sysinfo.dwNumberOfProcessors;
    
    cpuInfo.processorGroupCount = 1;
    uint64_t currentGroupAffinity;
    return cpuInfo.logicalProcessorCount;
}

static int GetCoreNumbers_WIN7V1() {
    cpuInfo = (CPUINFO){ 0, };

    InitCPUFuncs();
    if (NULL == cpuFunc.GetActiveProcessorGroupCount)   {
        return -1;
    }

    WORD procGroupNum = cpuFunc.GetActiveProcessorGroupCount();

    procGroupNum = procGroupNum > 10 ? 10 : procGroupNum;
    int numCore = 0;
    for (DWORD i = 0; i < procGroupNum; i++)   {
        int count= cpuFunc.GetActiveProcessorCount(i);
        numCore += count;
        cpuInfo.coreCountPerGrp[i] = count;
    }
    
    cpuInfo.processorGroupCount   = procGroupNum;
    cpuInfo.logicalProcessorCount = numCore;

    ULONG numaCount;
    cpuFunc.GetNumaHighestNodeNumber(&numaCount);
    cpuInfo.numaNodeCount = numaCount+1L;
    return cpuInfo.logicalProcessorCount; 
}

//https://docs.microsoft.com/en-us/windows/win32/api/sysinfoapi/nf-sysinfoapi-getlogicalprocessorinformation
static int GetCoreNumbers_WINXP()
{         
    cpuInfo = (CPUINFO){ 0, };

    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION buffer = NULL;
   

  
    PCACHE_DESCRIPTOR Cache;

    DWORD returnLength = 0;
    BOOL  done         = FALSE;
    while (!done)  {

        DWORD flagOK = GetLogicalProcessorInformation(buffer, &returnLength);
        if (FALSE == flagOK)  {
            if (GetLastError() == ERROR_INSUFFICIENT_BUFFER)   {
                if (buffer) {
                    free(buffer);
                }                
                buffer = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION) malloc(returnLength);
                if (NULL == buffer) { //_tprintf(TEXT("\nError: Allocation failure\n"));
                    return (-1);
                }
            }  else   {   //_tprintf(TEXT("\nError %d\n"), GetLastError());
                return (-1);
            }
        } else {
            done = TRUE;
        }
    }

    
    
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION ptr =  buffer;
    DWORD byteOffset = 0;
    while (byteOffset + sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION) <= returnLength)     {

        switch (ptr->Relationship)
        {
        case RelationNumaNode:    
            cpuInfo.numaNodeCount++;     // Non-NUMA systems report a single record of this type.
            break;

        case RelationProcessorCore:
            cpuInfo.processorCoreCount++;
            // A hyperthreaded core supplies more than one logical processor.
            cpuInfo.logicalProcessorCount += CountSetBits64(ptr->ProcessorMask);
            break;

        case RelationCache:
            // Cache data is in ptr->Cache, one CACHE_DESCRIPTOR structure for each cache. 
            Cache = &ptr->Cache;
            if (Cache->Level == 1) {
                cpuInfo.processorL1CacheCount++;
            } else if (Cache->Level == 2)   {
                cpuInfo.processorL2CacheCount++;
            } else if (Cache->Level == 3)   {
                cpuInfo.processorL3CacheCount++;
            }
            break;

        case RelationProcessorPackage:
            // Logical processors share a physical package.
            cpuInfo.processorPackageCount++;
            break;

        default:
           // _tprintf(TEXT("\nError: Unsupported LOGICAL_PROCESSOR_RELATIONSHIP value.\n"));
            break;
        }
        byteOffset += sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
        ptr++;
    }

    /*
    _tprintf(TEXT("\nGetLogicalProcessorInformation results:\n"));
    _tprintf(TEXT("Number of NUMA nodes: %d\n"),      numaNodeCount);
    _tprintf(TEXT("Number of physical processor packages: %d\n"),      processorPackageCount);
    _tprintf(TEXT("Number of processor cores: %d\n"),        processorCoreCount);
    _tprintf(TEXT("Number of logical processors: %d\n"),     logicalProcessorCount);
    _tprintf(TEXT("Number of processor L1/L2/L3 caches: %d/%d/%d\n"),     processorL1CacheCount,
        processorL2CacheCount,      processorL3CacheCount);
    */

    free(buffer);

    return cpuInfo.logicalProcessorCount;
}

static int GetCoreNumbers_WIN7V2() {

    #ifdef WIN64_OS

    cpuInfo = (CPUINFO){0, };
    

    InitCPUFuncs();

    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX buffer = NULL;
    DWORD returnLength = 0;
    BOOL  done         = FALSE;
    while (!done) {
         DWORD fOK = cpuFunc.GetLogicalProcessorInformationEx(
                             RelationAll, (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX) buffer, &returnLength);         
         if (FALSE == fOK) {
            if (GetLastError() == ERROR_INSUFFICIENT_BUFFER)  {
               if (buffer) 
                   free(buffer);
               buffer = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX)malloc(returnLength);

               if (NULL == buffer) {  // _tprintf(TEXT("\nError: Allocation failure\n"));
                    return (-1);
               }
            } else {       //_tprintf(TEXT("\nError %d\n"), GetLastError());
                return (-1);
            }
         }  else {
                done = TRUE;
         }

    }

     PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX  ptr = buffer;
 
     DWORD byteOffset = 0;
     while (byteOffset < returnLength) {

            switch (ptr->Relationship)
            {
            case RelationNumaNode:
                // Non-NUMA systems report a single record of this type.
                cpuInfo.numaNodeCount++;
                break;
                
            case RelationProcessorCore:
                cpuInfo.processorCoreCount++;
                /*
                printf("************\n");
                printf("%u\n",  ptr->Processor.Flags);
                for (int i = 0; i < ptr->Processor.GroupCount; i++)  {

                    printf("%u\n",   ptr->Processor.GroupMask[i].Group);
                    printf("%08x\n", ptr->Processor.GroupMask[i].Mask);
                }

                printf("------\n");
                */
                //A hyperthreaded core supplies more than one logical processor.
                cpuInfo.logicalProcessorCount += ptr->Processor.Flags + 1L;
                break;

            case RelationCache:
                // Cache data is in ptr->Cache, one CACHE_DESCRIPTOR structure for each cache. 

                break;

            case RelationProcessorPackage:
                // Logical processors share a physical package.
                cpuInfo.processorPackageCount++;
                /*
                printf("************\n");
                printf("%u\n", ptr->Processor.Flags);
                for (int i = 0; i < ptr->Processor.GroupCount; i++) {
                    printf("%u\n", ptr->Processor.GroupMask[i].Group);
                    printf("%08x\n", ptr->Processor.GroupMask[i].Mask);
                }
                printf("------\n");
                */
                break;
            case RelationGroup:
                /*
                printf("************\n");
                printf("%u\n", ptr->Group.ActiveGroupCount);
                for (int i = 0; i < ptr->Group.ActiveGroupCount; i++) {
                    printf("%u\n", ptr->Group.GroupInfo[i].ActiveProcessorCount);
                    printf("%08x\n", ptr->Group.GroupInfo[i].ActiveProcessorMask);
                }
                printf("------\n");
                */
                {
                    int activeCount = ptr->Group.ActiveGroupCount;
                    for (int i = 0; i < activeCount; i++) {
                          cpuInfo.coreCountPerGrp[cpuInfo.processorGroupCount + i] = ptr->Group.GroupInfo[i].ActiveProcessorCount;
                     }
                    cpuInfo.processorGroupCount += activeCount;
                }
                break;
            default:
                //_tprintf(TEXT("\nError: Unsupported LOGICAL_PROCESSOR_RELATIONSHIP value.\n"));
                break;
            }

            byteOffset += ptr->Size;
            ptr         = (char*)ptr + ptr->Size;
        }

     
     PROCESSOR_NUMBER procNum;
     cpuFunc.GetCurrentProcessorNumberEx(&procNum);
     cpuInfo.currentGroup = procNum.Group;
     cpuInfo.currentCoreNumber = procNum.Number;

     GROUP_AFFINITY grpAffinity;
     cpuFunc.GetThreadGroupAffinity( GetCurrentThread(), &grpAffinity);
     cpuInfo.currentThreadAffinity = grpAffinity.Mask;
     
     return cpuInfo.logicalProcessorCount;

  
    #endif

     return 0;
}

static void RankCPU() {

#ifdef WIN64_OS


    int nGrp   = cpuInfo.processorGroupCount;
    int curGrp = cpuInfo.currentGroup;
    int curNo  = cpuInfo.currentCoreNumber;

    int coreCounts = cpuInfo.coreCountPerGrp[curGrp];

    int idx = 0;
    for (int i = 0; i < coreCounts; i++) {
        if (i != curNo) {
            cpuInfo.cpuGroup[idx]      = curGrp;
            cpuInfo.numCPUCoresToUseber[idx] = i;
            idx++;
        }
    }
    idx = coreCounts - 1;
    cpuInfo.cpuGroup[idx]      = curGrp;
    cpuInfo.numCPUCoresToUseber[idx] = curNo;


    // Sort through the other processor groups
    idx = coreCounts;
    for (int grp = 0; grp < nGrp; grp++) {
        if (grp == curGrp)
            continue;        

        coreCounts = cpuInfo.coreCountPerGrp[grp];
        for (int i = 0; i < coreCounts; i++) {
            cpuInfo.cpuGroup[idx %256]      = grp;
            cpuInfo.numCPUCoresToUseber[idx %256] = i;
            idx++;
        }

    }
#endif
}

int GetCPUInfo() {

    InitCPUFuncs();

    int nCores;
    if (cpuFunc.GetLogicalProcessorInformationEx != NULL) {    
        nCores=GetCoreNumbers_WIN7V2();
        RankCPU();
    }  else  {
        nCores = GetCoreNumbers_WIN32();
    }
    PrintCPUInfo();
    return nCores;
}


void PrintCPUInfo() {

    r_printf("\nCPU Information:\n");
    r_printf(" - Number of NUMA nodes: %d\n", cpuInfo.numaNodeCount);
    r_printf((" - Number of physical processors (sockets): %d\n"), cpuInfo.processorPackageCount);
    r_printf((" - Number of processor cores: %d\n"), cpuInfo.processorCoreCount);
    r_printf((" - Number of logical processors: %d\n"), cpuInfo.logicalProcessorCount);
    r_printf(" - Number of processor groups: %d\n", cpuInfo.processorGroupCount);
    for (int i = 0; i < cpuInfo.processorGroupCount; i++) {
        r_printf(" -- Processor group #%d: %d cores\n", i, cpuInfo.coreCountPerGrp[i]);
    }
    r_printf((" - Number of processor L1/L2/L3 caches: %d/%d/%d\n"), cpuInfo.processorL1CacheCount,
        cpuInfo.processorL2CacheCount, cpuInfo.processorL3CacheCount);
    r_printf(" - Group ID of current thread: %d\n", cpuInfo.currentGroup);
    r_printf(" - Core ID of current thread: %d\n", cpuInfo.currentCoreNumber);
    r_printf(" - CPU affinity mask of current thread: %#x\n", cpuInfo.currentThreadAffinity);

}

 void CPU_ZERO(cpu_set_t* cpus) {
    //cpus->cpu is a PROCESSRO_NUMBER; the reserved field must be sset to 0s; otherwise
    //CreateRemoteTHreadEX fails.
    //PROCESSOR_NUMBER procNum = {0,};
    memset(cpus, 0, sizeof(cpu_set_t));
}
 void CPU_SET(int i, cpu_set_t* cpus) {     
    #ifdef WIN64_OS
    cpus->ProcNumber.Group  = cpuInfo.cpuGroup[i];
    cpus->ProcNumber.Number = cpuInfo.numCPUCoresToUseber[i];
    #endif
}


 // pthread_create has a conflict with the pthread lib if provided/
 // instead of exploring options with weak symbol, I make it static here
int  pthread_create0(pthread_t* tid,  const pthread_attr_t* attr, void* (*start) (void*), void* arg)
{

#ifdef WIN64_OS
    if (cpuFunc.CreateRemoteThreadEx == NULL) {
        r_printf("the CreateRemoteThreadEx funnction is not detected!\n");
        *tid = CreateThread(NULL, NULL, (LPTHREAD_START_ROUTINE) start, arg, 0, 0);
    }  else {
        //https://docs.microsoft.com/en-us/windows/win32/procthread/thread-stack-size
        *tid = cpuFunc.CreateRemoteThreadEx(
            GetCurrentProcess(),
            (LPSECURITY_ATTRIBUTES)NULL,
            (SIZE_T)attr->dwStackSize, // if 0, use the default stack size
            (LPTHREAD_START_ROUTINE)start,
            (LPVOID)arg,
            (DWORD)0,
            (LPPROC_THREAD_ATTRIBUTE_LIST)attr->lpAttributeList,
            (LPDWORD)NULL
        );

        //https://microsoft.public.win32.programmer.kernel.narkive.com/7QcTnq9M/createthread-fails-with-error-not-enough-memory
    }
#else
    * tid = CreateThread(NULL, attr->dwStackSize, (LPTHREAD_START_ROUTINE)start, arg, 0, 0);
#endif

    return 0;
}


 

#endif


#if defined(_WIN32) || defined(WIN64_OS)

    #if _WIN32_WINNT >= 0x0602
       //The function below is to get hte stackszie of the current thread.  
        int get_thread_stacksize() {
            ULONG_PTR lowLimit;
            ULONG_PTR highLimit; 
            GetCurrentThreadStackLimits(&lowLimit, &highLimit);
            return (highLimit - lowLimit);
        }
    #else
       // https://www.xsprogram.com/content/can-i-get-the-limits-of-the-stack-in-c-c.html
        int get_thread_stacksize() {
             NT_TIB* tib = (NT_TIB*)NtCurrentTeb();
             LPVOID   StackLimit = tib->StackLimit;
             LPVOID   StackBase  = tib->StackBase;
             return  StackLimit;
        }
    #endif

    // On Windows, Seems like that there is no functions to get the default stack size;
     //, except parseing the PE's header: https://stackoverflow.com/questions/42420727/default-stack-size
    // A sloppy way to get the default stack size: create a thread and run GetCurrentTrheadLimits
    static void * __getstacksize(void * arg) {  *((size_t*) arg)=get_thread_stacksize();  return 0; }
    int pthread_attr_getstacksize_win32(pthread_attr_t* attr, size_t* stacksize) {
        static int default_stacksize=0;
        if (attr->dwStackSize > 0) {
            *stacksize = attr->dwStackSize;            
        } else{
            if (default_stacksize == 0) {
                pthread_t tid =  CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)__getstacksize, stacksize, 0, 0);
                pthread_join(tid, NULL);
                default_stacksize = *stacksize;
            }
             
            *stacksize = default_stacksize;
        }       

        return 0;        
    }
    
#elif defined(LINUX_OS) 

//https://docs.oracle.com/cd/E88353_01/html/E37843/pthread-getattr-np-3c.html
int get_thread_stacksize() {

    pthread_attr_t attr;
    size_t   stksize;
    void*    stkaddr;
    //pthread_getattr_np not available in Solaris as reported by R-hub
    (void)pthread_getattr_np(pthread_self(), &attr);
    (void)pthread_attr_getstack(&attr, &stkaddr, &stksize);
    (void)pthread_attr_destroy(&attr);
    return stksize;
}

#else // For non-windows systems, define a dummy function
//pthread_getattr_np ; But not portable to Mac OS
//pthread_self() 
int get_thread_stacksize() {

    return  0;
    }
#endif
#include "abc_000_warning.h"