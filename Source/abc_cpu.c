#include <string.h>

#include "abc_000_warning.h"

#include "abc_common.h"
#include "abc_ide_util.h" // r_printf
#include "abc_cpu.h"

#include <stdio.h>

//https:// github.com/Mysticial/FeatureDetector
//https:// stackoverflow.com/questions/6121792/how-to-check-if-a-cpu-supports-the-sse3-instruction-set/22521619#22521619
/*
Mysticial's answer is a bit dangerous -- it explains how to detect CPU support but not OS support. You need to use _xgetbv to check
whether the OS has enabled the required CPU extended state. See here for another source. Even gcc
has made the same mistake. The meat of the code is:
*/

//https:// insufficientlycomplicated.wordpress.com/2011/11/07/detecting-intel-advanced-vector-extensions-avx-in-visual-studio/
/*
    // Checking for AVX requires 3 things:
    // 1) CPUID indicates that the OS uses XSAVE and XRSTORE
    //     instructions (allowing saving YMM registers on context
    //     switch)
    // 2) CPUID indicates support for AVX
    // 3) XGETBV indicates the AVX registers will be saved and
    //     restored on context switch
    //
    // Note that XGETBV is only available on 686 or later CPUs, so
    // the instruction needs to be conditionally run.
*/

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
    //#if _WIN32 

    /*********************************************/
    //        WINDOWS
     /*********************************************/
    #if defined(COMPILER_MSVC) 


        #define WIN32_LEAN_AND_MEAN
        #include <Windows.h>
        #include <intrin.h>

        #define _XCR_XFEATURE_ENABLED_MASK  0 //t should be already defined in intrinh, but Rtootls stsill compains

        void     cpuid(int32_t out[4], int32_t eax, int32_t ecx) {   __cpuidex(out, eax, ecx);      }
        uint64_t xgetbv(unsigned int x)                          {   return _xgetbv(x);             }
 
        //  Detect 64-bit - Note that this snippet of code for detecting 64-bit has been copied from MSDN.
        typedef BOOL(WINAPI* LPFN_ISWOW64PROCESS) (HANDLE, PBOOL);
        BOOL IsWow64(void)           {

            BOOL                bIsWow64         = FALSE;
            LPFN_ISWOW64PROCESS fnIsWow64Process = (LPFN_ISWOW64PROCESS)GetProcAddress(
                GetModuleHandle(TEXT("kernel32")), "IsWow64Process");

            if (NULL != fnIsWow64Process)       {
                if (!fnIsWow64Process(GetCurrentProcess(), &bIsWow64))     {
                    r_printf("Error Detecting Operating System.\n");
                    r_printf("Defaulting to 32-bit OS.\n\n");
                    bIsWow64 = FALSE;
                }
            }
            return bIsWow64;
        }

        int detect_OS_x64(void) {
            #ifdef _M_X64
                    return 1;
            #else
                    return IsWow64() != 0;
            #endif
        }
 
    /*********************************************/
    //        LINUX
    /*********************************************/
    #elif defined(__GNUC__) || defined(__clang__)

        #include <cpuid.h>
        void cpuid(int32_t out[4], int32_t eax, int32_t ecx) {
            __cpuid_count(eax, ecx, out[0], out[1], out[2], out[3]);
        }
        uint64_t xgetbv(unsigned int index) {
            uint32_t eax, edx;
            __asm__ __volatile__("xgetbv" : "=a"(eax), "=d"(edx) : "c"(index));
            return ((uint64_t)edx << 32) | eax;
        }
        #define _XCR_XFEATURE_ENABLED_MASK  0
        int detect_OS_x64(void) {
            //  We only support x64 on Linux.
            return 1;
        }

    #elif defined (OS_SOLARIS)
    
    // cpuid:  https:// docs.oracle.com/cd/E23824_01/html/821-1475/cpuid-7d.html

    // __asm and __asm__ are synonoous for Solaris CC
    // https:// docs.oracle.com/cd/E26502_01/html/E28387/gmacx.html#scrolltoc
    // Oracle® Solaris Studio 12.4: C User'sGuide

   //
        #define __cpuid(level, a, b, c, d)                                             \
          __asm__("cpuid\n\t" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "0"(level))

        #define __cpuid_count(level, count, a, b, c, d)                                \
          __asm__("cpuid\n\t"                                                          \
                  : "=a"(a), "=b"(b), "=c"(c), "=d"(d)                                 \
                  : "0"(level), "2"(count))

        void cpuid(int32_t out[4], int32_t eax, int32_t ecx) {
            __cpuid_count(eax, ecx, out[0], out[1], out[2], out[3]);
        }
        #define _XCR_XFEATURE_ENABLED_MASK  0
        uint64_t xgetbv(unsigned int index) {
            uint32_t eax, edx;
            __asm__ ("xgetbv" : "=a"(eax), "=d"(edx) : "c"(index));
            return ((uint64_t)edx << 32) | eax;
        }

        int detect_OS_x64(void) {
            //  We only support x64 on Linux.
            return 1;
        }
        #warning "Warning:No cpuid intrinsic defined for compiler."
    #else
        #define _XCR_XFEATURE_ENABLED_MASK  0
        void cpuid(int32_t out[4], int32_t eax, int32_t ecx) {
            out[0]=out[1]=out[2]= out[3]=0;
        }
 
        uint64_t xgetbv(unsigned int index) {            
            return 0;
        }
        int detect_OS_x64(void) {
            //  We only support x64 on Linux.
            return 1;
        }
        //#error  "No cpuid intrinsic defined for compiler."
        #warning "No cpuid intrinsic defined for compiler: a placeholder created!"
    #endif
#elif defined(cpu_ARM64) || defined (cpu_POWERPC64)
        //https://stackoverflow.com/questions/60588765/how-to-get-cpu-brand-information-in-arm64
        //https://stackoverflow.com/questions/23934862/what-predefined-macro-can-i-use-to-detect-the-target-architecture-in-clang
        #define _XCR_XFEATURE_ENABLED_MASK  0
        void cpuid(int32_t out[4], int32_t eax, int32_t ecx) {
            out[0] = out[1] = out[2] = out[3] = 0;
        }

        uint64_t xgetbv(unsigned int index) {
            return 0;
        }
        int detect_OS_x64(void) {
            //  We only support x64 on Linux.
            return 1;
        }
       // #warning "No cpuid intrinsic defined for ARM64: a placeholder created!"

// https://stackoverflow.com/questions/34488604/equivalent-for-sse-in-power-pc
// There is no SSE/AVX for POWEPC,Altivec is the POWERPC alternatives
#else
#   error "No cpuid intrinsic defined for processor architecture."
#endif

 
 
static void cpu_print(const char* label, uint8_t yes) {
    r_printf("%s%s\n", label, (yes ? "Yes" : "No"));
}
 
 
uint8_t detect_OS_AVX(void){
    //  Copied from: http://stackoverflow.com/a/22521619/922184

    Bool avxSupported = _False_;

    int  cpuInfo[4];
    cpuid(cpuInfo, 1, 0);

    Bool osUsesXSAVE_XRSTORE = (cpuInfo[2] & (1 << 27)) != 0;
    Bool cpuAVXSuport        = (cpuInfo[2] & (1 << 28)) != 0;

    if (osUsesXSAVE_XRSTORE && cpuAVXSuport)   {
        uint64_t xcrFeatureMask = xgetbv(_XCR_XFEATURE_ENABLED_MASK);
        avxSupported = (xcrFeatureMask & 0x6) == 0x6;
    }

    return avxSupported;
}
Bool detect_OS_AVX512(void){
    if (!detect_OS_AVX())
        return _False_;

    uint64_t xcrFeatureMask = xgetbv(_XCR_XFEATURE_ENABLED_MASK);
    return (xcrFeatureMask & 0xe6) == 0xe6;
}
void get_vendor_string(char * name){
    int32_t CPUInfo[4];
    //char name[13];

    cpuid(CPUInfo, 0, 0);
    memcpy(name + 0, &CPUInfo[1], 4);
    memcpy(name + 4, &CPUInfo[3], 4);
    memcpy(name + 8, &CPUInfo[2], 4);
    name[12] = '\0';
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void detect_host(struct cpu_x86 *cpu){

    memset(cpu, 0, sizeof(struct cpu_x86));

    //  OS Features
    cpu->OS_x64    = detect_OS_x64();
    cpu->OS_AVX    = detect_OS_AVX();
    cpu->OS_AVX512 = detect_OS_AVX512();

    //  Vendor
    char vendor[13];
    get_vendor_string(vendor);
    if (strcmp(vendor,"GenuineIntel") ==0){
        cpu->Vendor_Intel = _True_;
    }else if (strcmp(vendor, "AuthenticAMD") ==0){
        cpu->Vendor_AMD = _True_;
    }

    int info[4];
    cpuid(info, 0, 0);
    int nIds = info[0];

    cpuid(info, 0x80000000, 0);
    uint32_t nExIds = info[0];

    // r_printf("%d\n",nIds);

    //  Detect Features
    if (nIds >= 0x00000001){
        cpuid(info, 0x00000001, 0);
        cpu->HW_MMX    = (info[3] & ((uint32_t)1 << 23)) != 0;
        cpu->HW_SSE    = (info[3] & ((uint32_t)1 << 25)) != 0;
        cpu->HW_SSE2   = (info[3] & ((uint32_t)1 << 26)) != 0;
        cpu->HW_SSE3   = (info[2] & ((uint32_t)1 <<  0)) != 0;

        cpu->HW_SSSE3  = (info[2] & ((uint32_t)1 <<  9)) != 0;
        cpu->HW_SSE41  = (info[2] & ((uint32_t)1 << 19)) != 0;
        cpu->HW_SSE42  = (info[2] & ((uint32_t)1 << 20)) != 0;
        cpu->HW_AES    = (info[2] & ((uint32_t)1 << 25)) != 0;

        cpu->HW_AVX    = (info[2] & ((uint32_t)1 << 28)) != 0;
        cpu->HW_FMA3   = (info[2] & ((uint32_t)1 << 12)) != 0;

        cpu->HW_RDRAND = (info[2] & ((uint32_t)1 << 30)) != 0;
    }
    if (nIds >= 0x00000007){
        cpuid(info, 0x00000007, 0);
        cpu->HW_AVX2         = (info[1] & ((int)1 <<  5)) != 0;

        cpu->HW_BMI1         = (info[1] & ((int)1 <<  3)) != 0;
        cpu->HW_BMI2         = (info[1] & ((int)1 <<  8)) != 0;
        cpu->HW_ADX          = (info[1] & ((int)1 << 19)) != 0;
        cpu->HW_MPX          = (info[1] & ((int)1 << 14)) != 0;
        cpu->HW_SHA          = (info[1] & ((int)1 << 29)) != 0;
        cpu->HW_RDSEED       = (info[1] & ((int)1 << 18)) != 0;
        cpu->HW_PREFETCHWT1  = (info[2] & ((int)1 <<  0)) != 0;
        cpu->HW_RDPID        = (info[2] & ((int)1 << 22)) != 0;

        cpu->HW_AVX512_F     = (info[1] & ((int)1 << 16)) != 0;
        cpu->HW_AVX512_CD    = (info[1] & ((int)1 << 28)) != 0;
        cpu->HW_AVX512_PF    = (info[1] & ((int)1 << 26)) != 0;
        cpu->HW_AVX512_ER    = (info[1] & ((int)1 << 27)) != 0;
        //https://stackoverflow.com/questions/53566029/1-31-cannot-be-represented-by-type-int
        //Why does -fsanitize=undefined throw runtime error: left shift of 1 by 31 places cannot be represented in type 'int'
        cpu->HW_AVX512_VL    = (info[1] & ((uint32_t)1 << 31)) != 0;
        cpu->HW_AVX512_BW    = (info[1] & ((int)1 << 30)) != 0;
        cpu->HW_AVX512_DQ    = (info[1] & ((int)1 << 17)) != 0;

        cpu->HW_AVX512_IFMA  = (info[1] & ((int)1 << 21)) != 0;
        cpu->HW_AVX512_VBMI  = (info[2] & ((int)1 <<  1)) != 0;

        cpu->HW_AVX512_VPOPCNTDQ = (info[2] & ((int)1 << 14)) != 0;
        cpu->HW_AVX512_4FMAPS    = (info[3] & ((int)1 <<  2)) != 0;
        cpu->HW_AVX512_4VNNIW    = (info[3] & ((int)1 <<  3)) != 0;

        cpu->HW_AVX512_VNNI      = (info[2] & ((int)1 << 11)) != 0;

        cpu->HW_AVX512_VBMI2     = (info[2] & ((int)1 <<  6)) != 0;
        cpu->HW_GFNI             = (info[2] & ((int)1 <<  8)) != 0;
        cpu->HW_VAES             = (info[2] & ((int)1 <<  9)) != 0;
        cpu->HW_AVX512_VPCLMUL   = (info[2] & ((int)1 << 10)) != 0;
        cpu->HW_AVX512_BITALG    = (info[2] & ((int)1 << 12)) != 0;


        cpuid(info, 0x00000007, 1);
        cpu->HW_AVX512_BF16      = (info[0] & ((int)1 <<  5)) != 0;
    }
    if (nExIds >= 0x80000001){
        cpuid(info, 0x80000001, 0);
        cpu->HW_x64   = (info[3] & ((int)1 << 29)) != 0;
        cpu->HW_ABM   = (info[2] & ((int)1 <<  5)) != 0;
        cpu->HW_SSE4a = (info[2] & ((int)1 <<  6)) != 0;
        cpu->HW_FMA4  = (info[2] & ((int)1 << 16)) != 0;
        cpu->HW_XOP   = (info[2] & ((int)1 << 11)) != 0;
    }
}

void print_cpuinfo(struct cpu_x86 *cpu) {
    r_printf("CPU Vendor:\n");
    cpu_print("    AMD         = ", cpu->Vendor_AMD);
    cpu_print("    Intel       = ", cpu->Vendor_Intel);
    r_printf(" \n");
 
    r_printf("OS Features:\n");
#ifdef _WIN32
    cpu_print("    64-bit      = ", cpu->OS_x64);
#endif
    cpu_print("    OS AVX      = ", cpu->OS_AVX);
    cpu_print("    OS AVX512   = ", cpu->OS_AVX512);
    

    r_printf("\nHardware Features:\n");
    cpu_print("    MMX         = ", cpu->HW_MMX);
    cpu_print("    x64         = ", cpu->HW_x64);
    cpu_print("    ABM         = ", cpu->HW_ABM);
    cpu_print("    RDRAND      = ", cpu->HW_RDRAND);
    cpu_print("    RDSEED      = ", cpu->HW_RDSEED);
    cpu_print("    BMI1        = ", cpu->HW_BMI1);
    cpu_print("    BMI2        = ", cpu->HW_BMI2);
    cpu_print("    ADX         = ", cpu->HW_ADX);
    cpu_print("    MPX         = ", cpu->HW_MPX);
    cpu_print("    PREFETCHW   = ", cpu->HW_PREFETCHW);
    cpu_print("    PREFETCHWT1 = ", cpu->HW_PREFETCHWT1);
    cpu_print("    RDPID       = ", cpu->HW_RDPID);
    cpu_print("    GFNI        = ", cpu->HW_GFNI);
    cpu_print("    VAES        = ", cpu->HW_VAES);
 
    r_printf("\nSIMD: 128-bit\n");
    cpu_print("    SSE         = ", cpu->HW_SSE);
    cpu_print("    SSE2        = ", cpu->HW_SSE2);
    cpu_print("    SSE3        = ", cpu->HW_SSE3);
    cpu_print("    SSSE3       = ", cpu->HW_SSSE3);
    cpu_print("    SSE4a       = ", cpu->HW_SSE4a);
    cpu_print("    SSE4.1      = ", cpu->HW_SSE41);
    cpu_print("    SSE4.2      = ", cpu->HW_SSE42);
    cpu_print("    AES-NI      = ", cpu->HW_AES);
    cpu_print("    SHA         = ", cpu->HW_SHA);
 

    r_printf("\nSIMD: 256-bit\n");
    cpu_print("    AVX         = ", cpu->HW_AVX);
    cpu_print("    XOP         = ", cpu->HW_XOP);
    cpu_print("    FMA3        = ", cpu->HW_FMA3);
    cpu_print("    FMA4        = ", cpu->HW_FMA4);
    cpu_print("    AVX2        = ", cpu->HW_AVX2);
 

    r_printf("\nSIMD: 512-bit\n");
    cpu_print("    AVX512-F         = ", cpu->HW_AVX512_F);
    cpu_print("    AVX512-CD        = ", cpu->HW_AVX512_CD);
    cpu_print("    AVX512-PF        = ", cpu->HW_AVX512_PF);
    cpu_print("    AVX512-ER        = ", cpu->HW_AVX512_ER);
    cpu_print("    AVX512-VL        = ", cpu->HW_AVX512_VL);
    cpu_print("    AVX512-BW        = ", cpu->HW_AVX512_BW);
    cpu_print("    AVX512-DQ        = ", cpu->HW_AVX512_DQ);
    cpu_print("    AVX512-IFMA      = ", cpu->HW_AVX512_IFMA);
    cpu_print("    AVX512-VBMI      = ", cpu->HW_AVX512_VBMI);
    cpu_print("    AVX512-VPOPCNTDQ = ", cpu->HW_AVX512_VPOPCNTDQ);
    cpu_print("    AVX512-4FMAPS    = ", cpu->HW_AVX512_4FMAPS);
    cpu_print("    AVX512-4VNNIW    = ", cpu->HW_AVX512_4VNNIW);
    cpu_print("    AVX512-VBMI2     = ", cpu->HW_AVX512_VBMI2);
    cpu_print("    AVX512-VPCLMUL   = ", cpu->HW_AVX512_VPCLMUL);
    cpu_print("    AVX512-VNNI      = ", cpu->HW_AVX512_VNNI);
    cpu_print("    AVX512-BITALG    = ", cpu->HW_AVX512_BITALG);
    cpu_print("    AVX512-BF16      = ", cpu->HW_AVX512_BF16);
 

    r_printf("\nSummary\n");
    cpu_print("    Safe to use AVX:     ", cpu->HW_AVX && cpu->OS_AVX);
    cpu_print("    Safe to use AVX512:  ", cpu->HW_AVX512_F && cpu->OS_AVX512);
    r_printf("\n");
}

void detect_print_cpu(void) {
    struct cpu_x86 cpuinfo;
     detect_host(&cpuinfo);
     print_cpuinfo(&cpuinfo);
}

 
//https://stackoverflow.com/questions/12594208/c-program-to-determine-levels-size-of-cache
/*
To get L2 and L3 cache sizes, you need to use CPUID with EAX=4, and set ECX to 0, 1, 2, ... for
each caching level. The linked code shows this, and Intel's docs have details on which bits mean 
what.
*/

/*
When eax = 4 and ecx is the cache level,
Ways = EBX[31:22]
Partitions = EBX[21:12]
LineSize = EBX[11:0]
Sets = ECX
Total Size = (Ways + 1) * (Partitions + 1) * (Line_Size + 1) * (Sets + 1)
*/
void i386_cpuid_caches (Bool quiet) {
 
    for (int i = 0; i < 32; i++) {

        // Variables to hold the contents of the 4 i386 legacy registers
        uint32_t eax, ebx, ecx, edx; 

        eax = 4; // get cache info
        ecx = i; // cache id

        #if !defined(COMPILER_MSVC) && !defined(cpu_ARM64) && !defined (cpu_POWERPC64)
            __asm__ (
                "cpuid" // call i386 cpuid instruction
                : "+a" (eax) // contains the cpuid command code, 4 for cache query
                , "=b" (ebx)
                , "+c" (ecx) // contains the cache id
                , "=d" (edx)
            ); // generates output in 4 registers eax, ebx, ecx and edx 
        
        #else

            int32_t out[4];
            cpuid(out, eax, ecx);
            eax = out[0];
            ebx = out[1];
            ecx = out[2];
            edx = out[3];

        #endif

        // See the page 3-191 of the manual.
        int cache_type = eax & 0x1F; 

        if (cache_type == 0) // CODE_EOF of valid cache identifiers
            break;

        char * cache_type_string;
        switch (cache_type) {
            case 1: cache_type_string = "Data Cache"; break;
            case 2: cache_type_string = "Instruction Cache"; break;
            case 3: cache_type_string = "Unified Cache"; break;
            default: cache_type_string = "Unknown Type Cache"; break;
        }

        int cache_level = (eax >>= 5) & 0x7;

        int cache_is_self_initializing = (eax >>= 3) & 0x1; // does not need SW initialization
        int cache_is_fully_associative = (eax >>= 1) & 0x1;

        // See the page 3-192 of the manual.
        // ebx contains 3 integers of 10, 10 and 12 bits respectively
        unsigned int cache_sets = ecx + 1;
        unsigned int cache_coherency_line_size = (ebx & 0xFFF) + 1;
        unsigned int cache_physical_line_partitions = ((ebx >>= 12) & 0x3FF) + 1;
        unsigned int cache_ways_of_associativity = ((ebx >>= 10) & 0x3FF) + 1;

        // Total cache size is the product
        size_t cache_total_size = cache_ways_of_associativity * cache_physical_line_partitions * cache_coherency_line_size * cache_sets;

        if (!quiet)
            r_printf(
                "Cache ID %d:\n"
                "- Level: %d\n"
                "- Type: %s\n"
                "- Sets: %d\n"
                "- System Coherency Line Size: %d bytes\n"
                "- Physical Line partitions: %d\n"
                "- Ways of associativity: %d\n"
                "- Total Size: %zu bytes (%zu kb)\n"
                "- Is fully associative: %s\n"
                "- Is Self Initializing: %s\n"
                "\n"
                , i
                , cache_level
                , cache_type_string
                , cache_sets
                , cache_coherency_line_size
                , cache_physical_line_partitions
                , cache_ways_of_associativity
                , cache_total_size, cache_total_size >> 10
                , cache_is_fully_associative ? "true" : "false"
                , cache_is_self_initializing ? "true" : "false"
            );
    }
}

#include "abc_000_warning.h"