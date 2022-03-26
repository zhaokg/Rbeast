#pragma once

#include "abc_datatype.h" 
struct cpu_x86 {
    //  Vendor
    bool Vendor_AMD;
    bool Vendor_Intel;

    //  OS Features
    bool OS_x64;
    bool OS_AVX;
    bool OS_AVX512;

    //  Misc.
    bool HW_MMX;
    bool HW_x64;
    bool HW_ABM;
    bool HW_RDRAND;
    bool HW_RDSEED;
    bool HW_BMI1;
    bool HW_BMI2;
    bool HW_ADX;
    bool HW_MPX;
    bool HW_PREFETCHW;
    bool HW_PREFETCHWT1;
    bool HW_RDPID;

    //  SIMD: 128-bit
    bool HW_SSE;
    bool HW_SSE2;
    bool HW_SSE3;
    bool HW_SSSE3;
    bool HW_SSE41;
    bool HW_SSE42;
    bool HW_SSE4a;
    bool HW_AES;
    bool HW_SHA;

    //  SIMD: 256-bit
    bool HW_AVX;
    bool HW_XOP;
    bool HW_FMA3;
    bool HW_FMA4;
    bool HW_AVX2;

    //  SIMD: 512-bit
    bool HW_AVX512_F;
    bool HW_AVX512_CD;

    //  Knights Landing
    bool HW_AVX512_PF;
    bool HW_AVX512_ER;

    //  Skylake Purley
    bool HW_AVX512_VL;
    bool HW_AVX512_BW;
    bool HW_AVX512_DQ;

    //  Cannon Lake
    bool HW_AVX512_IFMA;
    bool HW_AVX512_VBMI;

    //  Knights Mill
    bool HW_AVX512_VPOPCNTDQ;
    bool HW_AVX512_4FMAPS;
    bool HW_AVX512_4VNNIW;

    //  Cascade Lake
    bool HW_AVX512_VNNI;

    //  Cooper Lake
    bool HW_AVX512_BF16;

    //  Ice Lake
    bool HW_AVX512_VBMI2;
    bool HW_GFNI;
    bool HW_VAES;
    bool HW_AVX512_VPCLMUL;
    bool HW_AVX512_BITALG;
};

extern void cpu_print_info();
extern void detect_host(struct cpu_x86* cpu);
extern void i386_cpuid_caches();

 