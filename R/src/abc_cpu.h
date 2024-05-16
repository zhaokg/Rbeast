#pragma once

#include "abc_datatype.h" 
struct cpu_x86 {
    //  Vendor
    Bool Vendor_AMD;
    Bool Vendor_Intel;

    //  OS Features
    Bool OS_x64;
    Bool OS_AVX;
    Bool OS_AVX512;

    //  Misc.
    Bool HW_MMX;
    Bool HW_x64;
    Bool HW_ABM;
    Bool HW_RDRAND;
    Bool HW_RDSEED;
    Bool HW_BMI1;
    Bool HW_BMI2;
    Bool HW_ADX;
    Bool HW_MPX;
    Bool HW_PREFETCHW;
    Bool HW_PREFETCHWT1;
    Bool HW_RDPID;

    //  SIMD: 128-bit
    Bool HW_SSE;
    Bool HW_SSE2;
    Bool HW_SSE3;
    Bool HW_SSSE3;
    Bool HW_SSE41;
    Bool HW_SSE42;
    Bool HW_SSE4a;
    Bool HW_AES;
    Bool HW_SHA;

    //  SIMD: 256-bit
    Bool HW_AVX;
    Bool HW_XOP;
    Bool HW_FMA3;
    Bool HW_FMA4;
    Bool HW_AVX2;

    // SIMD: 512-bit
    Bool HW_AVX512_F;
    Bool HW_AVX512_CD;

    //  Knights Landing
    Bool HW_AVX512_PF;
    Bool HW_AVX512_ER;

    //  Skylake Purley
    Bool HW_AVX512_VL;
    Bool HW_AVX512_BW;
    Bool HW_AVX512_DQ;

    //  Cannon Lake
    Bool HW_AVX512_IFMA;
    Bool HW_AVX512_VBMI;

    //  Knights Mill
    Bool HW_AVX512_VPOPCNTDQ;
    Bool HW_AVX512_4FMAPS;
    Bool HW_AVX512_4VNNIW;

    //  Cascade Lake
    Bool HW_AVX512_VNNI;

    //  Cooper Lake
    Bool HW_AVX512_BF16;

    //  Ice Lake
    Bool HW_AVX512_VBMI2;
    Bool HW_GFNI;
    Bool HW_VAES;
    Bool HW_AVX512_VPCLMUL;
    Bool HW_AVX512_BITALG;
};

extern void detect_print_cpu(void);
extern void detect_host(struct cpu_x86* cpu);
extern void i386_cpuid_caches(Bool quiet);

extern void print_cpuinfo(struct cpu_x86* cpu);