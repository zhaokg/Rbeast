#pragma once
#include "abc_datatype.h"
#include "abc_pthread.h"

extern VOID_PTR  GLOBAL_OPTIONS;
extern VOID_PTR  GLOBAL_RESULT;
extern VOID_PTR  GLOBAL_DATA0;

extern pthread_t       *thread_id;
extern pthread_mutex_t mutex;
extern pthread_cond_t  condVar;


extern volatile F32    PERCENT_COMPLETED;
extern volatile F32    REMAINING_TIME;
extern volatile I32    NEXT_PIXEL_INDEX;
extern volatile I32    NUM_OF_PROCESSED_PIXELS;
extern volatile I32    NUM_OF_PROCESSED_GOOD_PIXELS;

extern char GLOBAL_QUIET_MODE;


extern int beast2_main_corev4(void);
extern void beast2_main_corev4_gui(void);
extern int beast2_main_corev4_mthrd(void* dummy);
