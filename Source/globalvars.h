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

extern char GLOBAL_PRNT_WARNING;
extern char GLOBAL_PRNT_CPU;
extern char GLOBAL_PRNT_PARAMETER;
extern char GLOBAL_PRNT_PROGRESS;
extern char GLOBAL_IS_QUIET_MODE;

extern char GLOBAL_CPU_REQUEST;
extern  char GLOBAL_CPU_CURRENT;

//extern float   GLOBAL_MODEL_PRIOR_FACTOR ;

extern int  beast2_main_corev4(void);
extern void beast2_main_corev4_gui(void);
extern int  beast2_main_corev4_mthrd(void* dummy);
