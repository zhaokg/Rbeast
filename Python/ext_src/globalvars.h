#pragma once
#include "abc_datatype.h"
#include "abc_pthread.h"

extern VOID_PTR  GLOBAL_OPTIONS;
extern VOID_PTR  GLOBAL_RESULT;

extern pthread_t       *thread_id;
extern pthread_mutex_t mutex;
extern pthread_cond_t  condVar;


extern volatile F32    PERCENT_COMPLETED;
extern volatile F32    REMAINING_TIME;
extern volatile I32    NEXT_PIXEL_INDEX;
extern volatile I32    NUM_OF_PROCESSED_PIXELS;
extern volatile I32    NUM_OF_PROCESSED_GOOD_PIXELS;



extern int  IS_CPU_INSTRUCTON_SET;