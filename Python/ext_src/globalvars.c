#include "abc_000_macro.h"
#include "abc_000_warning.h"
//Added before globalvar.h bcz abc_pthread is included, which has a few static functions not already 
// used. We want to force -Wno-unused-function
#include "globalvars.h"

VOID_PTR  GLOBAL_OPTIONS;
VOID_PTR  GLOBAL_RESULT;
VOID_PTR  GLOBAL_DATA0;

pthread_t       *thread_id;
pthread_mutex_t  mutex;
pthread_cond_t   condVar;


volatile F32    PERCENT_COMPLETED;
volatile F32    REMAINING_TIME;
volatile I32    NEXT_PIXEL_INDEX;
volatile I32    NUM_OF_PROCESSED_PIXELS;
volatile I32    NUM_OF_PROCESSED_GOOD_PIXELS;
 


char GLOBAL_PRNT_WARNING ;
char GLOBAL_PRNT_CPU       ;
char GLOBAL_PRNT_PARAMETER;
char GLOBAL_PRNT_PROGRESS  ;
char GLOBAL_IS_QUIET_MODE  ;

char GLOBAL_CPU_REQUEST;
char GLOBAL_CPU_CURRENT = 0;


//float   GLOBAL_MODEL_PRIOR_FACTOR = 1.0;

#include "abc_000_warning.h"