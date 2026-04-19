#pragma once

#include "abc_000_macro.h"
#include "abc_datatype.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void Sleep_ms(int milliseconds);

// Warning: function declaration isn't a prototype [-Wstrict-prototypes]
// stackoverflow.com/questions/42125/warning-error-function-declaration-isnt-a-prototype 

extern void Timer_Start(void);
extern F64  Timer_ElapsedSecond(void);
extern void Timer_SetBreakPt(void);
extern F64  Timer_ElapsedSinceBeakPt(void);
extern U64  Timer_TickCount(void);

extern unsigned long long tic(void);
extern unsigned long long toc(void);


#ifdef __cplusplus
}
#endif