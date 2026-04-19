//Must be include twice in a unit. So #pragram once must be removed

#include "abc_000_macro.h"
//#pragma GCC diagnostic ignored "-Wunused-variable"

#ifndef WARNING_SWITCH
	#define WARNING_SWITCH
    DISABLE_MANY_WARNINGS
#else
    // For some reason,, on  MacOS, sign and warning defined in R. h expand into Rf_sign 
    // to mess up the pragma ingore statements. A sloppy solution is to simply undefine sign and wanring

    #ifdef sign
      #undef sign
    #endif

    #ifdef warning
      #undef warning 
    #endif

	ENABLE_MANY_WARNINGS 	
#endif
