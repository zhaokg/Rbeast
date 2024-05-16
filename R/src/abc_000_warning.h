//Must be include twice in a unit. So #pragram once must be removed

#include "abc_000_macro.h"
//#pragma GCC diagnostic ignored "-Wunused-variable"

#ifndef WARNING_SWITCH
	#define WARNING_SWITCH
    DISABLE_MANY_WARNINGS
#else
    // for MacOS, for some reason, sing and warning defined in R. h will expand into Rf_sign 
    // to mess up the pragram ingore statements.
    #undef sign    
    #undef warning 
	ENABLE_MANY_WARNINGS 	
#endif
